# Balancing k-mers across founders — brainstorm

**Created 2026-05-13.**

The +41% cactus h-bias in `cactus_em-v3` traces to per-founder k-mer-count
asymmetry in `cn_full_v3`. CV across the 231 founders' unique k-mer counts is
2.6% in v3 vs 0.9% in v2. The EM is a multiplicative-update simplex
algorithm — small per-row biases compound across hundreds of iterations into
the observed +41% h shift on cactus founders.

**Goal of this doc:** enumerate options for equalizing per-founder k-mer
contribution in `cn_full`, so the EM treats every founder's evidence with
equal weight regardless of how many unique k-mers their haplotype happens
to carry.

**Explicitly NOT considered:** Beagle imputation. Hard-rejected for SVs in
2026-05-03 LOO testing (−25 to −30pp concordance loss on small/medium SVs);
SVs are the variant class this project cares most about. See
`memory/feedback_no_beagle_solutions.md`.

## Why PanGenie's per-allele cap doesn't fix this

PanGenie indexes k-mers with three filters:

1. Internal var-dedup: k-mer unique within bubble's alleles
2. Internal ref-dedup: k-mer absent elsewhere in graph/reference
3. **Hard cap**: 16 k-mers per allele for biallelic, 32 per allele for multi-allelic bubbles

The cap is **per allele** (per bubble × per ALT), not per founder. The
asymmetry source is:

- Cactus founders have many private assembly variants → many bubbles where
  they carry distinguishing ALTs → many unique k-mers attributable to them.
- PG founders are locked into PanGenie short-read calls → fewer alleles
  uniquely theirs → fewer founder-discriminating k-mers per founder.

The cap controls **how many** k-mers each ALT gets; it does NOT control
**how that k-mer mass is distributed across founders.** That's where the
asymmetry lives, and that's what we need to balance.

## Five candidate approaches

### 1. Row-normalize cn_full at EM-load time ✗ TESTED 2026-05-13, FAILED

**Tested on `cov10_n200_g1` global mode (job 61985):**

| variant | R² | MAE | slope | cactus-specific mean residual | balanced | PG-specific |
|---|---:|---:|---:|---:|---:|---:|
| baseline (no row-norm) | 0.994 | 0.012 | +0.984 | +0.008 | +0.001 | −0.022 |
| row-norm + K_f correction | 0.955 | 0.035 | +1.039 | −0.019 | −0.006 | +0.040 |
| row-norm, NO correction  | 0.955 | 0.035 | +1.039 | −0.019 | −0.006 | +0.040 |

**Outcome: the bias direction FLIPS and grows.** Cactus residual went from
+0.008 to −0.019 (sign flip, 2.3× magnitude); PG residual went from −0.022 to
+0.040 (sign flip, 1.8× magnitude). Overall MAE is 3× worse. The K_f
post-EM correction barely matters because K_f only varies by 21% (max/min)
— the EM trajectory itself is what overshoots.

**Why it failed (post-hoc analysis):**

Row-normalizing replaces `cn[f,k] = 1` with `cn[f,k] = 1/K_f`. PG founders
have smaller K_f, so their per-k-mer rate coefficient `1/K_f` is LARGER. The
EM gives them more leverage per k-mer, pulling mass onto PG. Result: the EM
solution overshoots past the truth in the opposite direction.

In the original EM, the bias source is more subtle than "K_f imbalance" —
my prior derivation showed that for **founder-private k-mers** specifically,
the EM contribution is `c_k` regardless of cn value scaling. So row-norm
doesn't help with private-k-mer pinning at all; it only changes the
shared-k-mer redistribution math, in a direction that turns out to be
harmful.

**Decision:** the +41% h-bias is not predominantly K_f-driven. The dominant
mechanism is likely **founder-private-k-mer pinning** — assembly-specific
variants in cactus founders that each pin one founder's h above truth.
Approach #1 cannot address this. Skip to #2 or #3, or rethink the
mechanism with a dedicated diagnostic.



**Idea:** Replace `cn[f, k] = 1` (founder f carries k-mer k) with
`cn[f, k] = 1 / K_f`, where `K_f = Σ_k cn[f, k]` is founder f's total
unique-k-mer count.

**Effect:** Every founder's row sums to 1. Per-founder "evidence budget" is
equal by construction. Whichever founder is louder in raw counts has each of
their k-mers diluted to compensate.

**Math (M-step becomes):**
```
µ_k(h)      = Σ_f h_f · cn_norm[f, k]
h_new[f]    ∝ h[f] · Σ_k cn_norm[f, k] · c_k / µ_k(h)
```
λ cancels in the ratio as before; multiplicative-update simplex EM still
converges (Lee–Seung NMF generalizes to non-negative real matrices, not
just 0/1).

**Interpretation:** After convergence `h_f` is "founder f's mass × per-founder
k-mer-budget factor", not literal mass. To get literal mass for the AF
projection through `cn_var`, multiply `h_f` by `K_f` and renormalize the
result before computing `est_af = h_corrected · cn_var`.

**Pros:**
- Trivial implementation (~3 lines in `poolfreq/src/em_solver.py:load_cn`)
- Reuses existing `cn_full_v3` — no SLURM rebuild
- Keeps all signal — no k-mers dropped
- Deterministic (no random subsample)

**Cons:**
- `cn` matrix is no longer 0/1 — sparse CSR with int8 dtype changes to float
  (memory ~4× larger for dense format; sparse stays fine with float32)
- The h ↔ pool-mass mapping requires one extra step

**Effort:** ~1 day to implement + validate on n200_g1.

### 2. Per-founder k-mer subsample to min(K_f) at cn-build time

**Idea:** Find `K_min = min_f K_f`. For founders with `K_f > K_min`, randomly
subsample their k-mer rows down to `K_min`. Rebuild cn_full.

**Effect:** Each row literally has the same number of nonzero entries —
brute-force version of #1.

**Pros:**
- Matrix stays 0/1 (no interpretation change)
- Direct, physical equalization

**Cons:**
- Throws away k-mers (~`K_max - K_min` per cactus founder)
- Reproducibility depends on the random seed
- Requires a one-time cn-rebuild (~6-8h SLURM)

**Effort:** Modify `poolfreq/src/build_kmer_cn.py` to post-process the
output, then re-run cn_full build. Half-day code + overnight compute.

### 3. Per-bubble per-founder k-mer cap

**Idea:** PanGenie caps k-mers per ALLELE (16/allele biallelic, 32/allele
multi-allelic). Add a second cap: per founder per bubble (e.g. ≤ 4 k-mers
per founder per bubble). Operates at PanGenie's natural unit.

**Effect:** No single bubble can over-represent a single founder. Prevents
the case where a multi-allelic bubble's many cactus-private alleles all funnel
distinguishing k-mers into one cactus founder.

**Pros:**
- More targeted than #1/#2 — preserves rare bubbles where each founder gets
  few k-mers; only trims the cases of explosive per-founder accumulation
- Local fix, doesn't change global cn shape

**Cons:**
- Medium effort: need to walk through PanGenie-index output and identify
  which k-mers belong to which (bubble, founder) tuples before applying
  the cap
- Probably needs another full cn rebuild

**Effort:** 1-2 days code + overnight compute.

### 4. Drop singleton k-mers (k-mers carried by ≤ N founders) ✗ already tried

Filtered with N ∈ {2, 3, 5, 10} on cn_full_v3 in last session (filt2/3/5/10).
**None recovered v2 quality.** The very k-mers that distinguish individual
founders are also the ones that carry the bias signal. Removing them
removes both the bias *and* the identifiability that makes the EM work.

Not a viable path. Documented for completeness so we don't re-try it.

### 5. Intersection-only (k-mers with ≥1 cactus carrier AND ≥1 PG carrier)

**Idea:** Drop any k-mer that's exclusively cactus-carried or exclusively
PG-carried. Keep only the "shared" k-mers — those with at least one
carrier on each side of the panel split.

**Effect:** By construction, the group-specific signal that drives the
asymmetry is gone. Bias should disappear.

**Pros:** Conceptually clean — directly removes the asymmetry source.

**Cons:**
- Loses identifiability for cactus-specific variants entirely (no k-mer
  evidence pointing to those founders → EM can't recover their mass)
- Loses identifiability for PG-specific variants too
- Likely degrades R² on variants we currently get right

**Not recommended.**

## Implementation order

(Original ordering assumed #1 would be a clean win. Empirical test on
2026-05-13 showed #1 makes things worse — see the section above. Revised
ordering below.)

1. ~~**#1 (row-normalize)**~~ — TESTED, FAILED. Do not retry.

2. **Diagnose the actual mechanism more carefully before trying #2/#3.**
   Hypothesis: the bias is from founder-private k-mers, not K_f imbalance.
   Confirmation steps (no compute beyond a notebook):
   - Build a per-founder count of "k-mers carried ONLY by this founder"
     (true private k-mers, not just rare k-mers).
   - Plot cactus-vs-PG distribution of per-founder private-k-mer counts.
   - If cactus founders have many more private k-mers per founder, the
     pinning hypothesis is right.

3. If pinning is confirmed: **#2 (per-founder k-mer subsample to min K_f)**
   directly attacks it by reducing each cactus founder's k-mer count
   (including their private k-mers). Cost: throws away signal, but the
   signal it throws away is exactly the source of the bias.

4. Or **#3 (per-bubble per-founder cap)**: less aggressive than #2; targets
   the cases where a single bubble's multi-allelic structure assigns many
   k-mers to one founder.

5. If neither #2 nor #3 helps either, the bias is from a deeper source
   (simplex identifiability, recombination assumption violation, etc.)
   that k-mer-level rebalancing cannot fix.

## How to validate any of these

Plug the new cn into the existing eval pipeline:

1. Re-run `cactus_em` on the 4 cov10 regimes with the rebalanced cn_full.
2. Re-run `CARRIER_ASYMMETRY_DIAGNOSTIC.ipynb` — the binned-residual bars
   should flatten (no more monotone ramp from cactus-specific to PG-specific).
3. Stratified table in cell 5 should show ~zero mean residual in all three
   asymmetry classes.

The p82 control experiment (`HANDOFF_P82_CONTROL.md`) is the **prior** test:
it confirms the asymmetry IS the mechanism. If p82 streaks vanish AND a
rebalancing approach also kills the streaks on 231-founders, the mechanism
is closed.

## Open theoretical questions

- **Does the EM with row-normalized cn have the same identifiability as
  with binary cn?** Linear-algebraic answer: yes, identifiability is a
  property of the column space of cn, and row scaling preserves it (rank
  unchanged; null space transformed isomorphically). But it's worth
  confirming on a small synthetic before scaling up.

- **For #1, what's the right post-EM correction?** The projection
  `est_af = h · cn_var` should use the founder-mass interpretation of h.
  After row-norm EM, `h_f` is "mass × 1/K_f"; recovering literal mass
  requires multiplying by K_f and renormalizing. This is one extra step
  but well-defined.

- **What about the smoother and the anchor in ★★?** The Dirichlet anchor
  (`prior_h`) and the HMM smoother both operate on h directly; with #1
  they'd operate on the budget-factored h. Probably want to either
  (a) apply them before the budget correction, or (b) re-derive them in
  the rebalanced regime. Worth thinking through carefully.
