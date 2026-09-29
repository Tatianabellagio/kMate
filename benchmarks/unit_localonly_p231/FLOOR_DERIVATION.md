# Is `--min-kmers-per-block 200` arbitrary? — derivation + sweep

**Question.** The local-only / anchored window mode gates each block on a minimum
nonzero-k-mer count (`--min-kmers-per-block`, default **200**, **50** in the
bench): below it a block falls back to global-h (anchored) or abstains (local-only).
Is the *variable* (a k-mer count) and the *value* (200) principled, or arbitrary?

**Method** (`scripts/block_floor_diag.py` + `analyze_block_floor.py` + `plot_floor_curve.py`).
3 p231 sims spanning the founder-count and regime axes — `n50_g0_out`,
`n231_g1_self97`, `n50_g3_dom500nr` — fit local-only with **floor=1, no anchor,
no HMM smoothing** (every non-empty block fit on its own). One run gives the whole
sweep in post (a record is "called at floor T" iff its block has nnz ≥ T). Per
block we logged the raw count `nnz` plus principled resolvability metrics:
`effrank_design` (spectral-entropy rank of the founder-carriage matrix),
`effrank_fisher`/`cond` (the codebase's own `fisher_information_h` +
`_resolvability_from_J`), and `n_present`.

## Finding 1 — the *variable* is right: k-mer count predicts AF error best

Spearman ρ of per-block AF RMSE vs each predictor (negative = predicts lower error):

| pool / unit | **nnz** | effrank_design | effrank_fisher | n_present |
|---|---|---|---|---|
| n231 self · w10kb | **−0.52** | +0.29 | −0.18 | −0.42 |
| n231 self · dynld | **−0.61** | +0.09 | −0.42 | −0.38 |
| n50 g0 · w10kb | **−0.51** | +0.31 | −0.16 | −0.44 |
| n50 g3 dom · w10kb | **−0.46** | +0.14 | −0.17 | −0.33 |

`nnz` is the strongest, most consistent predictor — beating the identifiability
metrics. **Why, against intuition:** the scored estimand is per-record *AF*, not
per-founder *h*. AF cancels founder-collinearity (carriers in a block share their
`var_pa` value, so the non-identifiability bias cancels — the h-certainty result).
So block AF error is governed by **counting noise / total Fisher information ∝
nnz**, not by whether individual founders are separable. Gating on a k-mer count
is therefore the theoretically correct choice of variable. (The entropy-rank
`effrank_design` even flips positive: conditional on supply, higher local
ancestry diversity = harder = more error — it measures mixture complexity, not
precision.)

### Finding 1a — the low-nnz blocks are the centromere (not noise)

`block_floor_strip.png` (raw-data/swarm view, centromere blocks in red). The non-monotonic
bump in the low-nnz w10kb bins (e.g. 16–31) is the **Chr1 (peri)centromere**: those bins
are 50–75% centromeric (peaking at 16–31), and the worst windows sit at ~14.4–15.3 Mb (CEN1)
with high true AF. A 10 kb window with only 16–31 nonzero k-mers is k-mer-starved despite
normal size — repeat-rich sequence where panel k-mers are non-unique / repeat-guarded away.
So the floor is gating exactly the right blocks (pathological centromere/repeat windows).

### Finding 1b — it's k-mer COUNT, not DENSITY (why dynld ≈ 10 kb)

`block_floor_normalized.png` + `block_density_corr.tsv`. The error-vs-nnz curves collapse
dynld onto w10kb. Is total supply the sufficient statistic, or would a *density* (k-mers
per variant / per kb) explain more? Per-block *h* is fit **once** from all the block's
k-mers and every variant projects through that same *h*, so accuracy should track total
k-mers, not packing. Spearman ρ(RMSE, ·) on **dynld** (variable block size):

| nnz | k-mers/variant | k-mers/kb | n_variants |
|---|---|---|---|
| **−0.61** | −0.27 | −0.13 | −0.36 |

Raw count wins; density is the *worst* predictor. Within a fixed nnz band, density does
**not** reduce error (ρ≈0 or wrong-signed). So **at equal total k-mers a wide variant-rich
block and a tiny one have the same AF error** — which is exactly why the two units lie on
one curve vs nnz. Only the raw-count panel of `block_floor_normalized.png` collapses.

## Finding 2 — the *value* 200 is NOT where error plateaus

`block_floor_curve.png` — pooled median + p90 per-block AF RMSE vs nnz:

| nnz bin | median RMSE (w10kb) | p90 RMSE |
|---|---|---|
| [128, 256) | 0.102 | 0.191 |
| [256, 512) | 0.088 | 0.182 |
| [512, 1024) | 0.072 | 0.136 |
| [1024, 2048) | 0.052 | 0.086 |
| [2048, ∞) | 0.048 | 0.077 |

(distribution view: `block_floor_box.png` — boxplots of per-block RMSE per nnz bin, box=IQR,
whiskers=5–95 pct; the whole distribution, not just the median, collapses to the floor only
at nnz ≳ 1k.) Median block error keeps falling until **nnz ≈ 1000–2000**, where it hits the
accuracy floor (~0.048). At the current floors median error is ~2× the floor
(≈0.10 at nnz=200, ≈0.10–0.11 at 50) and the p90 tail stays bad (>0.18) until
nnz ≥ 1024. So **200 is not "enough to resolve the mixture" — it is a pragmatic
"minimum worth attempting."** This sits below, but is the small-sibling of, the
~7000-k-mer "full ecotype resolution" figure: AF needs less than full h
resolution, but still ~1–2k, far above 200.

## Finding 3 — there is no clean elbow; it's a coverage/accuracy dial

The sweep (`block_floor_sweep.png`) is smooth and concave: each doubling of the
floor buys ~+0.01–0.02 R² at a rising coverage cost. Marginal R²-per-coverage is
roughly flat up to ~200 then falls off above ~400. So 200 sits at a defensible
*soft knee*, and 50 is a valid choice that keeps ~6 pts more coverage for a small
accuracy cost — but neither is a discoverable optimum. To actually reach the
accuracy floor you'd need a floor (~1000) so high it abstains on ~70% of records
(windows just don't hold that many k-mers at 10× / 10 kb) — which is exactly why
the production answer is **global/anchored mode, not a higher local floor**.

## Finding 4 — panel & the private-k-mer filter (p80 vs p231)

`panel_compare_floor.png` (line), `panel_compare_strip.png` (strip distribution),
`panel_compare_joint.png` (joint density: hexbin nnz×RMSE by block count [reversed cmap so
the sparse tail shows] + nnz/RMSE marginals — makes the rightward supply shift of
*unfiltered* and the centromere RMSE tail legible; annotated with the **% of scored blocks
dropped** at min-kmers floors 100/200/500), job 35291298.

**Drop rate vs floor** (% of scored blocks with nnz below the floor → NaN): at **100**
p231 5.8% / p80-filt2 4.2% / p80-unfilt 5.0%; at **200** 8.6 / 7.1 / 7.1; at **500**
15.7 / 13.9 / 12.7. So even 200 costs <9% coverage, and *unfiltered* drops the fewest at
every floor (the supply shift again). **RMSE payoff is modest/diminishing**: the kept-block
mean falls 0.073→0.068→0.066→0.062 (p231, floor 0/100/200/500) while the *median* is nearly
flat 0.058→0.054 — the floor only removes the high-error low-nnz tail (felt by the mean, not
the median). p80-unfilt sits lowest throughout (μ 0.057→0.051).
The 231 panel uses **filt2inv** (drop `ac==1` private singletons — error/repeat-prone in a
*mixed* long+short-read panel — and `ac==F` invariants). Tested on the 80-founder all-long-
read panel in two variants. Centromere penalty at matched supply (nnz≥512):

| series | arm RMSE | cen RMSE | penalty |
|---|---|---|---|
| p231 filt2inv (231 mixed) | 0.053 | 0.080 | +0.027 |
| p80 filt2 (80 long-read) | 0.054 | 0.071 | +0.016 |
| **p80 unfiltered (keep private)** | **0.044** | **0.058** | **+0.014** |

1. **Keeping private k-mers helps** — and *at matched nnz* (a private ac=1 k-mer pins one
   founder, so it raises quality-per-k-mer, not just supply). Arm error drops to 0.044,
   below the filtered ~0.05 floor; it also shifts blocks to higher nnz (centromere median
   nnz 464→522, arm 1628→1899).
2. The long-read panel helps **specifically in the centromere** (penalty halved) but not the
   arm (coverage-limited either way). **Confound:** p80 also has fewer founders (80 vs 231).
3. The centromere penalty **persists even unfiltered** (+0.014) → bias (coverage/missingness/
   single-h), not uniqueness — consistent with Finding 1a.

> ⚠️ **Closed-loop caveat.** Sim reads are generated *from* the panel, so every panel private
> k-mer is real — no short-read assembly artifacts, the exact failure mode filt2inv targets.
> So this is the **best case** for private k-mers: they carry real signal *when the panel's
> private k-mers are trustworthy* (plausible for all-long-read p80; not for the mixed 231
> panel — hence the filter). Repeat guard was off (consistent across series).

## Takeaway for the paper

- Gating on a **k-mer count** is principled — it is the sufficient statistic for
  AF precision *because AF is collinearity-robust* (ties to the h-uncertainty
  framework). State this rather than leaving the floor unexplained.
- The value 200 is a coverage/accuracy operating point near a soft knee, **not**
  the resolution requirement. The resolution requirement for AF is ~1–2k k-mers/
  block (≪ the ~7k for full h), unreachable at 10× in 5–10 kb windows.
- This reinforces the standing decision: **global/anchored for production**;
  local-only only when a global prior would mask the signal, at a real cost.
