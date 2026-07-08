# EM non-identifiability, unit choice, and the "no prior" decision (2026-07-07)

Why kMate's per-block (`--unit ld`) estimator sprays AF in low-diversity regions,
why it is **not a bug**, and the resulting software guidance. Investigation on the
p231 benchmark, n231_g0 (equimolar) and n50_g0 (sparse), production in-house panel.

## TL;DR (software guidance)

The honest, **prior-free** hardening of kMate:

- **(a) `--unit chrom` is the robust default** — proven best on BOTH uniform
  (n231_g0: AF-MAE 0.0033) and sparse (n50_g0: 0.0028, eff_n≈true, ~7% mass on
  absent founders) pools. Pool all of a chromosome's k-mers into one `h`.
- **(b) Never fit low-diversity blocks in isolation.** Per-block (`--unit ld` at
  r²=0.1) isolates the Chr1 centromere (blk7, 14.3–17.3 Mb) into its own block,
  where founders are near-collinear and the fit collapses (AF-MAE 0.021 vs chrom
  0.003). For recombinant pools use **coarse / arm-scale blocks** so no block is
  centromere-only (this is what hapFIRE does).
- **(c) The centromere regions CAN be fixed with `anchor→global` (`prior_weight≈0.3`
  toward the chromosome-wide `h`) — and it's a flag, not a default.** On the n50
  centromere block it recovered AF 0.0079 (vs 0.0167) with only 0.14 mass on absent
  founders (vs uniform-anchor's 0.74). It is safe because it anchors to the pool's
  OWN aggregate composition, not to uniform. **But for a selfing species the
  ancestry is constant genome-wide, so `--unit chrom` already does this — no anchor
  needed.** Reserve `anchor→global` for genuinely recombinant pools (opt-in).
- **(d) NO uniform prior — ever.** On a sparse pool an anchor→uniform dumped **73%**
  of the mass onto founders NOT in the pool and made AF 7× worse. A prior cannot
  distinguish a truly-absent founder from a collapsed one.
- **(e) Surface convergence + eff_n per unit** — never silently project a
  non-converged fit. (Note `conv=False` alone is uninformative — see below.)

## The finding

For a selfing / equimolar pool the true `h` is genome-wide constant. `--unit chrom`
recovers it (eff_n ~205, AF-MAE 0.003). `--unit ld` (r²=0.1 blocks) fits an
independent `h` per block and **collapses in low-diversity blocks** (centromere
eff_n 33.7, pericentromere ~120), spraying AF for records in those blocks
(overall MAE 0.008, 0.6% outliers >0.10 vs chrom's 0.001%).

## Why it is NOT a bug (evidence)

1. **`solve_em_per_block` == standalone `solve_em`** on the centromere block:
   max|Δh| = 7.5e-8. No wiring bug in the block path.
2. **Block-scope kfw is correct.** Using the full-chromosome kfw instead of the
   block's kfw makes it far WORSE (blk0 eff_n 189→37). The per_founder normalizer's
   unit-scope (over the block) is right.
3. **The OLD multinomial normalization collapses HARDER** (centromere eff_n 21.6 vs
   the new per_founder 33.6). The new code did not introduce it.
4. **Every block has K_b = 231** distinct k-mer haplotype patterns (incl. the
   centromere, with 369k informative k-mers) — so it is NOT genuine
   non-identifiability by exact identity; it is **near-collinearity**: the
   centromere's k-mers are redundant (repetitive region), carrying ~30 *independent*
   constraints, so the EM cannot hold 231 apart.
5. **It is the flat-saddle EM behaviour.** Uniform `h` is a fixed point of the
   per_founder EM but an UNSTABLE saddle; the multiplicative EM slides off it toward
   a sparse-vertex optimum at a rate ∝ 1/(information). Whole-chromosome barely
   drifts (eff_n 205→202.6 at 3000 it); the centromere reaches its collapsed
   optimum (conv=True at 2838 it, eff_n 33). So kMate currently relies on
   **early-stopping (`max_iter=200`) as implicit regularization**, which is fragile
   (max_iter-dependent) and fails for low-information units.

## `conv=False` is universal, not a usable guard

At production settings (`tol=1e-7`, `max_iter=200`) **chrom AND all 18 blocks are
`conv=False`** — the multiplicative EM crawls the flat ridge and essentially never
hits `‖Δh‖<1e-7`. So `conv=False` fires everywhere (good blocks included) and cannot
by itself flag bad fits. Report it + `eff_n`, but don't gate on `conv` alone.

## Regularization benchmark (why "no prior" + HARP comparison)

n231_g0 (uniform) and n50_g0 (sparse), CHROM + centromere blk7. `mass_on_ABSENT` =
fraction of `h` on the 181 founders NOT in the n50 pool (flattening tripwire):

| method | n50 CHROM AF-MAE | n50 CHROM mass_absent | n50 blk7 AF-MAE |
|---|---|---|---|
| EM current            | 0.0028 | 0.075 | 0.0167 |
| anchor→uniform w0.3   | 0.0206 | **0.733** | 0.0218 |
| anchor→global w0.3    | 0.0028 | 0.075 | **0.0079** |
| Dirichlet α=1.5       | 0.0027 | 0.062 | 0.0166 |
| HARP-style NNLS       | 0.0065 | 0.059 | 0.0312 |

- **anchor→uniform is catastrophic on sparse pools** (73% mass on absent founders).
- **anchor→global is safe and fixes the centromere** (opt-in flag; see (c)).
- **HARP-style direct NNLS on k-mers is WORSE than the EM** — the direct-solve idea
  does not transfer; kMate's binary-k-mer evidence is the limiter, not the solver.

## How HARP / hapFIRE handle it (verified from source + binary)

- **`harp freq` is a plain multinomial ML EM — NO prior/penalty.** Its binary
  options: `em_iter`, `em_converge` (ML EM), plus prior-free robustness knobs
  `em_random_start_count` / `em_random_start_alpha` (symmetric Dirichlet used only to
  SAMPLE random INITIAL points — not an objective prior) / `em_random_start_seed`,
  keeping `em_run_best`; and `em_min_freq_cutoff` (a small-frequency floor).
- **hapFIRE's ecotype (founder) step** (`haplotype_generation.py`) is a prior-free
  constrained least-squares (`min‖Dᵀh−y‖ s.t. h≥0, Σh=1`, cvxpy/SCS) + length-weighted
  **block-averaging**; a `_lasso` (L1) variant exists but is not the default.
- So HARP's robustness is **richer evidence (per-base read likelihoods) + coarse
  blocks + averaging + multi-random-start**, NOT regularization. kMate should stay
  prior-free and lean on structural conditioning (chrom / coarse blocks).

## HARP's two prior-free mechanisms, tested on kMate — both fail (2026-07-08)

HARP's `harp freq` has two prior-free robustness knobs kMate lacks: multi-random-start
(best-of-likelihood) and `em_min_freq_cutoff`. Tested on the centromere block:

| n231_g0 centromere | eff_n | AF-MAE |
|---|---|---|
| baseline (uniform init) | 33.6 | 0.0197 |
| multi-start best-of-8 | 33.4 | 0.0197 |
| cutoff 1e-3 | 32.2 | 0.0203 |
| cutoff 5e-3 | 19.6 | 0.0293 |

- **Multi-start does nothing**: all 8 Dirichlet-random inits collapse to eff_n 33–34
  with ~identical likelihood. The surface is NOT multimodal — the sparse collapse is
  the genuine global ML optimum, reachable from any init.
- **`em_min_freq_cutoff` only hurts**: zeroing small founders removes more of the
  already-collapsed mass (AF-MAE degrades monotonically).

**Conclusion: the centromere collapse is the maximum-likelihood solution given kMate's
weak/redundant binary-k-mer evidence there. No prior-free mechanism (multi-start,
direct NNLS, min-freq cutoff) fixes it — it is the evidence, not the algorithm.** Only
a prior (anchor→global) recovers it, and that must stay an opt-in flag (never uniform).

## Provenance

Scripts/logs under `benchmarks/p231/scripts/` and `benchmarks/p231/logs/`:
`kb_per_block.py`, `centromere_bughunt.py`, `drift_battery.py`, `conv_map.py`,
`reg_benchmark.py`. Production GrENE-Net already ran `--unit chrom`, so this is
kMate-as-software hardening, not a production re-run.
