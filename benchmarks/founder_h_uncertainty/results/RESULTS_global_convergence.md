# ĥ error in GLOBAL mode is identifiability bias, not coverage

**Job 35275967** · `diag_convergence.py` · p80 panel, pool `cov10_n50_g0_s42_hotspots_p80_chr1`
(50-founder g0 = no-recombination mixture) · production EM (`ω_k = 1/m_b`).

## Setup (what was tested)

Clean **generative** test in **GLOBAL mode** (one EM over the entire Chr1 k-mer set,
8,108,672 k-mers, F=80). Counts are drawn straight from the model the EM assumes —
`c_k ~ Poisson(λ · μ_k(h_true))`, `μ_k = h_trueᵀ K` — so the model is *perfectly
specified* (no reads/mapping/overdispersion). This isolates the EM's intrinsic
recoverability of ĥ from any read-level confounder. The "noiseless" row uses
`c_k = μ_k(h_true)` exactly (coverage = ∞).

## Result — ‖ĥ − h_true‖₂ is flat across coverage AND precision

| coverage | precision | ‖ĥ−h_true‖₂ | mass off true support | n̂ support |
|---|---|---|---|---|
| ∞ (noiseless) | f32 | 0.0485 | 0.0140 | 52 |
| ∞ (noiseless) | f64 | 0.0485 | 0.0134 | 52 |
| 3× | f32 | 0.0485 | 0.0158 | 52 |
| 10× | f32 | 0.0482 | 0.0147 | 52 |
| 10× | f64 | 0.0485 | 0.0133 | 52 |
| 30× | f32 | 0.0488 | 0.0121 | 52 |
| 100× | f32 | 0.0486 | 0.0133 | 52 |

True support = 50 founders.

## Conclusions

1. **Not coverage.** Error is identical from 3× to **infinite** coverage (0.0485
   everywhere). With perfect noise-free data the EM is still 0.0485 off h_true. In
   global mode, coverage is **not** the bottleneck — the whole-chromosome k-mer set
   (8.1M k-mers) over-determines the 80-founder simplex regardless of depth.
2. **Not float32.** f32 ≡ f64 (0.0485). Precision is not the issue.
3. **Not under-iteration.** ‖Δh‖ plateaus ~1e-6 by iter 2000; the value is stable.
   The `converged=False` flag is **cosmetic**: production `tol=1e-7` on an 80-vector
   sits below the reachable float floor, so the loop runs to `max_iter` after ĥ has
   already stopped moving. (Worth relaxing the production tol or switching to a
   relative/￼value-stability criterion to stop wasting iterations.)
4. **It is identifiability / collinearity bias.** ĥ converges to a point on a
   near-flat likelihood ridge ~0.0485 from truth: ~1.3% of mass leaks onto absent
   founders (`n̂ support` 52 vs true 50), and `l2_on_true_support ≈ 0.047 ≈` the
   total error, i.e. most error is **mis-apportionment among collinear *present*
   founders** the k-mers cannot distinguish.

## Implication for uncertainty quantification

Variance-based SEs (observed Fisher, parametric Poisson bootstrap) describe scatter
*around* ĥ and **cannot capture this bias** — which is why the parametric bootstrap's
nominal 95% interval covered h_true only **6%** of the time at 3× (job 35271703).
Per-founder ĥ is **not identifiable** here (matches `BACKGROUND.md` limitation #4:
"per-ecotype recovery is noisy; aggregates accurate"). The identifiable, reportable
quantity is the AF functional `AF_r = ĥᵀv_r / ĥᵀu_r`, which sums over carriers and is
robust to AF-neutral swaps within the non-identifiable subspace.

## Note on global vs window mode

This flatness is **specific to global mode**, which pools k-mer evidence across the
entire chromosome. **Window mode** fits each ~10 kb window from only its local
(thousands of) k-mers, so (a) local identifiability is far worse and (b) coverage
*will* matter — low depth means windows lack the support to pin local h (and is the
mechanism behind low-info windows being NaN'd). Quantified separately in the
window-mode diagnostic.

*Data: `convergence_diag.tsv`, `convergence_traces.npz`. Useful as a paper supplemental
("recoverability of the founder mixture is identifiability- not coverage-limited in
global mode").*
