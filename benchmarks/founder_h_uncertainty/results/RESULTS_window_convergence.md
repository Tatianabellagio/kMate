# Window-mode ĥ: identifiability is far worse locally, and coverage now matters

**Job 35277452** · `diag_window.py` · p80, pool `cov10_n50_g0_s42_hotspots_p80_chr1`
(g0 ⇒ local ancestry = global, so per-window truth = global h_true) · 10 kb windows,
3043 windows (2886 with k-mers), per-window EM (`ω=1/m_b`), 3 Poisson reps/coverage.

## Result — per-window ‖ĥ_w − h_true‖₂

| condition | median | mean | p90 | frac > 0.1 |
|---|---|---|---|---|
| noiseless (∞) | **0.167** | 0.182 | 0.228 | 0.99 |
| 3× | 0.174 | 0.191 | 0.235 | 0.99 |
| 10× | 0.169 | 0.186 | 0.230 | 0.99 |
| 30× | 0.167 | 0.184 | 0.231 | 0.99 |

**Global-mode reference (whole chromosome): 0.0485, flat in coverage.**

Coverage effect on windows with ≥200 k-mers (identifiability subtracted):

| coverage | median(err_cov − err_noiseless) |
|---|---|
| 3× | **+0.007** |
| 10× | +0.002 |
| 30× | +0.001 |

`corr(log10 n_kmers, err_cov3) = −0.686` · 5% of windows below the 200-k-mer floor.

## Conclusions

1. **Per-window h is much less identifiable than global.** Even noiseless / infinite
   coverage, per-window error is **0.167** — ~3.5× the global 0.0485, with 99% of
   windows > 0.1. A 10 kb window's few thousand k-mers cannot pin 80 founders. So
   the window estimator trades global's identifiability for local resolution, and
   the per-window h is intrinsically noisy/biased.
2. **Coverage now matters — unlike global.** Lower depth monotonically inflates error
   (+0.007 at 3×, +0.002 at 10×, +0.001 at 30× over the noiseless floor). The
   mechanism is **support collapse**: at low λ a window's few k-mers get zero counts,
   shrinking effective support. This is exactly why low-coverage windows lack the
   evidence to fit local h (and is the mechanism behind NaN'ing low-info windows).
3. **Local k-mer support is the dominant lever.** `corr(log10 n_kmers, err) = −0.69`:
   windows with more local k-mers are far better resolved. Coverage acts largely
   *through* support (fewer nonzero k-mers at low depth), so a support floor
   (`--min-kmers-per-block`) is well-motivated, but note the floor controls
   identifiability, while coverage modulates how much of that support survives.

## Implication for uncertainty quantification

- A per-window h SE that captures only coverage/sampling variance would miss the
  dominant ~0.167 **identifiability** term (same blind spot as global, worse in
  magnitude). An honest per-window confidence must reflect **local identifiability**
  (k-mer support, local founder resolvability / Fisher condition number), not just
  depth. The coverage term (+0.001…+0.007) is real but secondary except where it
  drives support toward the floor.
- Whether this large per-window h error damages the **AF functional**
  (`AF_r = ĥ_wᵀv_r/ĥ_wᵀu_r`, which sums over carriers and may be AF-neutral to the
  non-identifiable swaps) is the next thing to measure — it determines whether the
  reportable uncertainty is an h-SE problem or an AF-SE problem.

*Data: `window_diag.tsv` (per-window n_kmers, local-true-founders, errors at each
coverage).*
