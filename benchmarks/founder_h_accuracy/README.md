# benchmarks/founder_h_accuracy/

**Question:** how accurately does kMate recover the founder mixture `h` itself
(as opposed to the per-record AF it projects from `h`)?

This is the benchmark most sensitive to the estimator internals — the per-founder
M-step normalization and haploblock collapse change the founder *decomposition* more
than they change AF, which is why ROADMAP_GLOBAL_REFRESH Phase 3 designates it the
highest-signal refresh target.

kMate-only: no competitor arm. For the kMate-vs-hapFIRE founder comparison see
`../vs_hapfire/` (h R²/RMSE, native panels) and `../vs_hapfire/ecotype_resolution/`
(set-resolution metrics).

## Scoring

`scripts/score_h_vs_truth.py --est <h_per_chrom.npz> --truth <sim>/pool_weights.tsv --label <tag>`
(run in the `basic` env). Metrics: per-founder h RMSE/MAE, **absorbed count**
(founders driven to ~0), cactus-vs-PanGenie mass balance, est-vs-true slope.

The absorbed count and the cactus/PG mass ratio are the diagnostics that caught the
pre-2026-07-06 founder-`h` collapse (19–49 of 231 founders driven to ~0). Under the
corrected estimator the expectation is **0 absorbed, cactus mass ratio ≈ 1.0,
slope ≈ 1.0** — if a run misses that, stop and investigate before trusting anything
downstream of it.

## Contents

    scripts/score_h_vs_truth.py       the scorer
    scripts/run_h_unit_chrom.sbatch   the --unit chrom runner
    results/BASELINE_n231_g0.csv      retained reference point for the old-vs-new delta
    results/ldr01_<regime>_*.{npz,json,csv}   --unit ld r²=0.1 arm (per-block / SNP-parity control)
    results/<regime>_<arm>_*.{json,csv}       the filt2mb / filt2u / filt2invu / raw arms

Arms scored here also include the earlier `filt2mb` / `filt2u` / `filt2invu`
k-mer-filter × weighting combinations, kept for the same delta purpose.
