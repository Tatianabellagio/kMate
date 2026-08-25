# benchmarks/ldblock_window_test/

**Question:** at what block granularity can window-mode kMate actually get a *local*
fit, and what does that cost in AF accuracy?

The block-size / k-mer-floor sweep behind the `dynld_K500` unit map. Its motivating
result is the **ghost-block problem**: on the old median-7-variant CLQ0.9 LD blocks
only ~42% of blocks got a local fit (≥50 observed k-mers) — ~30% coverage-limited,
~28% zero-panel-k-mer. p80 ≈ p231, so the cause is **block thinness, not panel
imbalance**. Growing each block along LD until it holds ≥500 panel k-mers
(`dynld_K500`) lifts the local-fit rate to ~80.5%.

## ⚠️ Status

**Window-mode work, and window mode is out of scope for the 2026-07-07 estimator
refresh** (ROADMAP_GLOBAL_REFRESH). Everything here predates the corrected estimator
(per-founder normalization + haploblock collapse), so its numbers are **not
comparable to refreshed `--unit chrom` results**. See the marked table in
`../accuracy_vs_competitors/BENCHMARK_DESIGN.md` §5, where both the "global" and
"dynld-window" columns are flagged pre-refresh for this reason.

Do not resurrect the "window beats global" conclusion from these numbers without
re-running both arms.

## Scoring

`score_ldblk.py` — two modes:

    --validate <out.tsv> <truth.gz> <hblocks.npz>   one run: % windows locally fit
                                                    (status==0) + MAE/R²/r on SNPs
                                                    and all records
    --sweep                                         score all sweep_s*_*.tsv,
                                                    aggregate per k-mer floor

The **% locally fit** number is the point of this benchmark; AF accuracy alone hides
whether a block was fit locally or fell back.

## Contents

    run_coarse_sweep.sbatch      the sweep driver
    score_ldblk.py               scorer (both modes)
    coarse_sweep_scores.csv      aggregated sweep result
    *_mkb50.{tsv,npz,log}        per-run window outputs at the mkb50 floor
    blockpanel_dynld_outcross_g3.png   the dynld-vs-global outcross figure

The `dynld_K500` map itself lives at
`analysis/grenenet_selection/blocks/results/blocks_mcf90/chr{N}_units_dynld_K500.tsv`.
