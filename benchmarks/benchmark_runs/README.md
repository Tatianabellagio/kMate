# benchmarks/benchmark_runs/

Shared **per-pool kMate output store** (~11 GB, gitignored) for the p80/p231
benchmark grid. Not a benchmark itself — it is the estimator-output layer that
several scorers read from.

    tsv/                 per-pool AF tables + h npz (232 files)
                           <pool>_global.tsv            --block-mode global era
                           <pool>_global.h_per_chrom.npz
                           <pool>_block_dynldK500.tsv   window/dynld era
                           <pool>_block_dynldK500.h_blocks_per_chrom.npz
    rows/                per-pool scored metric rows (58)
    rows_competitors/    competitor-side scored rows

## ⚠️ These outputs are PRE-REFRESH — check the timestamp before using any of them

The `*_global.tsv` files were last written **2026-07-07 08:44**. That is *after* the
per-founder normalization landed (`a8ba02d`, 07-06 22:06) but *before* three further
changes the same day:

| commit | time | change |
|---|---|---|
| `9669be7` | 07-07 10:22 | full-panel `Kf_w` (fixes observed-only survivorship bias) |
| `b4d6ce0` | 07-07 16:15 | haploblock collapse |
| `a9bf1c0` | 07-07 18:12 | estimator unified under `--unit` |

So they look current by date but are not. This is the trap that made the 4-tool
competitor tables stale — see ROADMAP_GLOBAL_REFRESH Phase 5.

**Refreshed kMate output lives elsewhere**, per the roadmap's "new output dirs, don't
overwrite" rule:

    ../p80/results/kmate_chrom_p80_{raw,filt2inv}/<regime>/
    ../p231/results/kmate_chrom_*/<regime>/

## Still read by

`../accuracy_vs_competitors/scripts/{plot_block_nofallback_af.py,
run_hapfire_pool.sbatch,stratify_vg_multiplicity.py}`. That live consumption is why
this dir is retained rather than archived — see the "not archived despite looking
stale" list in `../archive/README.md`.
