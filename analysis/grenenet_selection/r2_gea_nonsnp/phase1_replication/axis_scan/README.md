# axis_scan — site-level PC1 / bio-axis LFMM climate scan

Ported 2026-07-21 from the retired `gea_newpanel` snp/nonsnp fork
(`../../archive/gea_newpanel_snp_nonsnp_fork_retired/`). This is the one
genuinely-new capability that fork carried and that the main clq0.9 pipeline
does **not** cover.

## What it is (and how it differs from the main pipeline)

- **Unit:** 31 flower-weighted **site means**, LFMM **K=3**, site-MAF ≥ 0.05 —
  NOT the plot-level gen9 K=16 unit used by `run_lfmm_lastgen.R` / the
  `run_3class_*` drivers and `multiaxis/`.
- **Classes:** SNP / non-SNP / SV (3-panel Manhattans + a QQ grid showing
  per-axis×class raw-vs-GIF calibration).
- **Axes:** the signal-bearing subset — bio5, pc1, bio12, bio13, bio16, bio19.
- **Blocks:** uses `blocks_clq09.assign_clq09_blocks` (the older
  `grenenet_gea/blocks_recompute/*_clq0.9_blocks_clq0.9.tsv` partition), not
  `blocks_mcf90`. Kept as-is so the ported analysis reproduces faithfully; treat
  its block ids as NOT comparable to the main pipeline's mcf90 blocks.

## Files

- `build_site_matrices.py` → `results/lfmm_site/lfmm_{cls}_site_Y.f64`
- `build_power_inputs.py`, `power_experiments.sbatch`, `run_lfmm_both_p.R`,
  `exp{1,2,3}_*.sbatch` → the site LFMM (raw + GIF p) over the axes → `results/power/`
- `report_power.py` → aggregates the `*_bothp*.csv`
- `_build_axis_scan_nb.py [pval_gif|pval_raw]` → builds
  `notebooks/axis_scan_manhattan_qq*.ipynb`, writing figures + `sig_blocks_axis_scan*.csv`
  into `results/axis_scan/`.
- `blocks_clq09.py` → the older CLQ-0.9 block assigner (external BLOCK_DIR path).

Input data trees (`results/{lfmm_site,env_site,env_plot,power}/`) were moved here
from the fork and are gitignored.
