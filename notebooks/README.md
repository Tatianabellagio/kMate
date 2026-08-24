# notebooks/ — cross-cutting notebooks

Notebooks that don't belong to a single subproject. Anything owned by a subproject lives
with it (`analysis/<topic>/notebooks/`, `benchmarks/<topic>/`, `preprocess_qc/notebooks/`).

## Current contents

| Notebook | What it shows |
|---|---|
| `blocks_snp_sv_density.ipynb` | SNP and SV density across the LD-block partition. Feeds the block/unit definition work in `analysis/grenenet_gea/blocks/`. |

## What used to be here

The 35 panel-development notebooks and their `plots/` subdir moved to
`archive/notebooks_panel_dev_2026-05/` on 2026-08-24, and the 100 root-level `plots/` PNGs
they rendered moved to `archive/plots_panel_dev_2026-05/`.

They were all last touched 2026-05-28 / 2026-06-03 and built on the **pre-arch3**
`cn_full` / `cn_var` / v3qc panels — `archive/cn_full_superseded/` and
`archive/deprecated_v3qc/` already record those panels as superseded. None was referenced
from any doc. See `archive/README.md` for the full rationale; nothing was deleted.

If you're looking for a `CN_FULL_COMPOSITION_*`, `AF_TRUTH_VS_ESTIMATE_v3qc_*`,
`FILT2_RESULT`, `FINAL_RESULTS*`, `RECOMB_SWEEP_RESULTS*`, `V3QC_SANITY_CHECK`, a
carrier/bubble diagnostic, or the h-imbalance set — it's in
`archive/notebooks_panel_dev_2026-05/`.
