# Retired 2026-08-25 — the in-house per-site GWAS presentation layer

Retired when the per-garden GWAS moved onto GEMMA and the cross-site meta was dropped
(scope decision: per-garden only, the all-sites question belongs to the GEA track).

| retired | why | replaced by |
|---|---|---|
| `notebooks/class_gwas_multitrait.ipynb` | **obsolete** — it is entirely the Bolormaa JOINT/GLOBAL/CLIMATE meta, which is no longer computed | nothing; the question moved to `r2_gea_nonsnp/` |
| `notebooks/class_gwas_persite.ipynb` | superseded — built on the in-house EMMAX/P3D `class_gwas_*.npz` | `notebooks/persite_gwas.ipynb` |
| `code/_build_class_gwas_multitrait_nb.py` | builder for the obsolete meta notebook | — |
| `code/_build_class_gwas_persite_nb.py` | builder for the superseded per-site notebook | `r3_persite_gwas/_build_persite_gwas_gemma_nb.py` |
| `code/plot_class_gwas_pngs.py` | rendered 319 PNGs, mostly multitrait contrasts that no longer exist | the notebook builder renders grid figures to `results/gemma_gwas/plots/` |

The 319 rendered PNGs moved to `results/varexp/gwas_plots_retired_2026-08-25/` (results are
gitignored, so they are on disk only).

Per LAYOUT.md rule 8, paths inside these files are NOT rewritten — they were written against the
layout and the estimator of their time, and repointing them would falsify the record.
