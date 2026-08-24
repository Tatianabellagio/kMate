# Retired GWAS notebooks (2026-07-08)

Retired during the post-Kf_w notebook cleanup. All are superseded — the live
SNP-vs-non-SNP-vs-SV story now lives in three notebooks under `notebooks/`:
`class_gwas_multitrait.ipynb`, `class_gwas_persite.ipynb`, and the all-marker
`multisite_founder_gwas.ipynb`.

| File | Why retired | How to rebuild if needed |
|---|---|---|
| `class_split_gwas.ipynb` | Overview/peak-overlay of SNP vs non-SNP (2-way only). Its peak-overlap question was folded into the *end* of `class_gwas_multitrait.ipynb` and `class_gwas_persite.ipynb` as **two significance-threshold overlap tables** (Bonferroni 0.05 and BH-FDR q<0.05), 3-way (snp/nonsnp/sv), at the clq0.9-block level — peak-specific, unlike the original whole-genome 20 kb window Spearman (which was dominated by ~5.5k mostly-null windows). | its builder is here: `_build_class_split_gwas_nb.py` |
| `multisite_founder_gwas_clq50.ipynb` | Block-clustering (r²≥0.5) variant of the all-marker multisite GWAS; stale (pre-Kf_w) and redundant with the marker-level `multisite_founder_gwas.ipynb`. | `run_multisite_downstream.sh clq50 _clq50` |
| `multisite_founder_gwas_clq90.ipynb` | Same, r²≥0.9 blocks. | `run_multisite_downstream.sh clq90 _clq90` |
| `multisite_founder_gwas_clq90_pc1.ipynb` | Same, r²≥0.9 blocks + bioclim-PC1 climate axis. PC1 is already a column in `class_gwas_multitrait.ipynb`. | `run_multisite_downstream.sh clq90 _clq90` with the PC1 axis env |
| `founder_gwas231.ipynb` | Single-site (site 4 only) all-231-founder scan; legacy, no class split. | `_build_foundergwas_nb.py` (kept in place under `analysis/grenenet_gea/`) |

Rebuild tooling kept live: the parameterized `_build_multisite_gwas_nb.py`
(builds the main notebook by default) and `_build_foundergwas_nb.py` remain in
`analysis/grenenet_gea/`. The compute script `class_split_gwas.py` (produces the
`class_gwas_{snp,nonsnp,sv}.npz` that the live notebooks read) is unchanged and
also stays live.
