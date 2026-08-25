# Notebooks — index

Every rendered notebook in the tree lives here, centralized 2026-08-24. Before
that, results 1 and 3 were in this directory while result 2's sat three levels
down in `phase1_replication/multiaxis/`, which made "where is the notebook for
X" unanswerable.

Notebooks are **outputs**. Each is emitted by a `_build_*_nb.py` builder that
lives in the section that owns the analysis — so the code stays with its
section, and the rendered result is browsable here. To regenerate one, run its
builder (listed below), not the notebook.

Figures written by these notebooks go to `<section>/results/*/plots/`.

---

## Result 1 — negative selection on SVs
Section: `r1_sv_negative_selection/` · detail: that folder's `README.md`

| notebook | what it shows |
|---|---|
| **`temporal_s_consolidated.ipynb`** | **the headline** — per-variant selection coefficient `s`, SV vs p0-matched SNP baseline, plus the climate gradient. Consolidated 2026-07-15 from seven separate notebooks |
| `temporal_s_nofilter.ipynb` | no-filter robustness twin of the above (matched smoothing bandwidth). Deliberate twin, not a stale copy |
| `sv_parallelism_climate.ipynb` | per-site SV−matched-SNP parallelism excess vs climate (bio1 +0.44, bio18 −0.68) — corroborates the climate-slope β arm via an independent statistic. Reduced 2026-08-25 from the former `parallelism_picmin.ipynb`; its bulk-parallelism and PicMin sections are archived |
| `sfs_shift_by_site.ipynb` | folded-SFS shift per site: mean shift null, but **SV extinction rate elevated** — robust across every aggregation, and not a variance artifact |

## Result 2 — GEA hits visible only in non-SNP data
Section: `r2_gea_nonsnp/` (production pipeline in `phase1_replication/`)

| notebook | what it shows |
|---|---|
> **2026-08-25 — the live track is raw per-record LFMM.** WZA block aggregation and the
> quasi-binomial / kendall models were retired to `../archive_gea/` (see its README for
> why: WZA is ~170× anti-conservative here, and honest recalibration collapses ~7,200
> block-hits to 6). 21 notebooks moved with them and are **no longer in this directory**.
> Everything listed below is live and reads the raw per-record tables.

| notebook | what it shows |
|---|---|
| `raw_manhattan_snp_vs_{sv,nonsnp}_lfmm_tile.ipynb` (+ `_site_tile`) | **the headline figures** — per-variant mirror Manhattans on the gap-free tiling partition, 22 climate axes (bio1-19 + pc1-3), each with a QQ panel reporting λ |
| `gif_manhattan_snp_vs_{sv,nonsnp}_lfmm_tile.ipynb` | the same, GIF-corrected (λ divided out where λ>1) |
| `newpeak_dotgrid_lfmm_{sv,nonsnp,smallindel}.ipynb` | new-peak dot-grids — significant blocks per class across climate axes (raw, block-lead) |
| `nonsnp_peak_dotgrid_lfmm_tile.ipynb` | non-SNP peak dot-grid + the raw per-record Bonferroni/FDR count table, tiling partition |
| `snp_vs_nonsnp_peaks_viz.ipynb` | which peaks the non-SNP scan finds that the SNP scan misses (GWAS-side) |
| `persite_new_peaks_{nonsnp,sv}.ipynb` | the same question per garden, on the production tiling blocks |
| `persite_new_peaks.ipynb` | ⚠ **superseded** — the pre-tiling version; its partition silently dropped ~28% of variants to inter-block gaps. Kept only as the predecessor |
| `block_gene_significance_overlap.ipynb`, `nonsnp_block_characterization.ipynb` | block→gene attribution and characterization of the non-SNP blocks |
| `lfmm_k_calibration.ipynb` | the CAM5 replication (`cam5_replication/`) |
| `10_lfmm_k_selection.ipynb` | why K=16 latent factors — the standing K decision, and still load-bearing for the live LFMM |

## Result 3 — per-site GWAS on the ecotype-selection trait
Section: `r3_persite_gwas/`

| notebook | what it shows |
|---|---|
| `class_gwas_persite.ipynb` | the 31 individual per-garden LOCO-EMMAX scans, SNP vs non-SNP block overlap |
| `class_gwas_multitrait.ipynb` | Bolormaa multi-trait meta across sites (JOINT / GLOBAL / CLIMATE) |
| `sv_indel_tagging_masked.ipynb` | SV/indel→SNP tagging r² (masked, MAC-floored) — the "is the non-SNP layer redundant" control |
| `nonsnp_only_genes.ipynb` | annotated genes under blocks the non-SNP scan flags and the SNP scan misses |

## Method sections

**`wza/`** — ⚠ **RETIRED 2026-08-25, moved to `../archive_gea/`.** The whole section
(`wza_script.py` + `investigation/`) and its notebooks — `wza_investigation.ipynb`,
`wza_manhattan_cap_vs_nocap.ipynb`, `wza_sd_fix_test.ipynb`, `wza_sd_fit_audit.ipynb`,
`cap_poly_decision.ipynb`, `fit_inspection_{kendall,lfmm,binomial}.ipynb`,
`pc1_manhattan.ipynb` — now live under `archive_gea/`. `manhattan_pc1.ipynb` stayed
(it reads no WZA output). Also retired there: the kendall / quasi-binomial model
notebooks (`kendall_fix_compare`, `binomial_fix_compare`, `manhattan_*_wza`,
`manhattan_{3models_deg7cap2000,clq90_deg2}`, `nonsnp_bonf_overlap`,
`snp_vs_nonsnp_new_peaks`, `sv_polarity_enrichment_clq90`).

**`blocks/`** — what is a test unit?
`block_coherence_clqcut.ipynb`, `blocks_units_decision.ipynb`,
`block_breakage_window_vs_global.ipynb` (retired window-mode subsystem, kept as
record), `haploblock_r20{10,20}_eps0_validation.ipynb`.

**`qc/`** — QC of this analysis.
`qc_coverage_audit.ipynb`, `panel_stats_arch3.ipynb`,
`panel_overlap_grenenet.ipynb`, `founder_9977_kmer_space.ipynb`,
`seedmix_kmate_vs_hapfire.ipynb`,
`kmate_founder_fix_results.ipynb`, `kmate_kfw_fullpanel_results.ipynb`.

**`extras/`** — `sv_selection_haplotype_audit.ipynb` (haplotype-unit SV
enrichment; the apparent signal dies under a genome-rotation null).

---

## Notes

- **A notebook with no obvious builder is usually not an orphan.** Several
  builders construct the filename dynamically (`f"newpeak_dotgrid_lfmm_{cls}.ipynb"`),
  so a plain grep for the filename finds nothing. Check the section's
  `_build_*_nb.py` before concluding a notebook is unreproducible.
- `seedmix_kmate_vs_hapfire.ipynb` and `block_breakage_window_vs_global.ipynb`
  genuinely have no live builder — the first has none in the repo, the second's
  is inside `archive/window_hapfreq_retired/`.
- Numbered notebooks (`06_`, `10_`) are survivors of the original gen-3 pilot
  series; `01`–`05` and `07`–`09` were retired with that generation on
  2026-08-24 (`archive/retired_2026-08-24/gen3_pilot/`).
