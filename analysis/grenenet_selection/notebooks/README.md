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
| `parallelism_picmin.ipynb` | replicate-plot parallelism + PicMin + per-site parallelism vs climate |
| `sfs_shift_by_site.ipynb` | folded-SFS shift per site: mean shift null, but SV extinction rate elevated |
| `sfs_time_site4.ipynb` | AF spectrum across generations 0→3 at site 4, by variant class |

## Result 2 — GEA hits visible only in non-SNP data
Section: `r2_gea_nonsnp/` (production pipeline in `phase1_replication/`)

| notebook | what it shows |
|---|---|
| `newpeak_dotgrid_lfmm_{sv,nonsnp,smallindel}.ipynb` | the new-peak dot-grids — significant blocks per class across climate axes |
| `wza_manhattan_snp_vs_{sv,nonsnp}_lfmm_final.ipynb` | block-level WZA Manhattans on the settled regime (clq0.9 tiling · isotonic SD · deg-5 clamped mean · no cap) |
| `raw_manhattan_snp_vs_{sv,nonsnp}_lfmm_tile.ipynb` (+ `_site_tile`) | per-variant Manhattans on the gap-free tiling partition |
| `gif_manhattan_snp_vs_{sv,nonsnp}_lfmm_tile.ipynb` | the same, GIF-corrected |
| `nonsnp_peak_dotgrid_lfmm_tile.ipynb` | non-SNP peak dot-grid, tiling partition |
| `snp_vs_nonsnp_new_peaks.ipynb`, `snp_vs_nonsnp_peaks_viz.ipynb` | which peaks the non-SNP scan finds that the SNP scan misses |
| `persite_new_peaks_{nonsnp,sv}.ipynb` | the same question per garden, on the production tiling blocks |
| `persite_new_peaks.ipynb` | ⚠ **superseded** — the pre-tiling version; its partition silently dropped ~28% of variants to inter-block gaps. Kept only as the predecessor |
| `block_gene_significance_overlap.ipynb`, `nonsnp_block_characterization.ipynb`, `nonsnp_bonf_overlap.ipynb` | block→gene attribution and characterization of the non-SNP blocks |
| `cap_poly_decision.ipynb`, `wza_sd_fit_audit.ipynb` | how the WZA correction regime was chosen |
| `manhattan_3models_deg7cap2000.ipynb`, `manhattan_clq90_deg2.ipynb`, `manhattan_lastgen_wza.ipynb` | earlier-regime Manhattans (deg-7/deg-2); superseded by the `_final` set above for citation |
| `kendall_fix_compare.ipynb`, `binomial_fix_compare.ipynb` | before/after checks on the kendall and binomial fixes |
| `sv_polarity_enrichment_clq90.ipynb` | SV enrichment / insertion-deletion polarity across the three models |
| `lfmm_k_calibration.ipynb`, `manhattan_{byclass,combined}_wza.ipynb` | the CAM5 replication (`cam5_replication/`) |
| `10_lfmm_k_selection.ipynb` | why K=16 latent factors — the standing K decision |

## Result 3 — per-site GWAS on the ecotype-selection trait
Section: `r3_persite_gwas/`

| notebook | what it shows |
|---|---|
| `class_gwas_persite.ipynb` | the 31 individual per-garden LOCO-EMMAX scans, SNP vs non-SNP block overlap |
| `class_gwas_multitrait.ipynb` | Bolormaa multi-trait meta across sites (JOINT / GLOBAL / CLIMATE) |
| `sv_indel_tagging_masked.ipynb` | SV/indel→SNP tagging r² (masked, MAC-floored) — the "is the non-SNP layer redundant" control |
| `nonsnp_only_genes.ipynb` | annotated genes under blocks the non-SNP scan flags and the SNP scan misses |

## Method sections

**`wza/`** — was WZA the right block-aggregation choice?
`wza_investigation.ipynb` (the main writeup), `wza_manhattan_cap_vs_nocap.ipynb`,
`wza_sd_fix_test.ipynb`, `fit_inspection_{kendall,lfmm,binomial}.ipynb`,
`pc1_manhattan.ipynb`, `manhattan_pc1.ipynb`.

**`blocks/`** — what is a test unit?
`block_coherence_clqcut.ipynb`, `blocks_units_decision.ipynb`,
`block_breakage_window_vs_global.ipynb` (retired window-mode subsystem, kept as
record), `haploblock_r20{10,20}_eps0_validation.ipynb`.

**`qc/`** — QC of this analysis.
`qc_coverage_audit.ipynb`, `panel_stats_arch3.ipynb`,
`panel_overlap_grenenet.ipynb`, `founder_9977_kmer_space.ipynb`,
`06_sv_support_filter.ipynb`, `seedmix_kmate_vs_hapfire.ipynb`,
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
