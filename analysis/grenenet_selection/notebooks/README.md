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
| **`sv_founder_confound.ipynb`** | **is the climate gradient real?** — measured as **kb of sequence absent from Col-0** (median 121 kb). Cold-origin ecotypes carry more inserted sequence (ρ=−0.542) and those ecotypes decline as gardens warm. Four confounds tested: panel half (ruled out for kb, **real for counts**), provenance (real but plateaus), Col-0 reference bias (ruled out), panel construction (ruled out — graph representation, assembly size, and **sequencing platform** all fail to explain it; holds at full strength in the 151 founders with **no assembly**). Deletions null throughout; maps content in geography + origin-climate space |

| **`sv_insertion_content.ipynb`** | **what is in the insertions?** — all 172,220 SV insertions (no MAC floor, 388.6 Mb). Insertions are **5.2× depleted in CDS** (log2 obs/exp −2.37), a clean purifying gradient through UTR → intron → intergenic. 89.7% lift back into a carrier's assembly; cargo read from the precomputed Helixer / Liftoff / TRASH tracks. Organellar contamination screen clean (0.14%). Carrier-climate split is **frequency-matched** — TE overlap is consistently higher for cold-origin insertions (5/5 strata), CDS share is not (flips 3/5). §4 aligns the inserted sequence back to TAIR10 (dc-megablast): **52.9% TE-derived** across 314 families (ATREP/Helitron-led), 28.3% duplicated Col-0 sequence, and `novel` is an alignment-power floor that collapses 64.1% (<100 bp) → **1.0%** (>10 kb) — read cargo **within a size bin**. The cold-origin excess is a **TE excess** (42.3% vs 32.7% at n_alt 8-20, cold > warm 5/5 strata) but **broad-based** — only 16/33 families cold-skewed. §5: one outlier, **ATCOPIA78 (ONSEN)** at cold:warm 2.3×, the only Bonferroni survivor and the youngest cargo present (99.4% identity to TAIR10 vs 94.6%, median 4,961 bp). The Helixer TE-ORF caveat was tested and **runs the other way** — TE-derived has the *lowest* de novo gene rate (22.4% vs 28.0% overall) |

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
| `class_signal_gain.ipynb` | significant blocks SNP / non-SNP / SV and their overlap, GEA **and** per-garden GWAS on one block rule (re-assigned unmerged tiling), split shared / SNP near-miss / SNP-blind, per axis, per garden, GEA vs GWAS. **§9 = the raw-vs-GIF arithmetic**: union over 22 axes, GIF (χ²/λ, only where λ>1) leaves **0** Bonferroni blocks for snp/nonsnp/smallindel and **4** for SV (BH: 9 SV, 3 nonsnp, from sv·bio6, sv·pc1, nonsnp·bio6 only), against raw Bonferroni 1,430 snp / 1,312 nonsnp / 273 sv — while raw BH goes the other way and calls 66–69% of blocks. Builder `genes/convergence/_build_class_signal_nb.py`, numbers `genes/convergence/class_signal_counts.py` + the GEA production table `raw_block_significance_lfmm{,_union}.csv` |
| **`flc_case.ipynb`** | **a worked case: is our FLC hit real, or made by the panel?** Written for someone outside the project. **§1** the locus in both scans, p-values AND the mean Δp behind them: ours is Bonferroni-significant on bio6 (best −log10 p 9.45 vs a 7.15 bar; 6 non-SNP + 10 SNP records), the phase-1 scan — **same model** (`lfmm_ridge` K=16), same 355/352 last-gen pools, same climate — reaches 5.01 against a 7.32 bar. Mean Δp per record agrees only at r = 0.63 over the 225 shared SNPs. **§2** rules out ascertainment, frequency (p0 MAF 0.22–0.45), imputation (76% directly called in raw 1001G), genotype error (98.6% vs the assemblies) and calibration (their λ 1.99 > our 1.81). **§3** the two pipelines' **founder mixtures agree — median r = 0.991** over 352 pools, so the fault is in the *projection*, not the mixture. **§4 the decisive test**: kMate's OWN `h` with the 153 missing genotypes filled from the phase-1 panel → agreement with hapFIRE **r 0.68 → 0.997**, and the garden-level bio6 gradient falls from **−0.673 (p 3e−5) to −0.547 (p 1.5e−3)**, onto phase-1's −0.538. **The 78-founder restriction is what makes it genome-wide significant** — all three estimates still show a real gradient at p≈1e−3, so the effect is plausibly real, the *significance* is not. **§5** genome-wide, the 78-only stratum runs ~2× the non-SNP and up to 12× the SV hit rate of well-called records within the same MAF band (λ 2.1–2.3 vs 1.7–1.8); ~**12% of non-SNP and 22% of SV hits** are excess → a **call-rate floor of 0.50** (0.90 would cost 56%/82% of tests for no measurable gain). Builder `genes/convergence/_build_flc_case_nb.py`, data `genes/convergence/flc_case_data.py` + `r2_gea_nonsnp/phase1_replication/multiaxis/callrate_calibration.py` |
| `lfmm_k_calibration.ipynb` | the CAM5 replication (`cam5_replication/`) |
| `10_lfmm_k_selection.ipynb` | why K=16 latent factors — the standing K decision, and still load-bearing for the live LFMM |

## Result 3 — per-site GWAS on the ecotype-selection trait
Section: `r3_persite_gwas/`

| notebook | what it shows |
|---|---|
| **`persite_gwas.ipynb`** | **the results view** — 30 per-garden GEMMA scans: λ + hit counts, QQ grids, the MAC-stratified QQ, Manhattan grids, and the hit list with MAC attached. Builder: `r3_persite_gwas/_build_persite_gwas_gemma_nb.py` |
| ~~`class_gwas_persite.ipynb`~~ | ⚠ **retired 2026-08-25** → `archive/retired_2026-08-25/notebooks/`. Built on the in-house EMMAX/P3D arm; superseded by `persite_gwas.ipynb` |
| ~~`class_gwas_multitrait.ipynb`~~ | ⚠ **retired 2026-08-25** → same place. It was entirely the Bolormaa JOINT/GLOBAL/CLIMATE meta, which is no longer computed (per-garden-only scope; the all-sites question is the GEA track's) |
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

| notebook | what it shows |
|---|---|
| **`pangenie_accuracy.ipynb`** | **how good is the PanGenie side of the 231-founder panel** — four figures. **Coverage**: of the 8.49 M records the 78 assembly founders reveal, PanGenie calls **90.3%** in >=1 of the 153 short-read founders (mean call rate 80.6%), falling with size — **90.7% -> 70.1% callable** and **81.8% -> 52.8% call rate** from SNPs to >5 kb SVs. Genome-wide that concentrates in 131 regions, 3.11 Mb (**2.60% of the genome**, 445 genes), **80% in the chromosome arms** not centromeres, TE-dense; FLC is one. **Accuracy** (leave-one-out, 73 founders — 5772/9947 excluded as dropped duplicate assemblies): read **nRD**, not GC (GC counts the ref/ref majority and sits at 99.5%+ regardless). The harness scores RAW output; production applies het->missing + GQ>=20, which **thirds** the error — median nRD **1.7% SNP / 1.9% small indel / 0.94% small SV / 0.41% medium SV / 0.15% large SV**. PanGenie is an order of magnitude **more** accurate on large SVs than on SNPs. The **old GrENE-Net SNP panel** (independent, Beagle-imputed) agrees: SNP nRD 4.5% -> 1.5% -> **1.33%**. **Accuracy does not sag near the dropout holes** (<=0.6 pp), so the two failure modes are separate. Resolves the open question in `RESULTS_LOG.md` 2026-04-29: the LOO outliers are 96-99% concordant vs GrENE-Net and only 81-92% vs "their" assembly — **mislabelled assemblies, not bad genotyping** (7164, 6939, 9606 are new candidates), so these rates are a lower bound. Builders `analysis/panel_qc/pg_dropout/pg_dropout_regions.py` + `analysis/panel_qc/pangenie_accuracy/{bin_loo_concordance.sh,loo_qc_stages.py,loo_qc_stages_grenenet.py,_build_pangenie_accuracy_nb.py}` |



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
