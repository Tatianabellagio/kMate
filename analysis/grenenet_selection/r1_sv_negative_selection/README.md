# Result 1 — negative selection on SVs

**The question.** Comparing SVs and SNPs *at the same initial frequency*, is
there more negative selection on SVs, and does it track climate?

Everything here is per-variant and temporal: no haploblock unit, no GWAS. The
estimand is a selection coefficient `s` from the allele-frequency trajectory
across generations, with replicate plots within a site as the replication.

Outputs → `results/` (`sv_adaptive/`, `site_temporal/`); figures → `results/*/plots/`.

**Read the result in:** `../notebooks/temporal_s_consolidated.ipynb` — the
consolidated notebook (2026-07-15) that replaced seven separately-named ones.
Supporting: `parallelism_picmin.ipynb`, `sfs_shift_by_site.ipynb`,
`sfs_time_site4.ipynb`, `temporal_s_nofilter.ipynb`.

Narrative + caveats: `SV_TEMPORAL_PURGING_SUMMARY.md` (this folder).

> ⚠ **Do not retire the `_temporal_s_*` / `_compute_s_*` family on name
> similarity.** They look like successive versions but were audited 2026-07-15
> and each tests a genuinely different statistic. The docs around them contain
> *retractions* of a superseded label, not declarations of one.

---

## A. The `s` estimator and the SNP-vs-SV comparison

| script | what it does | writes |
|---|---|---|
| `_temporal_selection_snp_vs_nonsnp.py` | the origin: pure per-variant, per-site temporal test | `temporal_selection_snp_vs_nonsnp_summary.csv` |
| `_temporal_s_snp_vs_nonsnp.py` | switches Δp → selection coefficient `s`, multi-generation | `temporal_s_snp_vs_nonsnp.csv` |
| `_temporal_s_plots_snp_vs_nonsnp.py` | `s` with **plots as replicates** — the power-gaining refinement | `temporal_s_plots_snp_vs_nonsnp.csv` |
| `_compute_s_dist_by_stratum.py` | **the headline computation**: distribution of `s` by class *within initial-frequency strata*, per site | `s_dist_by_stratum.csv`, `..._sitemeta.csv` |
| `_temporal_s_enrich_initqty.py` | SNP-null-free variant: matches on initial quantity instead of a SNP null | `temporal_s_enrich_initqty.csv` |
| `_nonsnp_temporal_category.py` | asks the same of non-SNP as a whole category, not SVs alone | — |
| `_temporal_sel_drift_maf.py` | robustness: explicit drift and MAF controls | `temporal_sel_drift_maf.csv` |
| `_build_temporal_s_consolidated_nb.py` | builds the consolidated notebook + its figures | `temporal_s_consolidated.ipynb` |

## B. Climate

| script | what it does | writes |
|---|---|---|
| `_compute_s_climate_slope.py` | the two climate-differential views, constructed to cancel the founder confound | `s_climate_slope*.npz`, `s_climate_slope_sign_by_site.csv` |

## C. Artifact controls — what makes section A defensible

Each answers a specific "could this be an artifact?" challenge. Keep them with
the result; they are the reason it survives.

| script | the challenge |
|---|---|
| `_audit_s_classes.py` | are the SNP/indel/SV `s`-distributions matching because of a bug, or really? |
| `noise_check_sv_s.py` | is the SV > SNP \|s\| excess a k-mer-support (noise) artifact? |
| `_founder_load_test.py` | founder-level test of the "selection on SVs" conclusion (audit 2026-07-06) |
| `_sv_temporal_direct.py` | do common SVs actually fall in frequency over generations? |
| `_sv_callqual_artifact.py` | calling-confidence proxy for the insertion-polarity artifact-vs-biology question |
| `_extract_vcf_callqual.sh` | pulls the `F_MISSING` / `CONFLICT` / `MA` fields the above needs |

## D. Parallelism / PicMin arm

Does the same signal appear as parallel change across replicate plots and sites?

| script | what it does | writes |
|---|---|---|
| `_compute_parallelism.py` | per-variant parallelism across replicate plots (AF-vapeR / PicMin philosophy) | `parallelism.npz` |
| `_picmin.py` | PicMin proper (Booker et al. 2024), lineages = sites | `picmin.npz` |
| `_site_parallelism.py` | per-site parallelism of *founder* frequency change | `site_parallelism.csv` |
| `_plot_site_parallelism.py` | per-site purging vs per-site parallelism — the confound scatter | `site_parallelism_vs_purging.png` |
| `_build_parallelism_picmin_nb.py` | builds the consolidated replicate-arm notebook | `parallelism_picmin.ipynb` |

## E. Site-frequency-spectrum arm

A different signature from the central-tendency test: the SFS *mean* shift is
null, but SV **extinction rate** is elevated.

| script | what it does | writes |
|---|---|---|
| `_compute_sfs_shift_by_site.py` | folded-SFS shift per site, per class, per p0 decile | `sfs_shift_by_site.csv`, `..._sitemeta.csv` |
| `_build_sfs_shift_nb.py` | builds the SFS-shift notebook | `sfs_shift_by_site.ipynb` |
| `_compute_sfs_time_site4.py` | AF spectrum over generations 0→3 at site 4, by class | `sfs_time_site4.npz`, `..._summary.csv` |
| `_build_sfs_time_site4_nb.py` | builds the site-4 SFS-over-time notebook | `sfs_time_site4.ipynb` |

## F. Per-site block enrichment (site-4 pilot → cross-site)

The one block-unit arm in r1: are temporally-selected blocks SV-enriched?

| script | what it does | writes |
|---|---|---|
| `site_variant_temporal_scoef.py` | per-variant `s` at one site, by class | `site{N}_scoef_*.npz` |
| `site_sv_enrichment.py` | SV enrichment in selected clq0.9 blocks at one site | `site{N}_sv_enrichment.json`, `site{N}_clq90_blocks.csv.gz` |
| `run_site_enrichment.sbatch` | SLURM driver for the above across sites | — |
| `aggregate_sites_enrichment.py` | cross-site synthesis | `cross_site_enrichment.csv/.png` |
| `_enrich_threshold_sweep.py` | is the enrichment an artifact of the top-fraction cut? size-matched permutation null | stdout only |
| `_render_site4_enrichment_fig.py` | renders the site-4 figure from precomputed plot data | `site4_sv_enrichment.png` |

---

## Filed here but arguably belonging elsewhere

Both are **builders whose notebooks are about a different question** — flagged,
not moved, pending a decision:

- `_build_sv_selection_audit_nb.py` → `sv_selection_haplotype_audit.ipynb`. This
  is the haplotype-unit SV-enrichment audit; its compute scripts live in
  `../extras/`. Only the builder stayed behind.
- `_build_sv_polarity_manhattan_nb.py` → `sv_polarity_enrichment_clq90.ipynb`.
  Content is a WZA/clq0.9 GEA-block enrichment (Kendall + LFMM K=16 + binomial
  vs bio1), i.e. r2 material.

## Also here

- `BAYPASS_TEMPORAL_PLAN.md` — a **plan, not completed**; its own banner marks
  its haploblock-frequency inputs as stale (the retired Pipeline-B chain).
