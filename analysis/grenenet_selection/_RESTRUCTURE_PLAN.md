# GEA restructure — proposed file map (REVIEW DRAFT, disposable)

Nothing has been moved. This is the proposal only. Delete this file once the
restructure is done or rejected.

Organising principle: the **three results sections** as stated by the user.

1. **r1_sv_negative_selection** — SVs vs SNPs matched on *initial frequency*;
   more negative selection on SVs; correlated with climate.
2. **r2_gea_nonsnp** — GEA hits that are invisible in SNPs and only appear
   when non-SNP data is used.
3. **r3_persite_gwas** — the same question via per-site GWAS, with the
   ecotype-selection coefficient as the trait.

Everything that is not one of the three is either shared input-building
(`common/`), the LD-block/unit definition machinery (`blocks/`), or method QC
(`qc/`).

---

## Target layout

```
analysis/grenenet_selection/
  README.md
  lib.py                      # stays at top level — 85 scripts import it
  common/                     # shared inputs: AF store, gen/pool matrices, p0, h cache
  blocks/                     # LD-block + unit definition, coherence, tiling
  r1_sv_negative_selection/
  r2_gea_nonsnp/
  r3_persite_gwas/
  qc/                         # panel / seedmix / coverage validation
  archive/                    # existing, unchanged
```

Derived-output dirs (gitignored, not code) stay where they are and are not
touched: `af_store/ gen_matrices/ pool_matrices/ lfmm/ gea/ varexp/ fitness/
site_temporal/ sv_adaptive/ blocks_recompute/ blocks_mcf90*/ sv_snp_ld*/
structure/ ecotype_fitness/ logs/ verify_logs/ rerun_logs/ __pycache__/`.

---

## ⚠ Mechanical blocker — must be handled *during* the move

**86 top-level scripts** contain exactly:

```python
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lib
```

That resolves to the script's own directory. It works today only because
`lib.py` is a sibling. Moving a script one level down silently breaks
`import lib` unless the line becomes:

```python
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
```

which is already the pattern 22 scripts in `phase1_replication/` and
`wza_investigation/` use. So the move is: `git mv` + patch that one line,
per file, in the same commit. Scripted, not hand-edited.

Also present and **safe, verified**: 6 scripts hardcode
`/global/scratch/users/tbellg/kmate/...` — that path is a *symlink* to the
current tree (`/global/scratch/projects/fc_moilab/...`), same inode content,
so it is not a stale second copy. Left alone.

5 scripts use the cwd-relative `sys.path.insert(0, "analysis/grenenet_selection")`,
which only works when run from the repo root. Worth normalising while we are
in there, but it is a pre-existing condition, not caused by the move.

---

## r1_sv_negative_selection/ (31 scripts)

Core frequency-matched comparison:
- `_compute_s_dist_by_stratum.py` — s by class *within initial-frequency strata*, per site
- `_compute_s_climate_slope.py` — the climate-differential views
- `_temporal_selection_snp_vs_nonsnp.py`, `_temporal_s_snp_vs_nonsnp.py`,
  `_temporal_s_plots_snp_vs_nonsnp.py`, `_temporal_s_enrich_initqty.py`,
  `_temporal_sel_drift_maf.py`, `_nonsnp_temporal_category.py`
- `site_variant_temporal_scoef.py`, `_sv_temporal_direct.py`
- `_build_temporal_s_consolidated_nb.py`, `_build_temporal_s_nofilter_nb.py`

SFS / distributional twins:
- `_compute_sfs_shift_by_site.py`, `_build_sfs_shift_nb.py`
- `_compute_sfs_time_site4.py`, `_build_sfs_time_site4_nb.py`

Parallelism / PicMin arm:
- `_compute_parallelism.py`, `_picmin.py`, `_site_parallelism.py`,
  `_plot_site_parallelism.py`, `_build_parallelism_picmin_nb.py`

Per-site enrichment arm:
- `site_sv_enrichment.py`, `aggregate_sites_enrichment.py`,
  `_render_site4_enrichment_fig.py`, `_enrich_threshold_sweep.py`

Audits / artifact checks (these are what make the section defensible — keep with it):
- `_audit_s_classes.py`, `noise_check_sv_s.py`, `_sv_callqual_artifact.py`,
  `_founder_load_test.py`, `_build_sv_selection_audit_nb.py`,
  `_build_sv_polarity_manhattan_nb.py`

## r2_gea_nonsnp/ (32 scripts + 5 existing subdirs)

Models:
- `build_kendall.py`, `_build_kendall_nb.py`, `_kendall_site_collapsed_test.py`
- `build_lfmm_input.py`, `_build_lfmm_nb.py`, `build_lfmm_scree.py`,
  `_build_kselection_nb.py`, `build_struct_pca.py`
- `build_mixedmodel.py`, `_build_mixedmodel_nb.py`
- `build_two_stage_gea.py`, `_build_twostage_nb.py`, `build_two_stage_pooled.py`,
  `build_lmm_scoef.py`

Block aggregation:
- `build_wza.py`, `wza_script.py`, `_build_wza_nb.py`, `build_significant_blocks.py`

Non-SNP-only hit → gene layer (the actual section-2 claim):
- `nonsnp_only_genes.py`, `nonsnp_only_genes_describe.py`,
  `_build_nonsnp_only_genes_nb.py`, `go_enrichment_nonsnp.py`,
  `genes_from_regions.py`, `check_lea_elip_persite.py`,
  `_build_snp_vs_nonsnp_viz_nb.py`,
  `_build_persite_new_peaks_nb.py`, `_build_persite_new_peaks_sv_nb.py`

Locus figures:
- `_build_candidate_nb.py`, `plot_candidates.py`, `plot_phase1style.py`,
  `plot_gi_locus.py`, `plot_micropangenome.py`

Existing subdirs that fold in **as-is** (moved wholesale, contents untouched):
- `phase1_replication/` (90 tracked files — incl. `multiaxis/`)
- `wza_investigation/` (50)
- `cam5_replication/` (50)
- `genes_expl/` (24, **currently untracked** — needs committing)
- `driver_passenger/` (5) — see open question Q1

## r3_persite_gwas/ (41 scripts)

Trait + genotype construction:
- `build_selection_trait.py` — the per-founder per-site selection-coefficient trait
- `founder_genotype.py`, `extract_shortread_geno.py`
- `build_class_grms.py`, `build_untagged_grm.py`

The GWAS itself:
- `class_split_gwas.py`, `_build_class_gwas_persite_nb.py`,
  `_build_class_gwas_multitrait_nb.py`, `plot_class_gwas_pngs.py`
- `ecotype_fitness.py`, `ecotype_fitness_gwas.py`, `ecotype_fitness_enrichment.py`
- `founder_persite_sv_enrichment.py`

Variance partition:
- `varexp_selection.py`, `varexp_untagged.py`, `varexp_bioclim.py`, `lasso_bioclim.py`

SNP-tagging (the "is the non-SNP layer redundant" control):
- `build_nonsnp_tagging.py`, `build_tagging_masked.py`, `_build_tagging_masked_nb.py`,
  `compute_tagging_null.py`, `verify_tagging_permutation_null.py`

Haplotype-unit SV enrichment + its nulls:
- `ecotype_block_sv_enrichment.py`, `ecotype_hap_sv_enrichment.py`,
  `ecotype_hap_sv_regional.py`, `ecotype_hap_sv_rotation.py`,
  `ecotype_hap_sv_rotation2.py`
- `_sv_hap_context.py`, `_sv_hap_freqrobust.py`, `_sv_haplotype_axes_sweep.py`,
  `_sv_hap_rotationnull.py`
- `_sv_founder_direction.py`, `_sv_founder_mechanism.py`, `_sv_passenger_test.py`,
  `_sv_winning_genetics.py`, `_audit_sv_fitness.py`

"Winners" arm:
- `_winners_ratio.py`, `_winners_sv_depletion.py`, `_plot_winner_ratio_allsites.py`,
  `_plot_winner_ratio_vs_purging.py`, `_plot_winners_sv_depletion.py`

## common/ (7 scripts + drivers)

- `build_af_store.py` (+ `build_af_store_array.sh`)
- `build_gen_matrices.py`, `build_gen_matrix.py` (+ `build_gen_matrices_array.sh`)
- `build_pool_matrix.py`
- `build_p0.py`
- `build_sample_h_cache.py`
- `build_hap_membership.py`

`lib.py` deliberately stays at the top level.

## blocks/ (23 scripts + 2 subdirs)

- `blocks_tiling.py`, `recompute_blocks.py`, `dynamic_ld_blocks.py`, `make_coarse_blocks.py`
- `block_cluster_pc1ve.py`, `blockcoherence_data.py`, `block_founder_vs_evolved.py`,
  `block_haplotype_counts.py`, `block_kmer_coverage.py`, `block_missing_sensitivity.py`,
  `block_panel_support_tag.py`, `block_unit_frontier.py`
- `coherence_vs_floor.py`, `delta_p_coherence.py`, `eval_block_coherence.py`,
  `verify_founder_ve.py`, `_check_clq90_blocks.py`
- `unit_distributions.py`, `plot_unit_distributions.py`
- `build_sv_landscape.py`, `plot_sv_landscape.py`
- `_build_blockcoherence_nb.py`, `_build_blocks_units_nb.py`
- subdirs: `hap_blocks/` (23), `bigld_env/` (3)

Note `blocks_tiling.py` is imported as `bt` by other sections — moving it needs
the same import patch on its consumers.

## qc/ (10 scripts + 1 subdir)

- `seedmix_identifiability.py`, `seedmix_kmer_identifiability.py`, `_build_9977_nb.py`
- `compare_window_vs_global_af.py`
- `plot_coverage_distribution.py`, `_build_qc_notebook.py`
- `_compute_panel_overlap_grenenet.py`, `_build_panel_overlap_grenenet_nb.py`
- `_build_panel_stats_nb.py`, `_build_support_nb.py`
- subdir: `seedmix_validation/` (40)

## Misfiled — belongs outside the GEA tree

- `_build_floor_derivation_nb.py` → builds `benchmarks/localonly_p231/FLOOR_DERIVATION.ipynb`.
  Move to `benchmarks/localonly_p231/`.

---

## Open questions (need your call — not assumed)

**Q1. `driver_passenger/`** — is "is the SV the driver or a passenger" part of
section 1 (SV selection) or section 2 (GEA hits)? Provisionally placed in r2.

**Q2. The haplotype-SV-enrichment cluster** (`ecotype_hap_sv_*`, `_sv_hap_*`,
`_winners_*`, ~14 scripts). Conceptually it is "are SVs under selection"
(section 1), but the *unit* is the founder/haplotype GWAS (section 3). I put it
in r3 to keep the unit together. Say if you want it in r1.

**Q3. `notebooks/` (53 tracked files)** — currently one flat dir spanning all
three sections. Split per section, or keep as a single `notebooks/` at top
level? Splitting is more consistent but touches the most files.

**Q4. `results/`, `rerun_kfw_hb/`** — ~~move under `common/`?~~ **Resolved by
size: `rerun_kfw_hb/` is 1.6 TB.** It is the production AF data behind a 2-file
tracked README/map. Leave it exactly where it is; do not move. Same for
`results/` (2.1 MB, keep at top level).

**Q7 (new, raised by the size scan). Code and bulk output are interleaved in
the same tracked dirs.** Measured:

| dir | size | tracked code files |
|---|---|---|
| `phase1_replication/` | **115 GB** | 90 |
| `wza_investigation/` | 952 MB | 50 |
| `driver_passenger/` | 231 MB | 5 |
| `genes_expl/` | 24 MB | 0 (untracked) |
| `notebooks/` | 69 MB | 53 |

`git mv` on one filesystem is a rename, so size is **not** a cost barrier to the
restructure. But it means each of these dirs is a code/output hybrid, and the
outputs are gitignored subpaths *inside* a tracked dir. Optional second axis of
cleanup: push derived outputs to a sibling `out/` tree so code dirs stay small
and browsable. **Not doing this unless asked** — it is a bigger change than the
section regrouping and would touch every hardcoded output path.

Other bulk (all gitignored output, untouched): `af_store/` 52 GB,
`gen_matrices/` 43 GB, `pool_matrices/` 35 GB, `archive/` 50 GB,
`hap_blocks/` 1.2 GB, `site_temporal/` 1.1 GB, `sv_snp_ld/` 2.9 GB.

**Q5. Untracked work.** `genes_expl/` (24 files) plus ~25 loose untracked
notebooks/build scripts are uncommitted. Commit them as part of this, or
handle separately?

**Q6. `_RESTRUCTURE_PLAN.md`** (this file) and the deleted-but-uncommitted
`wza_investigation/RESULTS.md` — confirm the RESULTS.md deletion should be
committed.
