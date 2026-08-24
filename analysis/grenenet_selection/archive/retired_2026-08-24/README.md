# Retired 2026-08-24

Derived output retired during the analysis-tree cleanup. Everything here was
**moved, not deleted** — restore by moving it back one level up.

Each item was checked for inbound references across the live code dirs before
being moved. Nothing here had a reader.

## Directories

| dir | why retired |
|---|---|
| `blocks_mcf90_p80/` | The p80-control BigLD run never finished: the directory holds only two `.btmp` scratch files (`chr1_clq0.9_40_geno_matrix.btmp`, `chr1_clq0.9_40_snpINFO.btmp`) and no `.tsv` block map. Zero references in any `.py`/`.sh`/`.ipynb`. Reproducible by re-running `blocks/blocks_recompute_p80.sbatch`. |
| `verify_logs/` | 7 ad-hoc logs from the 2026-07-15 parallelism/PicMin verification session. No programmatic reader. The findings they verified are written up in `r1_sv_negative_selection/SV_TEMPORAL_PURGING_SUMMARY.md`. |
| `rerun_logs/` | 37 logs from the 2026-07-08 Kf_w rerun. No programmatic reader. The rerun itself is documented in `common/rerun_kfw_hb/RERUN_MAP.md`. |
| `structure/` | Single file `struct_pca_gen3.npz` (pool PCA scores/eigvals) written by `r2_gea_nonsnp/build_struct_pca.py`. A repo-wide grep for `struct_pca` found the writer and two inventory mentions, but **no reader**. ⚠ `common/rerun_kfw_hb/RERUN_MAP.md` still lists `build_struct_pca.py` as a step-2 rebuild target — if the PCA-covariate arm is revived, restore this. |

## `sidecars/`

Pre-correction backup snapshots (~198 MB), superseded by the live files that
sat beside them. `r3_persite_gwas/VAREXP_SELECTION_HANDOFF.md` records that the
current versions were regenerated and verified on 2026-07-21.

| file | superseded by |
|---|---|
| `gwas_z.npz.preRerun0721` (166 MB) | `ecotype_fitness/gwas/gwas_z.npz` |
| `sample_global_h.npz.preKfw`, `.preFix_multinomial` | `ecotype_fitness/sample_global_h.npz` |
| `sample_genome_h.npz.preKfw`, `.preRerun0721` | `fitness/sample_genome_h.npz` |
| `selection_s_matrix.npz.preKfw`, `.preFix_multinomial` | `varexp/selection_s_matrix.npz` |
| `candidate_block_cache.npz.preKfw` | `gea/candidate_block_cache.npz` |
| `cross_site_enrichment.mac2.csv.bak` | `site_temporal/cross_site_enrichment.csv` |

## Deliberately NOT retired

Two directories look stale but are load-bearing:

- **`blocks_recompute/`** — superseded as the *production* block map (`blocks/diag_genomewide.sbatch` calls it "the OLD 0.5 map"), but `blocks/_build_blockcoherence_nb.py` still globs its `chr*_panel_support_tag.csv` and `chr*_missing_sensitivity.csv` for two panels of `notebooks/block_coherence_clqcut.ipynb`.
- **`sv_snp_ld/`** — the name collides with the retired `archive/sv_snp_ld_unmasked_superseded/`, but this dir holds *inputs* (the GrENE-Net VCF + extracted short-read genotypes), read by `r3_persite_gwas/build_tagging_masked.py`. The superseded tagging *outputs* were archived separately, earlier.
