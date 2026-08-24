# Downstream re-run map — after the rerun_kfw_hb cohort (2026-07-07)

The founder mixtures `h` changed (full-panel Kf_w fix; evolved compositions are
more concentrated, min founding p0 9977: 4e-6 → 1.5e-4). Everything derived from
`h` or from the per-record AF store must be regenerated. Ordered by dependency.

Roots to repoint (in `analysis/grenenet_selection/lib.py`):
- `OUT`     `rerun_perfounder/evolved` → `rerun_kfw_hb/evolved`
- `SEEDMIX` `rerun_perfounder/seedmix` → `rerun_kfw_hb/seedmix`
(`GEA` and `AF_STORE` unchanged; AF_STORE is *rebuilt in place*.)

---

## Step 0 — INVALIDATE stale caches (critical: these silently serve old h/AF)
Delete before re-running, else scripts reuse pre-fix values:
- `ecotype_fitness/sample_global_h.npz`  (ecotype_fitness.py cache — the big trap)
- `fitness/sample_genome_h.npz`           (build_sample_h_cache.py)
- `group_means.npz`                       (lib.load_group_means cache)
- `af_store/sv_support_cache.npz`         (SV support cache)
- `gea/candidate_block_cache.npz`         (candidate notebook cache)
- `varexp/selection_s_matrix.npz`         (regenerated in Branch A)
- (existing `*.preFix_multinomial` snapshots are a *previous* fix's backups — leave or archive separately)

## Step 1 — REBUILD roots (feed both branches)
| script | output | notes |
|---|---|---|
| `build_af_store.py` (init + convert, array `build_af_store_array.sh`) | `af_store/` | reads the 2168+8 new per-sample TSVs; SNP/non-SNP split store. HEAVY → sbatch. |
| `build_p0.py` | `p0_seedmix_all.pkl` (+ founding AF) | founding p0 from the new seedmix |
| `rebuild_group_means.sbatch` | `group_means.npz` | per-group Δp foundation |

---

## Branch A — ECOTYPE SELECTION COEFFICIENT (h-dependent)
Order matters; each reads the previous.
1. `build_sample_h_cache.py` → `fitness/sample_genome_h.npz`  (per-sample genome h)
2. `ecotype_fitness.py` → `ecotype_fitness/{ecotype_fitness.csv, founder_site_dh.npz, sample_global_h.npz}`  (per-founder fitness axes + the h cache)
3. `build_selection_trait.py` → **`varexp/selection_s_matrix.npz`**  ← the per-founder×site selection coefficient S (the headline "selection on ecotypes"). Note: `P0_FLOOR=1e-4` is now a no-op (9977 = 1.5e-4).
4. `build_fitness_table.py` → fitness table (census/relative)
5. Founder GWAS on the traits:
   - `ecotype_selection_site.py` (per-site), `founder_gwas_231.py` (single-site)
   - `founder_gwas_multisite.py` → `hapfreq/multisite_founder_gwas.{csv,npz}` (JOINT/GLOBAL/CLIMATE)
   - `class_split_gwas.py` → `varexp/class_gwas_{snp,nonsnp,sv}.npz` (reads selection_s_matrix)
6. Variance partition: `varexp_selection.py`, `varexp_untagged.py` (read selection_s_matrix)
7. GWAS downstream: `run_multisite_downstream.sh`, `multisite_climate_perm.py`,
   `ecotype_fitness_gwas.py`, `ecotype_hap_sv_enrichment.py`, `sv_adaptive/*`,
   `founder_persite_sv_enrichment.py`, `derive_climate_axis.py`

## Branch B — AF-DEPENDENT GEA (AF-store-dependent)
Driven off the rebuilt `af_store/` (Step 1). Key drivers (each has many followers):
1. `build_gen_matrix.py` / `build_gen_matrices.py` (+ array) → per-generation AF matrices  ← feeds most of B
2. `build_pool_matrix.py`, `build_struct_pca.py`, `build_significant_blocks.py`
3. Per-marker temporal/climate stats: `build_kendall.py`, `build_lmm_scoef.py`,
   `build_mixedmodel.py`, `site_variant_temporal_scoef.py`
4. Two-stage SV climate-GEA: `build_two_stage_gea.py`, `build_two_stage_pooled.py`
5. WZA: `build_wza.py` → `gea/wza/wza_*.csv`
6. LFMM: `build_lfmm_input.py` + `phase1_replication/run_lfmm_gea.sh`
7. Candidates: `_build_candidate_nb.py` (delete candidate_block_cache.npz first)
8. New-panel GEA suite (RETIRED 2026-07-21 into
   `archive/gea_newpanel_snp_nonsnp_fork_retired/`; unique analyses preserved
   there): RDA `rda_prda/rda_pcreg`, perm-nulls `permnull_site_gea`/`permnull_precip_gea`,
   `shape_gea*`, `ksweep_hits`. The site-level PC1/bio-axis LFMM scan was ported
   into `phase1_replication/axis_scan/`; the production quasi-binomial into
   `phase1_replication/run_quasibinom.py`.
9. SV / temporal analyses: `_sv_*`, `_temporal_*`, `noise_check_sv_s.py`,
   `_nonsnp_temporal_*`, `_compute_s_climate_slope.py`
10. Prediction / winners: `predict_ecotype_performance.py`, `cross_site_winners*.py`,
    `compare_analyses_site.py`, `cross_chrom_agreement.py`

## Branch C — NOTEBOOKS / FIGURES (rebuild last, after A+B)
`_build_*_nb.py`: multisite_gwas, class_gwas_multitrait, class_gwas_persite,
candidate, wza, support, twostage, svenrich, sv_selection,
climate_locality, ccagree — plus the `FINAL_RESULTS_*` benchmark notebooks.

Retired 2026-07-08 (see `analysis/grenenet_selection/archive/gwas_notebooks_retired/`):
`class_split_gwas` (its peak-overlap question folded into the *end* of class_gwas_multitrait
+ class_gwas_persite as two Bonferroni/FDR block-overlap tables, 3-way snp/nonsnp/sv), `founder_gwas231`, and the `multisite_founder_gwas_clq{50,90,90_pc1}`
block-clustering variants. Rebuild tooling kept (`_build_foundergwas_nb.py`,
parameterized `_build_multisite_gwas_nb.py` / `run_multisite_downstream.sh`).

---

## Suggested minimal critical path (headline results first)
Step 0 (invalidate) → Step 1 (`af_store`, `p0`, `group_means`) →
A1–A3 (**selection coefficient** `selection_s_matrix.npz`) → A5 (`founder_gwas_multisite`, `class_split_gwas`) →
B1 (`gen_matrix`) → B4/B5 (two-stage + WZA) → C (the two headline notebooks).
Everything else (gea_newpanel sweep, the many `_sv_*`/`_temporal_*` explorations)
is secondary and can follow.

## Compute notes
- `build_af_store` and `build_gen_matrix` are the heavy ones → `sbatch` (see the
  `*_array.sh` wrappers). Selection/GWAS steps are light (run off caches).
- Sanity gate before trusting Branch A/B: diff a few `selection_s_matrix` columns
  and a `gen_matrix` slice new-vs-stale to confirm the shift is the expected
  Kf_w concentration, not a plumbing break.
