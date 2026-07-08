# Re-run checklist — downstream of the per_founder EM fix

> **STATUS (2026-07-08):** This rerun is **COMPLETE** and has been superseded by the
> full-panel Kf_w rerun. Production outputs now live in
> `results/grenenet_gea/rerun_kfw_hb/{evolved,seedmix}`; the earlier pre-Kf_w
> `rerun_perfounder` directory has been **deleted**. Paths below are updated to
> `rerun_kfw_hb`; keep this doc as the checklist of what depends on what.

**Trigger (2026-07-06).** kMate per-sample outputs are being regenerated under the fixed EM
(`--normalize per_founder --kmer-weight uniform`, filt2inv, **global** mode) — see
`docs/FOUNDER_NORMALIZATION_FIX.md`. This changes **(a)** every sample's per-record SNP/SV
`alt_freq` and **(b)** every sample's founder mixture `h` (and thus the seed-mix `p0` and all
`p0`-anchored selection). Everything downstream that consumes those is stale.

**New per-sample outputs (the new inputs):**
`results/grenenet_gea/rerun_kfw_hb/{seedmix,evolved}/<SAMPLE>.tsv`
(+ per-chrom `<SAMPLE>_Chr{N}.tsv` and `<SAMPLE>_Chr{N}.h_per_chrom.npz`). Global mode → each
run emits BOTH the founder `h` and the per-record AF.

Cohort launch: seed-mix `35549607` (8) + evolved `35549608/609/647` (2168), global + per_founder
+ uniform + filt2inv, count-once DB on node-local `/dev/shm`.

---

## QC & known drops (2026-07-07)
- **Sample QC:** 17 samples with Chr1 usable-panel-k-mer fraction (`nzfrac`) < 0.10 are excluded
  globally (`lib.qc_excluded()`, list in `data/qc_lowcov_exclude.txt`). 5 of them are dead/
  contaminated libraries (nzfrac<0.01 — normal sequencing depth but ~0% panel k-mers ⇒ off-panel
  DNA, NOT low depth; depth is a poor QC signal here, corr(depth,nzfrac)≈0.57). Full audit +
  plots: `results/grenenet_gea/qc_coverage_audit.{csv,ipynb}`.
- **Site 33 dropped (selection trait: 31 → 30 sites).** Site 33 has 7 usable-cohort samples but
  only ONE at gen1 (plot 1, `MLFH330120180607`, nzfrac 0.088); its other 6 are gen2-only in plots
  with no gen1 anchor. The selection slope needs a gen1 anchor within a plot, so site 33's whole
  estimate rests on that single low-usable-data sample — which the nzfrac<0.10 QC excludes → site
  33 has no anchorable plot → dropped. Site 33 has long been anomalous (it is the single
  lowest-mean-nzfrac site, 0.170), so dropping it is accepted/expected, not a concern. (A laxer
  nzfrac<0.05 cut would retain site 33, but only on that one borderline sample.)
- **p0 floor:** `build_selection_trait.py` floors the estimated seed-mix p0 at 1e-4 so no ecotype
  disappears (keeps all 231 analyzable, NOT forced to 1/231). Only founder 9977 hits the floor
  (p0~4e-6 under uniform ω; it is diffusely non-identifiable — k-mer-poor, no private k-mers, max
  Jaccard 0.68 ≈ panel median — so it carries ~no signal and reads s≈0).

## Step 0 — repoint + clear caches  ✅ DONE (2026-07-06)
- Repointed `analysis/grenenet_gea/lib.py` `OUT`/`SEEDMIX`, `build_af_store.py` `OUT`, and
  `_build_support_nb.py` `OUTBASE` → `results/grenenet_gea/rerun_kfw_hb/{evolved,seedmix}`.
  (Old multinomial outputs remain at `results/grenenet_kmate_arch3` / `seedmix_kmate_arch3`.)
- Moved stale on-existence caches aside (suffix `preFix_multinomial`, reversible — they would
  otherwise silently return old data): `af_store/`, `group_means.npz`, `p0_seedmix_all.pkl`,
  `pilot_qc.csv`, `ecotype_fitness/sample_global_h.npz`, `varexp/selection_s_matrix.npz`.

> Do NOT start Groups A–C until the cohort run finishes (watcher tallies 8 + 2168 final TSVs).

> **⚠ ORDERING CORRECTION (2026-07-07):** `build_af_store.py` (Group B step 1) must run **FIRST**,
> before Group A — `ecotype_fitness.py` calls `lib.pool_table()`, which reads
> `af_store/samples.npy`. So the true start is: **build_af_store.py → ecotype_fitness.py →
> build_selection_trait.py**, then the rest of B, then C. (`build_selection_trait.py` FLOOR
> already revised to 2e-4 → keeps 230 founders, estimated p0, no 1/231 forcing.)

---

## GROUP A — founder h / p0 / selection (MOST affected; the decomposition itself moved)
Run in order:
1. `ecotype_fitness.py` → `ecotype_fitness/{sample_global_h.npz, ecotype_fitness.csv, founder_site_dh.npz}` — builds the `sample_global_h.npz` cache the trait reads. **Run first.**
2. `seedmix_identifiability.py` + notebook `notebooks/seedmix_kmate_vs_hapfire.ipynb` — **validation:** confirm the new seed-mix p0 (more founders called, ~1/231, fewer collapsed) before trusting anything downstream.
3. `build_selection_trait.py` → `varexp/selection_s_matrix.npz`. **⚠ FIX DURING RERUN:** the `FLOOR = 1e-3` / `analyzable = p0 > FLOOR` cutoff currently drops to **212 founders**. The fix makes far fewer founders collapse, so this now discards legitimate rare-start founders — revise to keep ~230 / anchor p0 at 1/231. Must run AFTER (1), BEFORE any varexp/class-GWAS.
4. `varexp_selection.py`, `varexp_untagged.py` (read selection_s_matrix + class_grms).
5. `class_split_gwas.py` → `varexp/class_gwas_{snp,nonsnp,sv}.npz`. **Doubly stale:** trait changed AND genomic control was removed this session (uncommitted).

**NOT stale (don't rerun):** `build_class_grms.py`, `build_untagged_grm.py` (panel-derived GRMs);
`varexp_bioclim.py`, `lasso_bioclim.py` (climate phenotype + panel GRMs — no kMate AF/selection).

---

## GROUP B — per-record AF (SNP/SV-level GEA/GWAS/WZA). Dependency-ordered:
1. `build_af_store.py` → `af_store/` (rebuild from new TSVs; everything in B depends on it).
2. `build_p0.py` → `af_store/p0_{snp,nonsnp}.npy` (Δp reference).
3. `build_gen_matrices.py`, `build_pool_matrix.py`, `build_struct_pca.py`.
4. `build_kendall.py` → `build_two_stage_gea.py` → `build_wza.py` → `build_significant_blocks.py`.
5. `build_lmm_scoef.py`, `build_mixedmodel.py`, `build_lfmm_input.py` (+LFMM run).
6. `baypass_build_inputs.py`.
7. SV-adaptive / temporal family that reads `p0_nonsnp.npy` + pool matrices (rerun the ones in use):
   `_sv_passenger_test.py`, `_sv_founder_direction.py`, `_sv_founder_mechanism.py`,
   `_temporal_*`, `site_variant_temporal_scoef.py`, `_compute_s_climate_slope.py`,
   `_compute_parallelism.py`, `nonsnp_only_genes*.py`, etc.

---

## GROUP C — figures / tables / notebooks (regenerate after their A/B inputs)
`plot_class_gwas_pngs.py`; `varexp/manhattan_*.png`, `persite_lambda_3way.png`;
`plot_candidates.py`, `plot_varlen_manhattan.py`, `top_block_genes.py`, `genes_from_regions.py`;
`lib.build_qc_table` (`pilot_qc.csv`) + `01_pilot_qc.ipynb`; `lib.build_group_means`
(`group_means.npz`) + `02_*`; notebooks `03_kendall`…`10_lfmm_k_selection`, the `s_*` family,
`class_gwas_*`, `sv_enrichment`, `sv_selection_currency` (via their `_build_*_nb.py` generators);
refresh the `EXPORT_MANIFEST.md` / `push_gea_to_drive.sh` export set.

---

## GROUP D — window-mode chain  🅧 STALE, NOT BEING UPDATED (decision 2026-07-07)

> **DECISION (2026-07-07):** window mode is **not in production use** (the GrENE-Net cohort is
> analyzed in GLOBAL mode — heavy selfing). The entire window-mode chain below is therefore left
> **AS-IS and is STALE**: its inputs (the separate `grenenet_kmate_window` / `_seedmix` run,
> per-block `*_Chr{N}.h_blocks_per_chrom.npz`) were produced by the **old multinomial EM** and are
> NOT being regenerated under the per_founder fix. Any output derived from it —
> `hapfreq/`, `founder_gwas_multisite.py` (+ its `cross_site_winners*`, `derive_climate_axis` PC1,
> `multisite_climate_perm`, `founder_persite_sv_enrichment`), the whole `block_ld_lmm*` family, and
> all of `sv_adaptive/` — is **stale and should not be trusted / used** until (if ever) window mode
> is revisited. Do NOT spend effort rerunning it. If it is ever revived, it needs its own
> window-mode cohort rerun (per_founder default) + repointing its hardcoded `grenenet_kmate_window*`
> paths. (`class_split_gwas.py` / `founder_gwas_multisite.py` also had genomic control removed this
> session — for the GLOBAL path that lives in Group A and IS rerun; the window-mode
> `founder_gwas_multisite.py` is part of this stale chain and is not.)
>
> **Archived (2026-07-07):** the raw `grenenet_kmate_window` / `_seedmix` / `_smoke` result
> directories (~1.6 TB) have been moved to `results/archive/` given this decision — any script in
> this stale chain that still points at `results/grenenet_kmate_window*` will need repointing to
> the archived location if this chain is ever revived.

<details><summary>Stale window-mode chain (archived detail — not being rerun)</summary>
The global rerun does **NOT** regenerate the separate **window-mode** run
(`grenenet_kmate_window` / `_seedmix`, per-block `*_Chr{N}.h_blocks_per_chrom.npz`).
- **Regardless of scope:** `founder_gwas_multisite.py` is stale NOW from the uncommitted
  genomic-control removal → must be regenerated (its `_meta.json` λ_GC self-check will fire
  otherwise). Its consumers follow: `cross_site_winners*`, `derive_climate_axis.py` (PC1 →
  regenerate first; feeds the clq90_pc1 run), `multisite_climate_perm.py`,
  `founder_persite_sv_enrichment.py`, and all of `sv_adaptive/`.
- **DECISION:** is the window-mode founder-GWAS / `block_ld_lmm*` / `sv_adaptive` family still a
  live analysis? If **yes**, it needs its own **window rerun** (per_founder default; then rebuild
  `build_hapfreq_matrix.py`, `founder_gwas_231.py`, `ecotype_selection_site.py`, the
  `block_ld_lmm*` family, etc.), and its many hardcoded `grenenet_kmate_window*` paths must be
  repointed. If **superseded** by the global founder-level path, Group D is stale only from the
  GC removal (regenerate `founder_gwas_multisite.py` + consumers, no window EM refit).

</details>

---

## Hard ordering
`0 → 1(af_store) → 2(p0)/3(ecotype_fitness) → 5(build_selection_trait, with floor fix)
→ 7A(varexp/class-GWAS)`; `ecotype_fitness` before `build_selection_trait`;
`build_selection_trait` before every varexp/class-GWAS; `build_pool_matrix`/`build_p0` before
Group B stage 4+; `derive_climate_axis` PC1 before the clq90_pc1 founder-GWAS in Group D.
