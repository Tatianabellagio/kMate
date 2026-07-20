# Archived: window-unit / haploblock-frequency (hapfreq) diagnostics — RETIRED 2026-07-08

## What this is
Scripts from the retired **window-unit** and **haploblock-frequency (hapfreq)**
branch of the GrENE-net GEA. They estimated per-*window* founder frequencies and
projected them into per-haplotype-cluster frequencies, plus a family of
window-vs-global / cross-chromosome agreement diagnostics.

## Why it's retired
1. Production kMate now runs `--unit chrom` (GLOBAL mode; see
   `../../GLOBAL_MODE_DECISION.md`). The `--unit chrom` cohort produces **no
   per-block `h_blocks`** files.
2. These scripts read the **window stores**
   `results/grenenet_kmate_window` and `results/grenenet_kmate_window_seedmix`
   (per-sample `*_Chr{N}.h_blocks_per_chrom.npz`). **Those directories have been
   deleted from disk**, so every script here is broken/unrunnable as-is.
3. The window-vs-chrom question is settled: at r²=0.1 with eps=0 the haploblocks
   collapse to ~231 genome-wide (≈ chromosome-wise), so the window/haploblock
   reframe is "framing only" — it does not change the AF/`h` estimates. Earlier
   agreement reads (r≈0.99–0.995) were consistent with this but predate the
   full-panel `Kf_w` normalization fix and were partly circular (see
   `../../WINDOW_UNIT_VALIDATION.md`, `../../GLOBAL_MODE_DECISION.md`).

## Contents
- `cross_chrom_agreement.py`, `cross_chrom_by_generation.py` — cohort cross-chromosome
  `h`-agreement diagnostics (read the deleted window store).
- `block_residual_consistency.py`, `block_vs_global_h.py` — per-site block-vs-global
  `h` consistency diagnostics.
- `compare_analyses_site.py` — per-site window-vs-global comparison.
- `compare_win_vs_glob.sbatch` — driver (also referenced the old
  `results/grenenet_kmate_arch3` store, which is NOT retired and stays on disk).
- `gen9_window/` — the full window/gen9 "breakage" diagnostic script subtree
  (clustering, extraction, seedmix, haploblock/dynld founder-vs-evolved, notebook
  builders). Self-contained: no code outside this subtree imported it.
  **NOTE:** the associated DATA directory `results/grenenet_gea/gen9_window/`
  (kept for GEA) was **NOT** moved — only the scripts are archived here.

## Update 2026-07-08 (later same day) — Pipeline-B chain also retired, separately
The "NOT archived" chain originally described here — `build_hapfreq_matrix.py`,
`build_hapfreq_p0_seedmix.py`, `build_hap_trajectories.py`, `build_hap_gea.py`,
`build_hap_wza.py`, its drivers, `phase1_replication/build_hap_lastgen_matrix.py`, and
the dedicated consumers (`sv_adaptive/sv_enrichment_gea.py`, `sv_temporal_markermatched.py`,
etc.) — was itself retired (both SV threads it fed had resolved null) to
**`../pipelineB_hapfreq_retired/`** (a sibling archive dir, not this one). `build_hap_membership`
was kept in place (shared with the live SV-selection audit at the time).

## Update 2026-07-10 — the remaining maintainer decision made: fully retired
The "maintainer decision" flagged below as pending is resolved: the window-mode founder-GWAS /
`block_ld_lmm*` family is **not** a live analysis. Consolidated into this directory as part of
the `analysis/` vs `results/` centralization pass:
- `hapfreq/` — the DATA directory this README's original note said was "kept for GEA" (moved
  from `results/grenenet_gea/hapfreq/`, 1.8 GB). With the Pipeline-B chain above already retired,
  nothing live still reads it.
- `gen9_window/` — the matching DATA for the script subtree already archived here (moved from
  `results/grenenet_gea/gen9_window/`, 14 GB).
- `window_vs_global/` — data + its driver `aggregate_window_vs_global.py` (moved from
  `results/grenenet_gea/window_vs_global/` + `analysis/grenenet_gea/`, 9.3 MB).

**Correction (later same pass, 2026-07-10):** `sv_adaptive/` was initially archived wholesale too
— that was wrong for the *data*, half-right for the *code*. `analysis/grenenet_gea/sv_adaptive/`
was genuinely dead as a **code folder** (12 scripts, all confirmed reading `hapfreq/` or
`lib.multisite_gwas_raw()` — `audit_sv_passenger_emmax.py`, `build_sv_landscape.py`,
`climate_cluster_enrichment.py`, `locality_index.py`, `plot_sv_block_climate.py`,
`plot_sv_landscape.py`, `sv_block_drilldown.py`, `sv_block_haplotype_resolve.py`,
`sv_enrichment.py`, `sv_frac_markermatched.py`, `sv_joint_diagnosis.py`, `sv_joint_mechanism.py`
— correctly archived here, in `sv_adaptive/`). But its **output directory name** (`sv_adaptive/`,
reached via `lib.GEA`) was *also* the write target of a completely different, **live** family of
~25 top-level scripts (the s_\*/parallelism/picmin/temporal_\* SV-selection-vs-drift analysis,
`docs/RERUN_AFTER_FIX.md` Group B item 7) that share no code with the dead chain — they just
happened to write into a directory with the same name. Archiving the whole data directory broke
that live family's output path. Fixed: the live ~43 files were moved back out to a **recreated
live** `analysis/grenenet_gea/sv_adaptive/` (not this one), and the ~25 live scripts' path
references were reverted to point there. Only the genuinely-dead output files stayed here.

The consumer scripts of `hapfreq/` — `founder_gwas_multisite.py`, `cross_site_winners.py`,
`cross_site_winners_multisite.py`, `derive_climate_axis.py`, `multisite_climate_perm.py`,
`founder_gwas_231.py`, `ecotype_selection_site.py`, `build_fitness_table.py`,
`predict_ecotype_performance.py`, `winner_convergence_site.py`, `build_crosssite_climate.py`,
`omega_contrast.py`, `baypass_build_inputs.py`, `_sv_haplotype_enrichment.py`, `_sv_ld_check.py`,
`_sv_haplotype_check.py`, their notebook builders (`_build_foundergwas_nb.py`,
`_build_multisite_gwas_nb.py`, `_build_ccagree_nb.py`, `_build_predict_perf_nb.py`,
`_build_climate_locality_nb.py`), their generated notebooks (`multisite_founder_gwas.ipynb`,
`predict_ecotype_performance.ipynb`, `cross_chrom_agreement.ipynb`,
`climate_cluster_locality.ipynb`), and their sbatch/sh drivers (`run_cross_site.sbatch`,
`predict_ecotype_perf.sbatch`, `run_multisite_gwas.sbatch`, `run_multisite_downstream.sh`) —
all individually verified reading `hapfreq/` or `lib.multisite_gwas_raw()`, none importable
by live code — were moved here too, physically, not just repointed.

One shared helper, `genome_h()` (+ `CHROMS`), was pulled out of `ecotype_selection_site.py`
before archiving it: it's a generic per-sample founder-h loader with no `hapfreq` dependency,
and two live scripts (`_sv_winning_genetics.py`, `build_sample_h_cache.py`) import it. Promoted
to `lib.py` instead of leaving a shim behind.

**Left for a maintainer, not archived:** `_sv_founder_direction.py` and `_sv_founder_mechanism.py`
have a genuinely mixed dependency — live `af_store`/`pool_matrices` data AND the dead
`hapfreq/multisite_founder_gwas_clq90_pc1` as their candidate-locus input. They can't run as-is,
but whether that's "retire" or "needs a non-hapfreq locus source" is a scientific call, not a
filing one. Their two output CSVs stayed in `sv_adaptive/results/` here (not moved to the live
dir) pending that call.

**Second correction (same pass):** `build_sv_landscape.py` + `plot_sv_landscape.py` were
archived here too on first pass, but checked individually they read only panel block TSVs
(`{ch}_clq0.9_blocks_clq0.9.tsv`) — no `hapfreq`/`multisite_gwas_raw` dependency at all. Their
output `sv_landscape_clq0.9.csv` is a shared input consumed by *both* the genuinely-dead
scripts above *and* a separate live SV-haplotype audit (`_sv_hap_context.py`,
`_sv_hap_rotationnull.py`, `_sv_hap_freqrobust.py`, `_sv_haplotype_axes_sweep.py`,
`_build_sv_selection_audit_nb.py` — deliberately kept per
`../pipelineB_hapfreq_retired/README.md`, not part of this retirement). Reverted both scripts
+ their data back to live `analysis/grenenet_gea/sv_adaptive/`. See that same README for an
open question about whether the audit itself is still fully runnable — 3 of its 4 producer
scripts also call `lib.multisite_gwas_raw()` and read the dead Pipeline-B `hap_gea.csv`.

Date archived: 2026-07-08 (first script batch); 2026-07-10 (remaining data, corrected sv_adaptive
split, hapfreq/multisite-GWAS consumer family, final retirement decision).
