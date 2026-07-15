# GrENE-Net SV-GEA analysis (kMate)

Re-running the GrENE-Net phase-1 genotype–environment-association (GEA) pipeline on the
**structural-variant + indel** allele-frequency time-series produced by kMate, to ask whether
large SVs (and the non-SNP layer generally) carry adaptive signal that SNP studies miss.

> **Index updated 2026-07-06.** Findings summary below reflects the current state; the
> per-doc index points to the standing decisions, results, and plans. Superseded session
> logs were consolidated into the surviving docs (git `5cefcfa` snapshots the pre-cleanup state).

## Current findings (short version)

The answer depends on **which unit and which question**, and they do not all agree — state
them separately, don't collapse to one headline:

1. **Per-variant temporal selection — split verdict, audited 2026-07-15 (do not cite as a uniform
   null; the original "these all share one confound" dismissal was checked directly and does not
   hold for most of these arms).**
   Two genuinely different questions get asked under this heading, and they have different answers:
   - **Whole-distribution / median shift: verified null.** The de-trended metric
     (`s_distribution_by_site`, subtracting the per-p0-bin ALL-class median) shows **SV
     excess-vs-baseline ≈0: median +0.0011, negative at only 15/31 sites, sign test n.s.** No
     genome-wide, frequency-independent shift in typical SV behavior.
   - **Tail-specific, hot-site-concentrated purging: verified REAL, survives the same de-trending.**
     `s_histogram`, `s_vs_climate`, `ecdf`, `ecdf-difference`, `shiftfunction`, and climate-slope β
     all independently describe the *same* pattern in their own takeaways: not a whole-distribution
     shift, but extra SV mass in the purged tail, concentrated at hot sites. This was previously
     dismissed as sharing the median-shift's confound — checked directly (10th-percentile SV vs
     matched-SNP gap, de-trended the same way): hot-site median tail gap barely moves
     (−0.071→−0.055), correlation with bio1 *strengthens* (ρ −0.47→−0.51, p=0.0075→0.0032).
     Climate-slope β's own sign-excess number got the same direct test: bio1 ρ +0.38→+0.42, bio18
     ρ −0.54→−0.59 — also survives. **Neither is explained by the artifact that nulled the median.**
   - Parallelism/PicMin **shrank** on the Kf_w rerun (not re-verified null either way);
     parallelism-by-site's climate-gradient is unchanged but hasn't had this direct test run yet.
   - Note: climate-slope β's original dismissal additionally cited "+0.54/−0.16, weaker/mixed" —
     that was a **mislabeled, unrelated statistic** (general purging-intensity-vs-climate, not the
     SV-specific finding), not just unverified evidence.
   None of these views are redundant with each other — they test a different, still-standing
   question (tail behavior, not central tendency) and none of them have been explained away.
   **Consolidated 2026-07-15** into one notebook (`temporal_s_consolidated.ipynb`) with all four
   sections above, since keeping them as 7 separately-named notebooks was part of why they got
   wrongly assumed redundant in the first place; the old individual notebooks/build scripts were
   removed. → **`SV_TEMPORAL_PURGING_SUMMARY.md`** for the audited detail.

2. **Non-SNP layer adds ~nothing to the polygenic signal ("no kMate gain").**
   Genome-wide and SNP-untagged kinship, per-marker GWAS peaks, per-site and multi-site scans
   all say SNPs and non-SNP markers tell the same story (K_snp/K_nonsnp corr 0.998; joint LRT
   n.s.; same top peaks). SVs/indels behave as passengers at the *variance-partition / peak*
   level. → **`VAREXP_SELECTION_HANDOFF.md`**

3. **Block- and haplotype-unit SV "enrichment", and per-locus climate GEA, are null under
   proper nulls.** Apparent block/hap SV enrichment is a unit/size artifact that collapses
   under rotation/spatial nulls; per-locus climate association is null under site-permutation.
   (Valid at their unit; these motivated moving to the per-variant temporal test in #1.)
   → memory entries `gea-sv-enrichment-clq0.9`, `gea-climate-per-locus-null`.

**Net:** there is no SV-specific *climate-GEA* signal, no *block-enrichment* signal, and — after the
Kf_w / `--unit chrom` rerun — **no genome-wide per-variant *temporal* SV-insertion selection signal
either** (de-trended excess ≈ 0: median +0.0011, 15/31 sites, n.s.). Any residual purging is confined
to a few of the hottest gardens.

## Document index

**Standing decisions / methods**
| Doc | What it is |
|---|---|
| `GLOBAL_MODE_DECISION.md` | Why evolved AF is estimated in **GLOBAL mode** (window/block recombination-detection is circular, underpowered, and its sim floor doesn't transfer). Includes the methods sentence to cite. |
| `PIPELINE_B_POOLED_MODEL.md` | **⚠️ Superseded/dead** — pooled-trajectory model spec on a per-unit window-mode `h`; window mode is not being rerun, so this has no live input (see doc banner). Retained for the record. |
| `WINDOW_UNIT_VALIDATION.md` | **Block / dynld-unit definitions** (the GEA test units) + window-mode cohort run and window-vs-global AF validation. (Window mode as an AF *estimator* is superseded by the global-mode decision; the blocks are retained as test units.) |

**Results**
| Doc | What it is |
|---|---|
| `SV_TEMPORAL_PURGING_SUMMARY.md` | The per-variant temporal SV-insertion selection result (finding #1) — estimator, climate-slope β, parallelism, PicMin, artifact checks, caveats. |
| `VAREXP_SELECTION_HANDOFF.md` | SNP vs non-SNP variance partition + class-split GWAS (finding #2) — kinship redundancy, SNP-untagged GRM null, peak concordance, non-SNP-only candidate genes. |

**Plans / ops**
| Doc | What it is |
|---|---|
| `BAYPASS_TEMPORAL_PLAN.md` | **Plan (not completed** — core model smoke-tested only): Omega-conditioned per-haploblock temporal C2 contrast, kinship-controlled. Includes the method-comparison research appendix. |
| `EXPORT_MANIFEST.md` | Phase-1 GEA WZA results exported to the lab shared Drive (for the advisor's talks). |

## The target (from the phase-1 SNP paper)

The headline phase-1 signature: **an allele whose frequency rises in hot sites and falls in
cold ones** (climate-dependent directional selection). The phase-1 analyses that establish it,
which we mirror for SVs:

| Phase-1 (SNP) | What it computes | SV equivalent here |
|---|---|---|
| **Δp** (`delta_p_*.csv`) | per-SNP AF change from p0, per `site_gen_plot` | per-SV Δp matrix |
| **AF gradient vs temperature** | slope of Δp on site temperature | per-SV gradient |
| **Climate GEA** | Kendall-τ, LFMM (bio1–19, K=16), aggregated by **WZA** | same on SV AF (per-locus null under permutation) |
| **Temporal / selection** | E&R / selection detection on freqs | per-variant plot-replicate logit-slope `s` (finding #1) |

## Data

- kMate per-sample AF: `analysis/grenenet_gea/rerun_kfw_hb/{seedmix,evolved}/`
  (production AF dir, regenerated under `--unit chrom` + full-panel Kf_w; the old-panel
  `grenenet_kmate_arch3` was deleted 2026-07-08). Columns `chrom pos ref_len alt_len alt_freq info n_called se`. SVs = `ref_len!=1 | alt_len!=1`.
- Founding p0 (gen 0): mean alt_freq over the 8 SEEDMIX kMate reps (`lib.build_p0`).
- Sample → site/plot/generation/coverage: Table_S5 (`samples_data_fix57.csv`).
- Analysis unit = `site_gen_plot` pool; timepoints merged by flower-weighted mean AF.
- Climate (bio1–19) per site: `…/grenenet-phase1/climate_gwas/`; founder-origin climate under
  `/global/scratch/users/tbellg/gea_grene-net`.

## Environments
- **Notebooks / plotting**: `basic` env (matplotlib/statsmodels **hang** in the `plotting`
  env — compute → npz/JSON, render in `basic`). `lib.py` needs only numpy/pandas.
- **Compute**: `kmate`. Always `hostname` first — compute only on `n*.savio*`, never `ln00X`
  (hook-enforced). Heavy jobs via **sbatch**, not background.
- **GEA methods**: `lfmm_env` (LFMM), `baypass` (BayPass), `r_env` (MCMCglmm / WZA / BigLD).

## Layout
- `lib.py` — loaders, SV filter, rec_key, p0, block collapse, eff_n_founders.
- `phase1_replication/` — phase-1 kendall / lfmm / binomial + WZA replication on kMate AF.
- `notebooks/` — QC, per-variant temporal, parallelism/PicMin, class-split GWAS figures.
- derived outputs → `analysis/grenenet_gea/` (large, gitignored; rebuildable).
