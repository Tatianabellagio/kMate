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

1. **Per-variant temporal selection — signal did NOT survive the Kf_w / `--unit chrom` rerun.**
   The earlier headline (a real, small, SV-specific climate-graded purging excess, "three methods
   agree") **does not survive** regeneration under `--unit chrom` + full-panel Kf_w. On the rigorous
   frequency-de-trended metric the **SV excess-vs-baseline is ≈0: median +0.0011, negative at only
   15/31 sites, sign test n.s.** — no genome-wide SV purifying excess; any residual is confined to a
   few of the hottest gardens. The raw/single-tail methods (climate-slope β, parallelism, PicMin,
   histograms) still look purged only because they keep the frequency / founder-projection confound;
   the de-trended `s_distribution_by_site` is authoritative. → **`SV_TEMPORAL_PURGING_SUMMARY.md`**
   (see its superseded banner).

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
| `PIPELINE_B_POOLED_MODEL.md` | **Frozen spec** for the pooled-trajectory model: pool plots → site freq → weighted within-site logit-slope `s_g` (variance-components SE) → random-effects climate meta-regression + site-permutation null; block WZA. |
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

- kMate per-sample AF: `results/grenenet_gea/rerun_kfw_hb/{seedmix,evolved}/`
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
- derived outputs → `results/grenenet_gea/` (large, gitignored; rebuildable).
