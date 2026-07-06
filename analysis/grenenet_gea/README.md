# GrENE-Net SV-GEA analysis (kMate)

Re-running the GrENE-Net phase-1 genotype–environment-association (GEA) pipeline
on **structural-variant** allele-frequency time-series produced by kMate, to find
SVs under climate-dependent selection — and ask whether SVs carry adaptive signal
that SNP studies miss.

## The target (from the phase-1 SNP paper)

The headline signature: **an allele whose frequency rises in hot sites and falls
in cold ones** (climate-dependent directional selection). The phase-1 analyses
that establish it, which we mirror for SVs:

| Phase-1 (SNP) | What it computes | SV equivalent here |
|---|---|---|
| **Δp** (`delta_p_*.csv`) | per-SNP AF change from p0, per `site_gen_plot` | per-SV Δp matrix |
| **AF gradient vs temperature** (Fig B) | slope of Δp on site temperature → "AF gradient (°C)" | per-SV gradient |
| **Cross-site Δp scatter** (Fig B) | hot vs cold Δp — anti-correlated for adaptive loci | site-4 vs site-54 Δp |
| **Climate GEA** | Kendall-τ, LFMM (bio1–19, K=16), aggregated by **WZA** | same on SV AF |
| **Selection model** (`stabilizing_selection_model.R`) | MCMCglmm `log(p1/p0) ~ climate`, random site+ecotype | same on SVs |
| **adaptest** | E&R selection detection on hapFIRE freqs | same on kMate SV freqs |

**The new question:** do SV-GEA hits overlap SNP-GEA hits, or do large SVs reveal
adaptive signal SNP-tagging can't see?

## Data

- kMate per-sample AF: `results/grenenet_kmate_arch3/<MLFH...>.tsv`
  (`chrom pos ref_len alt_len alt_freq info n_called se`). SVs = `ref_len!=1 | alt_len!=1`.
- Founding p0 (gen 0): mean alt_freq over the 8 SEEDMIX kMate reps (`lib.build_p0`).
- Sample → site/plot/generation/coverage: Table_S5.
- Climate (bio1–19) per site: `…/grenenet-phase1/climate_gwas/` (for the full cohort).

## Scope: pilot vs full cohort

The full climate-GEA needs **all ~30 sites** (a 2-site pilot can't fit a climate
gradient). The **pilot (site 4 hot / Spain, site 54 cold / Cologne)** is for QC +
building & validating the framework, and a 2-climate preview:
- per-SV Δp in hot vs cold → SVs on the anti-diagonal (up-in-hot / down-in-cold)
  are candidate adaptive SVs;
- their generation trajectories (gen 1→2→3) in each climate.

## Environments
- **Notebooks / plotting**: the `plotting` conda env (the `kmate` pipeline env has
  no matplotlib/jupyter). Run: `conda activate plotting`. `lib.py` only needs
  numpy/pandas, so it imports in either.
- **GEA methods** (full cohort, later): `lfmm_env` (LFMM), `baypass` (BayPass),
  `r_env` (the MCMCglmm selection model / WZA in R) — already on this system.

## Key decisions
- **`GLOBAL_MODE_DECISION.md`** — why evolved AF is estimated in GLOBAL mode (not
  window/block): the recombination-detection test is circular on projected AF,
  underpowered (~3 gens, 97% selfing), and its sim noise floor doesn't transfer
  (panel-incompleteness confound). Includes the methods sentence to cite.

## Layout
- `lib.py` — loaders, SV filter, rec_key, p0, eff_n_founders.
- `notebooks/01_pilot_qc.ipynb` — per-sample QC + data-readiness (run first).
- `notebooks/02_pilot_sv_dp_hotcold.ipynb` — SV Δp, trajectories, hot-vs-cold (next).
- derived outputs → `results/grenenet_gea/`.

## Pilot design (from Table_S5)
Longitudinal, with attrition (the paper's extinction signal):
site 4 (hot): 24→22→11 plots over gen 1→2→3; site 54 (cold): 59→35→24.
Coverage ~6–7.5× (some to 2.8× → coverage filter matters). 2.74M SVs/sample
(small indels to 294 kb).
