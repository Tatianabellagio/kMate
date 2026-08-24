# BayPass TEMPORAL (founding→evolved) on haploblocks — execution plan

> **⚠️ STALE INPUTS (2026-07-08).** This plan's haploblock-frequency inputs come from the
> hapfreq/Pipeline-B chain, now RETIRED (archived under `archive/pipelineB_hapfreq_retired/`;
> it read the deleted `results/grenenet_kmate_window` store and was built on the pre-Kf_w window
> `h`). `hapfreq/hapfreq_matrix.npy` + `hapfreq_registry.csv` are therefore stale. Before executing,
> repoint to a haploblock frequency table regenerated under the corrected `--unit chrom` + Kf_w
> cohort (or retire this plan). `baypass_build_inputs.py` carries the same warning.

> **Status (2026-07-06): PLAN, not a completed analysis.** Only the CORE (Omega) model was
> smoke-tested at site 4 (`analysis/grenenet_gea/archive/window_hapfreq_retired/hapfreq/baypass_temporal_site4/core_baypass.log`);
> **no C2 contrast results exist yet.** The method-research doc
> `KINSHIP_TEMPORAL_METHODS_RESEARCH.md` is merged into this file (Appendix below). Its full raw
> deep-research claims dump (~400 lines) was dropped in the 2026-07-06 consolidation but is preserved
> in git (commit `5cefcfa`) if ever needed.

**Goal:** kinship/structure-controlled PER-HAPLOBLOCK temporal selection test = BayPass C2
contrast (founding vs evolved), Omega-conditioned. Temporal FIRST; defer spatial (climate
`-efile`) escalation. Rationale + method comparison in the Appendix below.

## 0. Environment — RESOLVED ✅
Binary `~/bin/g_baypass` = **BayPass v2.41** (supports `-contrastfile` C2 + `-poolsizefile` pool model).
It needs `libgfortran.so.3`, which lives in the existing **`baypass` conda env**. Working recipe
(from the user's 2024 jobs, confirmed runs 2026-06-25):
```bash
source ~/miniforge3/etc/profile.d/conda.sh
export LD_LIBRARY_PATH=/global/home/users/tbellg/miniforge3/envs/baypass/lib:$LD_LIBRARY_PATH
conda activate /global/home/users/tbellg/miniforge3/envs/baypass
~/bin/g_baypass -gfile ... [runs]
```

## 1. Why this is fast now (the old run was slow)
Old run = SPATIAL STD model (`-efile env_3rdgen` = climate) on **100k+ hapFIRE SNPs**, split into
~1000 loci-partitions × 3 chains (`individual_gfiles/partition_*.txt`). Slow purely from SNP count.
Haploblocks → ~40k one-vs-rest markers (k-1 per block, drop reference) → **no partitioning needed**,
one whole-genome run per chain.

## 2. Design (temporal contrast)
- **Populations** = founding (8 seedmix reps, gen 0) + evolved gen-3 pools.
  - START SMALL: site 4 → 8 seedmix + 11 gen-3 plots = **19 pops** (fast smoke test).
  - SCALE: all 355 evolved gen-3 pools + founding (matches their omega 355×355).
- **Contrast (C2):** founding = **-1**, evolved = **+1** (one line, one value per pop, 0 = ignore).
  Tests per-haploblock founding→evolved differentiation CONDITIONED on Omega (the genome-wide
  founder relatedness / whole-genome winnowing) → block-specific temporal selection.

## 3. Input files (build from hapfreq_matrix + window h)
BayPass gfile = one ROW per locus; per population two integer counts `(n_allele1 n_allele2)`,
space-separated (their old gfile was comma-sep — confirm delimiter; manual says space).
- **gfile:** rows = haploblock markers (k-1 per dynld block, drop max-panel-freq reference, MAC filter
  — reuse the `founder_gwas_site.py` keep-logic). Per pop: `n1 = round(freq·Neff)`, `n2 = round((1-freq)·Neff)`.
  `freq` = per-population haplotype frequency (hapfreq_matrix value for evolved pools; seedmix reps for founding).
  `Neff` = effective depth per pop (use coverage, or n_called, or pool size).
- **poolsizefile:** one line, haploid pool size per pop (their `pool_sizes_*` format). Engages the pool model.
- **contrastfile:** one line, +1/-1/0 per pop (founding -1, evolved +1).
- **omegafile:** from a CORE-model run first (below).

## 4. Workflow (adapt their working command style)
Their working flags: `-gfile -efile -omegafile -outprefix -seed -pilotlength 1000 -nval 50000 -nthreads 8 -npilot 20`.
1. **Estimate Omega (core model)** — no `-efile`/`-contrastfile`:
   `g_baypass -gfile gfile -poolsizefile psize -outprefix core -nthreads 8 -npilot 20 -pilotlength 1000 -nval 50000 -seed N`
   → `core_mat_omega.out` (the coancestry matrix). Run 2-3 chains, check Omega reproducibility (FMD/correlation).
2. **Contrast / C2** — reuse Omega:
   `g_baypass -gfile gfile -poolsizefile psize -omegafile core_mat_omega.out -contrastfile contrast -outprefix temporal -nthreads 8 -npilot 20 -pilotlength 1000 -nval 50000 -seed N`
   → outputs:
   - `temporal_summary_contrast.out` = per-locus **C2** statistic + calibrated p (`log10(1/pval)`); C2 ~ χ²₁ under null.
   - `temporal_summary_pi_xtx.out` = **XtX\*** (overall per-locus differentiation, calibrated).

## 5. Read results
Per haploblock: C2 (= founding-vs-evolved, structure-controlled), p = 10^(-log10(1/pval)), BH-FDR.
λ_GC of C2 should be ~1 (Omega calibrates). Map hits → genes via `genes_from_regions.py`.
Sanity vs our findings: blocks beyond the genome-wide winnowing background = the block-specific
temporal selection that the founder-GWAS (kinship-on-analysis-2) could not localize at site 4.

## TO VERIFY against saved manual (tool-results webfetch PDF) or `g_baypass` help once lib fixed
- contrastfile exact format (single line, +1/-1/0) and whether multiple contrasts allowed.
- gfile delimiter (space vs their comma) + that `-poolsizefile` auto-engages the Poisson-binomial pool model.
- exact column order of `summary_contrast.out` (C2_std, M_C2, log10(1/pval)).
- MCMC: nval 50000 was their choice; for a smoke test drop to nval 5000 / npilot 10.

---

# Appendix — kinship/structure control for per-haploblock TEMPORAL selection (methods research)

**Date:** 2026-06-25 · **Source:** deep-research workflow `wf_a5bc871b-9f8` (5 search angles, ~15
sources). *(The workflow's "all claims refuted" notification was a FALSE ALARM — the verifier agents
hit a session token limit and defaulted every vote to refuted; verification never ran. The claims are
standard, well-established methods.)*

## Question
Add a kinship/relatedness (population-structure) control to a PER-HAPLOBLOCK analysis of TEMPORAL
allele-frequency change (replicated E&R / pool-seq, known founder panel), keeping temporal+replicate
variance and ALLOWING block-specific residuals — i.e. NOT collapsing to a founder-level GWAS. Setup:
GrENE-net, ~200 inbred founders, near-uniform founding, ~11 plots/site, gens 0–3; ~0 within-experiment
recombination so whole ecotype genomes hitchhike and per-block selection is smeared by genome-wide
founder winnowing.

## Verdict (method comparison)

| Method | "kinship" object | temporal+replicate | per-block test | fit |
|---|---|---|---|---|
| **BayPass Omega / C2 / XtX\*** (Gautier 2015) | founder-coancestry **Omega** (allele-freq kinship analog) | via populations (encode site×gen) | **YES — calibrated C2 ~ χ²₁** | **WINNER** |
| Buffalo & Coop temporal covariance (`cvtkpy`) | time×time covariance Q (empirical) | native | NO — genome-wide variance decomposition only | drift-null + genome-wide sanity check |
| ACER (Spitzer/Pelizzola 2020, R) | scalar drift variance (Ne, supplied) | native (repl) | YES per-locus | controls DRIFT not STRUCTURE; baseline |
| LFMM2 (Caye & François, R `lfmm`) | latent factors | NO (static, not Δp) | YES | wrong data model (not temporal/freq-change) |
| haplovalidate (Otte/Schlötterer) | — | yes | NO (post-hoc) | confirms haplotype-block is the right unit |
| GWAS LMM (EMMAX) | founder GRM K | no | yes | = the founder GWAS we already ran |

## Recommended formulation (Omega-conditioned per-block temporal contrast)
1. **Omega** = among-population coancestry covariance from genome-wide block freqs (the structure to
   remove = whole-genome ecotype winnowing).
2. **Populations = site×generation (or plots×gen)** + seedmix(gen0) → puts time in the design.
3. Per block, standardize the population allele-freq vector and form the **C2 contrast** (founding vs
   evolved) conditioned on Omega → calibrated **χ²₁** per-block p/q. Mechanism = whiten freqs through
   Cholesky of Omega (GLS transform vs shared demography); MVN drift null `freq ~ N(p·1, p(1−p)·Omega)`.
4. **Drift scale** from among-replicate variance (ACER-style) feeds the null.
5. **Genome-wide check** with `cvtkpy`: fraction of Δp variance from linked selection (Buffalo-Coop got
   17–37% in E&R) + across-replicate convergence correlation.

This is the FORMAL, calibrated version of the ad-hoc `local_h − global_h` residual idea: Omega-whitening
= subtracting the genome-wide founder background, but with a proper per-block null, so a block beyond
its coancestry expectation still scores (allows block-specific signal).

## Software & key references
- **BayPass** (core) — Gautier 2015 Genetics; C2 contrast: Olazcuaga et al. 2020 MBE (worked template).
  Manual: gensoft.pasteur.fr/docs/baypass/. Origin: Günther & Coop 2013 Genetics (Omega/XtX).
- **ACER** (R) — Spitzer, Pelizzola & Futschik 2020 AOAS; github.com/MartaPelizzola/ACER.
- **cvtkpy** (Python) — Buffalo & Coop 2019 Genetics + 2020 PNAS; github.com/vsbuffalo/cvtkpy.
- **LFMM2** (R `lfmm`) — Caye & François 2019 MBE. **haplovalidate / poolSeq** — Otte & Schlötterer.
