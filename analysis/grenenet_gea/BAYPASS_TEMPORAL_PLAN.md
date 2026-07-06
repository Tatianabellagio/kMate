# BayPass TEMPORAL (founding→evolved) on haploblocks — execution plan

**Goal:** kinship/structure-controlled PER-HAPLOBLOCK temporal selection test = BayPass C2
contrast (founding vs evolved), Omega-conditioned. Temporal FIRST; defer spatial (climate
`-efile`) escalation. See [[gea-kinship-temporal-methods-research]] and the curated report
`KINSHIP_TEMPORAL_METHODS_RESEARCH.md`.

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
