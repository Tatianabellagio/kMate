# Per-garden GWAS on the ecotype selection coefficient — methods of record

Everything needed to write, referee, or reproduce the Result-3 GWAS: what the trait is, every
filter and its threshold, the model, the software, and what is *not* claimed.

**Scope, set 2026-08-25.** This analysis is **30 independent per-garden scans and nothing
else.** There is no cross-site meta-analysis. The all-sites / climate question is answered by
the GEA track (`r2_gea_nonsnp/`), so the earlier Bolormaa et al. 2014 JOINT/GLOBAL/CLIMATE meta
was dropped rather than carried forward. The practical consequence is that **every number this
pipeline reports now comes out of a single community-standard tool (GEMMA)** — there is no
in-house statistic left in the path.

---

## 1. The trait

Per-founder, per-garden **selection coefficient `s`** = the logit-slope of a founder's
genome-wide frequency over generations 0→3.

| item | value |
|---|---|
| built by | `build_selection_trait.py` → `results/varexp/selection_s_matrix.npz` |
| shape | `S[30 gardens × 231 founders]` |
| founder frequency `h` | kMate GLOBAL-mode (`--unit chrom`) chrom-averaged, from `results/ecotype_fitness/sample_global_h.npz` |
| generation 0 reference `p0` | mean over the 8 SEEDMIX replicates (**estimated, not forced uniform 1/231**) |
| generations used | 0, 1, 2, 3 — **gen-1 anchor required**, a plot without it is skipped |
| within a garden | per-plot OLS logit slope, then **mean over plots** (plots averaged, not pooled, so a low-flower plot cannot dominate) |
| timepoints → pool | `site_gen_plot`, collapsed by **flower-weighted** mean `h` |
| logit clip | `EPS = 1e-4`, i.e. `p` clipped to [1e-4, 1−1e-4] |
| founder filter | **none** — all **231/231** retained |
| distribution | mean −0.866, sd 0.649, **skew −0.88, excess kurtosis +1.40** |

**Why a logit-slope and not raw `h` or Δ`h`.** Raw founder frequencies are unusable as a GWAS
trait: skew **+12**, excess kurtosis **+237**, and ~10 founders carry **52%** of the
across-founder variance — a "231-founder" GWAS on `h` is really a ~10-clade GWAS. The
logit-slope is near-symmetric (see table) and spreads variance across founders.

**Why estimated `p0` rather than uniform 1/231.** Founders that kMate's EM cannot separate read
near-zero at *founding and every generation alike*, so the identifiability bias largely cancels
in a slope. The seed-mix mean is exactly 1/231 by closure; the spread around it is real and
reproducible across all 8 replicates.

**No p0 threshold, and no floor.** All 231 founders are analyzable. A `P0_FLOOR = 1e-4` guard
was carried until 2026-08-25 and removed: it floored **0 of 231** founders (lowest estimated p0
is founder 9977 at **1.54e-4**), and since `logit()` clips at the same `EPS` it could not have
changed a slope regardless. The historical `p0 > 1e-3` **drop**-threshold (which kept only 212
founders) was retired earlier; on current post-Kf_w data it would drop only 3.

**No reliability / cross-chromosome-SD weighting** — chromosome-averaging of `h` already
regularizes.

**Gardens:** 30, mean annual temperature (bio1) **6.0–21.6 °C**; plots per garden min 5,
median 12, max 12.

### Trait transform entering the model
Rank-inverse-normal **within each garden**: `y = Φ⁻¹((rank(s) − 0.5) / 231)`. Applied per
garden, so between-garden differences in the spread of `s` cannot drive the association.

---

## 2. Genotypes

Founder panel **arch3**, `panel/arch3/chr{1..5}/var_pa_231_arch3_chr{N}.*`.

- **Coding is haploid 0/1** (presence of the alternate allele in an inbred founder), *not* the
  diploid 0/1/2 dosage. This matters downstream — see §4.
- **Missing genotypes** are imputed to the marker's own alternate-allele frequency `p`
  (mean-imputation), before any modelling.
- **Variant classes**, from the record's own `ref_len`/`alt_len`:
  - `snp` — `ref_len == 1 and alt_len == 1`
  - `sv` — `|alt_len − ref_len| > 50 bp`
  - `indel` — everything else
  - `nonsnp` — indel ∪ sv (the primary 2-way contrast is snp vs nonsnp; sv-only is the
    secondary breakdown)

---

## 3. Filters (all of them)

### Test markers — the markers actually scanned

| filter | threshold | rationale |
|---|---|---|
| minor allele count | **MAC ≥ 5** of 231 → MAF ≥ **2.16%** | matches `build_class_grms.py`, so the GWAS and the variance-partition are on the same marker set |
| call rate | **≥ 90%** | |
| MAC definition | `min(n_alt, N − n_alt)` — **symmetric** | excludes monomorphic (MAC 0) and near-fixed (MAC 1) markers at *both* tails |

Markers passing, per class:

| class | markers scanned |
|---|---|
| snp | 1,752,846 |
| nonsnp (indel + SV) | 525,043 |
| sv (strict, >50 bp) | 12,789 |

**Low-frequency composition — state this in the paper.** The MAC ≥ 5 floor is permissive, and
the classes are not equally exposed:

| class | MAC 5–11 (MAF < 5%) | MAC < 23 (MAF < 10%) |
|---|---|---|
| snp | 36.6% | 58.2% |
| nonsnp | 35.8% | 57.6% |
| **sv** | **55.1%** | 74.8% |

Over half the SV panel sits below MAF 5%, so SV hit lists are the most exposed to
low-frequency instability. Planned check: whether Bonferroni hits are enriched at MAC < 12.

### GRM markers — the kinship correction

| filter | threshold |
|---|---|
| minor allele count | **MAC ≥ 12** |
| call rate | **≥ 90%** |
| classes | **all pooled** (snp + indel + sv together) |

Markers per chromosome: chr1 378,432 · chr2 222,580 · chr3 250,141 · chr4 230,633 ·
chr5 366,518 — **genome total M = 1,448,304**.

**One shared GRM across all three classes, deliberately.** The kinship correction is built from
all classes pooled and is identical for the snp, nonsnp and sv scans, so any difference between
class results is attributable to the **test** markers and not to also varying the correction.
This is justified by the variance-partition result that `corr(K_snp, K_nonsnp) = 0.998`.

---

## 4. Model

Univariate linear mixed model, one garden at a time:

```
y = μ + xβ + g + e        g ~ N(0, σ²_g · K_LOCO)     e ~ N(0, σ²_e · I)
```

- `y` — rank-inverse-normalized `s` for that garden (n = 231 founders)
- `x` — the test marker
- **`K_LOCO`** — leave-one-chromosome-out kinship. Standardized `Z = (g − p)/√(p(1−p))`,
  `K_chr = ZZᵀ`, and for chromosome *c*: `K_LOCO(c) = (Σ_chr K_chr − K_c) / (M − m_c)`.
  LOCO avoids proximal contamination (the tested marker contributing to its own correction).

**Software: GEMMA 0.98.5**, `-lmm 1` (Wald test), with `K_LOCO` supplied via `-k`. GEMMA
re-estimates the variance component **per marker** (exact Wald), which is the stricter, more
accurate form.

### GEMMA invocation and why its own QC is disabled

```
gemma -g <bimbam> -p <pheno> -k <K_LOCO> -lmm 1 -n <garden> \
      -maf 0 -miss 1 -r2 1 -hwe 0
```

GEMMA **does** filter by default (`-maf 0.01 -miss 0.05 -r2 0.9999`). We disable all of it
because QC is done upstream and **GEMMA's MAF filter is wrong for haploid coding**:

- GEMMA assumes diploid dosage, so it computes `maf = mean(geno)/2` — i.e. `p/2` for our 0/1
  coding. Verified: predicted vs actual drops were 2,921 vs 2,927 at `-maf 0.02` and 5,635 vs
  5,634 at `-maf 0.05`.
- That makes it **one-sided**: it can never exceed 0.5, so it cannot catch a near-fixed marker.
  A probe with **230/231 founders carrying ALT (true MAC = 1) passed GEMMA's default QC**, read
  as maf 0.4978.
- A fully monomorphic probe *was* caught — but by `-r2 0.9999` (constant vector collinear with
  the intercept), not by `-maf`, and `-r2 1` disables that too.

Our upstream symmetric `MAC ≥ 5` is strictly stronger than anything GEMMA would apply, so
nothing degenerate reaches it. (At GEMMA's default `-maf 0.01` it would drop 0 of our markers
anyway — MAC 5 lands at maf 0.0108 — but that is a 0.0008 margin, i.e. a coincidence, so we do
not depend on it.)

### Implementation note — marker chunking
GEMMA 0.98.5 aborts with `Enforce failed for N>0 in src/fastblas.cpp … fast_cblas_dgemm` on any
run of **≥ 20,000 markers**, *after* completing the scan but *before* writing output. Verified
boundary: 19,999 succeeds; 20,000 / 20,001 / 21,000 / 30,000 / 45,000 all fail; independent of
`OPENBLAS_NUM_THREADS` and of the QC flags. Markers are therefore chunked below the limit and
the per-chunk association output concatenated. This is a tool bug, not a data property, and has
no effect on results — the repeated per-run cost is only a 231×231 eigendecomposition.

---

## 5. Multiple testing

Corrections are computed **within each garden × class**, because each garden is an independent
scan:

- **Bonferroni** 0.05 / M (M = markers tested in that garden × class)
- **Benjamini–Hochberg** q < 0.05

**What is not corrected, and must be said plainly:** the 30 gardens × 3 classes = 90 scans are
**not** corrected against each other. A marker significant in one garden is a per-garden result.
Comparing gardens (e.g. "hit in 5 gardens") is **descriptive, not a tested contrast** — the
cross-garden test was exactly the meta that this analysis deliberately dropped.

**Genomic-control λ is reported per garden × class as a diagnostic and is not used to rescale
the statistics** — the LOCO kinship term is the structure correction.

---

## 6. Validation of the estimator

Before switching, the in-house EMMAX/P3D scan this replaces was benchmarked against GEMMA on
identical GRM, filters and founder order (`gemma_validation.py`; chr1 × snp, 20,000 markers ×
3 gardens):

- **pearson(z) = 0.9975**, spearman 0.9974
- the in-house arm was very slightly **conservative** — mean(|z_ours| − |z_gemma|) = −0.013 at
  |z| > 2, more significant 47.7% of the time
- that direction is the expected P3D signature: the in-house scan fixed δ under the null per
  garden × chromosome, GEMMA re-estimates per marker

So the switch to GEMMA is not expected to change conclusions, only to make them citable — and
if anything to yield slightly *more* hits.

---

## 7. Provenance

| stage | script | output |
|---|---|---|
| trait | `build_selection_trait.py` | `results/varexp/selection_s_matrix.npz` |
| pooled per-chrom GRMs | `gemma_validation.py::build_pooled_grms` | `results/gemma_validation/pooled_grms.npz` |
| per-garden scan | `gemma_persite_gwas.py scan` + `run_gemma_gwas.sbatch` (15-task array) | `results/gemma_gwas/parts/{chrom}_{class}.npz` |
| assembly + significance | `gemma_persite_gwas.py assemble` + `run_gemma_assemble.sbatch` | `results/gemma_gwas/persite_gwas_{class}.npz`, `persite_gwas_summary.json` |
| estimator validation | `gemma_validation.py` | `results/gemma_validation/{concordance.csv,summary.json}` |

Environments: `kmate` (python) + `gwas_tools` (the `gemma` binary; kept separate so the solver
cannot bump the OpenBLAS that `kmate`'s numpy links against). Compute nodes only.

### Methods paragraph, short form

> Per-garden association was tested with a univariate linear mixed model in GEMMA v0.98.5
> (`-lmm 1`, Wald), fitting the rank-inverse-normalized per-founder selection coefficient in
> each of 30 gardens across 231 inbred founders. Population structure was controlled with a
> leave-one-chromosome-out genomic relationship matrix built from 1,448,304 markers of all
> variant classes pooled (MAC ≥ 12, call rate ≥ 90%), identical across the SNP, non-SNP and SV
> scans so that class differences reflect the tested markers rather than the correction. Markers
> were tested at MAC ≥ 5 (MAF ≥ 2.16%) and call rate ≥ 90%, giving 1,752,846 SNP, 525,043
> non-SNP and 12,789 SV markers; missing genotypes were imputed to the marker mean. Significance
> was assessed within each garden by Bonferroni and Benjamini–Hochberg; gardens were not
> corrected against one another and cross-garden recurrence is reported descriptively.

### Citations
- **GEMMA** — Zhou & Stephens 2012, *Nat Genet* 44:821–824.
- **LMM / LOCO framing** — Kang et al. 2010 (EMMAX), *Nat Genet* 42:348–354.
- **kMate allele frequencies** — this work; see `ALGORITHM.md`.

---

## 8. Known limitations to disclose

1. **The trait is derived, not measured.** `s` is a logit-slope of kMate-estimated founder
   frequencies, so founder-frequency estimation error propagates into the phenotype. Founders
   that are poorly identifiable in the EM (e.g. 9977) carry near-zero, near-flat slopes.
2. **Low-frequency markers dominate the SV class** (55.1% below MAF 5%) — see §3.
3. **No cross-garden inference.** By design; see §5.
4. **231 founders is a small n for an LMM.** Per-garden power is limited and effect estimates at
   low MAC are unstable.
5. **The 30 gardens are not independent replicates** — they share the same founder panel and the
   same `p0`, so their scans are correlated. This is precisely why cross-garden claims need the
   meta that was dropped, and must not be smuggled in via hit-counting.
