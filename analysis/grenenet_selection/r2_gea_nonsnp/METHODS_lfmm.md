# Methods — every step and filter producing the LFMM p-values

GrENE-Net climate GEA, generation 9. Written 2026-08-25 for the manuscript. Every
claim below was checked against the code and, where a filter could be inert or
load-bearing, against the data. `file:line` citations are the implementation.

**Scope.** This document ends at the per-record p-value. Block aggregation and hit
calling are downstream and are summarised in §9 only so the p-value's role is clear.
WZA block aggregation, the quasi-binomial model and the Kendall model were retired
2026-08-25 (`../archive_gea/README.md`); this is the sole live track.

---

## Summary of the chain

```
kMate arch3 panel  ->  per-sample AF store  ->  flower-weighted pool matrices
   -> class split (snp / sv / smallindel / nonsnp)   [build_class_matrices.py]
   -> MAF filter                                     [the ONLY record filter that removes anything]
   -> Y = Δp = AF(pool) − p0(founder)                [build_lfmm_input.py]
   -> X = z-scored climate axis                      [22 axes: bio1–19, pc1–3]
   -> LFMM ridge K=16, raw p                         [run_lfmm_nogif.R]
   -> per-record p-value                             <- THIS DOCUMENT ENDS HERE
   -> clq0.9 tiling block assignment                 [reblock_blockdef.py]
   -> Bonferroni per (class × axis)                  [raw_block_significance.py]
```

Final record counts entering LFMM at gen 9, 352 pools across 31 field sites:

| class | definition | records |
|---|---|---|
| `snp` | `ref_len == 1 & alt_len == 1` | 2,016,071 |
| `sv` | `abs(alt_len − ref_len) > 50` | 27,586 |
| `smallindel` | non-SNP, `abs(alt_len − ref_len) <= 50` | 676,203 |
| `nonsnp` | `sv` + `smallindel` pooled | 703,789 |

`nonsnp` is exactly `sv + smallindel` (27,586 + 676,203 = 703,789), so it is a pooled
re-test of the same records, **not** an independent class. Do not treat `nonsnp` and
its two components as three independent tests.

---

## 1. Analysis unit — the pool

The unit is a **flower-weighted `site_gen_plot` pool**, matching the phase-1 SNP paper.
At gen 9 there are **352 pools over 31 sites** (`results/class_matrices/gen9.pools.csv`).
Pools are not independent: ~11 plots per site share a climate value, so the effective
number of distinct climate observations is 31, not 352. This matters for interpreting
the inflation in §7.

Source matrices: `common/results/pool_matrices/pool_gen9_{snp,nonsnp}_af.npy`
(`[pools × records]`, float32, NaN = missing) with row metadata
`pool_gen9_{snp,nonsnp}.meta.csv`.

## 2. Class split — `build_class_matrices.py:46-52`

Two source matrices exist: an all-SNP one and an all-non-SNP one. `snp` and `nonsnp`
take their whole source (`class_mask` returns all-True); `sv` and `smallindel` are
carved from the non-SNP source by indel length:

```python
SV_MIN_BP = 50
dlen = abs(alt_len - ref_len)
sv:         dlen  > 50
smallindel: dlen <= 50
```

The 50 bp SV threshold is the only length parameter. Note it is applied to
`abs(alt_len − ref_len)`, i.e. **net length change**, not alignment span.

## 3. The MAF filter — `build_class_matrices.py:100-113`

This is **the only filter that removes records.**

```python
finite  = np.isfinite(sub)
n_finite = finite.sum(axis=0)
p_bar   = nansum(...) / max(n_finite, 1)      # NaN-aware mean AF across gen-9 pools
maf     = min(p_bar, 1 - p_bar)
keep    = (maf >= 0.05)
```

- `p_bar` is the **contemporary** (gen-9) mean allele frequency across pools, not the
  founding frequency. MAF is therefore defined on the same generation being tested.
- Threshold `--maf-min`, default **0.05**, applied as `>=`.
- Invariant records (`p_bar` 0 or 1) fall out automatically, since MAF = 0 < 0.05.

**A second MAF filter is applied downstream** at analysis time (`MAF > 0.05`, strict).
Verified inert: **zero records in any class sit exactly at MAF = 0.05**, so the
`>=` / `>` mismatch removes nothing and the two filters are equivalent in practice.
Record counts are identical before and after.

### The coverage filter is disabled and inert

`build_class_matrices.py:113` carries a commented-out term:

```python
keep = (maf >= maf_min)  # & (n_finite >= ceil(min_finite_frac * n_pools))
```

kMate projects allele frequencies from a per-chromosome founder-haplotype
reconstruction, so **every record is finite in every pool** — `n_finite == n_pools`
for all records, verified across all generations × classes (~8.5 M records each). The
`--min-finite-frac 0.5` option is retained, not deleted, so the term can be re-enabled
if the pipeline is reused on data that genuinely has missingness (e.g. short-read
pools). `n_finite` is still computed and written to `records.csv` so the claim stays
auditable.

**Consequence for the paper: there is no call-rate or coverage filter in effect.**
State this positively rather than implying one was applied.

## 4. The response variable Y — `build_lfmm_input.py:45-57`

```python
dp = af - p0[None, :]          # Δp, [pools × records]
```

- `af` = gen-9 pool allele frequency.
- `p0` = **founding** allele frequency per record, from
  `common/results/af_store/p0_{snp,nonsnp}.npy`, subset to the kept records by their
  stored source-column index. `sv`, `smallindel` and `nonsnp` all read `p0_nonsnp.npy`.
- Y is therefore the **allele-frequency change from the founding population**, matching
  phase 1. It is *not* raw frequency, and *not* a per-generation difference.

### NaN imputation is present but inert

```python
colmean = np.nanmean(dp, axis=0)               # per-record mean Δp across pools
dp[~isfinite(dp)] = take(colmean, ...)         # impute with that record's column mean
```

Because of §3, nothing is missing: the gen-9 runs report **`imputed 0 cells (0.0%)`
for all four classes** (352 × 2,016,071 for snp, and likewise sv / smallindel /
nonsnp). Describe imputation as implemented-but-unused, or omit it — do not describe
it as applied.

Y is written as raw float64, C-order `[n_pools × n_records]`, and the matrix is deleted
after the model runs (climate-independent and rebuildable).

## 5. The predictor X — the climate axes

22 axes are tested per class: **bio1–bio19** (WorldClim bioclim variables carried on the
pool metadata) and three **climate principal components**.

PCs are computed on the 19 standardized bioclim variables across the 352 pools
(`build_class_matrices.py:75-87`):

```python
Z = (m - m.mean(0)) / m.std(0, ddof=0)
U, S, Vt = svd(Z - Z.mean(0))
score_k  = U[:,k] * S[k]      then z-scored
```

| PC | variance explained | dominant loadings | orientation |
|---|---|---|---|
| pc1 | 46.8% | bio18 −.31, bio14 −.28, bio17 −.28, bio1 +.28 | anchored to **+bio1** |
| pc2 | 23.4% | bio7 −.39, bio2 −.38, bio4 −.37, bio13 +.32 | **+bio7** (continentality) |
| pc3 | 11.7% | bio10 +.38, bio5 +.36, bio2 +.33, bio1 +.31 | **+bio10** (summer heat) |

Cumulative 82.0%. The three PCs are mutually orthogonal by construction
(pairwise correlations ≤ 2e-16, verified).

**Sign convention.** pc1 keeps the historical +bio1 anchor. pc2/pc3 have no natural
anchor, so they are oriented such that the bioclim variable with the largest absolute
loading receives a positive loading — deterministic and reproducible rather than
chosen by hand. **Sign does not affect the p-value** (two-sided test on a single
predictor); it fixes only the direction in which an effect is read.

The PCA is computed **per pool**, so sites contribute in proportion to their plot count
rather than equally. This is inherited from the original pc1 definition and was kept
for consistency; the added pc1 column reproduces the stored one to 2.2e-16.

Finally, the chosen axis is **z-scored** before being handed to LFMM
(`build_lfmm_input.py:59`, `ddof=0`), the LFMM convention.

## 6. The model — `multiaxis/run_lfmm_nogif.R`

```r
mod <- lfmm_ridge(Y = Y, X = X, K = 16)
pv  <- lfmm_test(Y = Y, X = X, lfmm = mod, calibrate = "gif")
write.csv(data.frame(pval = as.numeric(pv$pvalue)), outp)   # RAW
writeLines(as.character(pv$gif), gifout)                    # lambda, recorded not applied
```

- **Latent-factor ridge LFMM, K = 16.** The K decision is documented in
  `notebooks/10_lfmm_k_selection.ipynb` and is unchanged.
- One climate axis per run (X is a single column), so each (class × axis) is a separate
  model fit. 4 classes × 22 axes = 88 fits.
- Y and X are both mean-centred inside `lfmm_test` (`scale(..., scale = FALSE)`).

### The reported p-value is RAW — verified, and the `calibrate` argument is required

`calibrate = "gif"` is passed but **the p-values are not GIF-corrected.** From the
installed package source, `lfmm_test` sets `hp$pvalue` unconditionally *before* the
calibration branch, and the branch only *adds* fields:

```r
hp$pvalue <- hp$pvalue[, 1:d, drop = FALSE]        # RAW, set unconditionally
...
else if (calibrate == "gif") {
    hp$gif               <- compute_gif(hp$score)
    hp$calibrated.score2 <- sweep(hp$score^2, 2, hp$gif, FUN = "/")
    hp$calibrated.pvalue <- compute_pvalue_from_zscore2(...)   # SEPARATE field
}
```

`hp$pvalue` is never overwritten. The argument must be passed to obtain `hp$gif` at
all; with `calibrate = NULL` the branch is skipped and λ could not be recorded.

Empirical confirmation: λ recomputed independently from the written p-values
(`median χ²(1df) / 0.4549`) matches lfmm's own `$gif` closely — snp bio1 2.438 vs
2.445, bio15 2.672 vs 2.681, pc2 0.760 vs 0.762, pc3 3.134 vs 3.146. Had `$pvalue`
been the calibrated one, λ derived from it would be ≈1.0 by construction.

### Why GIF is deliberately not applied

GIF divides the test statistics by λ. That is a valid **deflation** only when λ > 1.
On axes where LFMM K=16 over-corrects (λ < 1), applying it would **inflate** the
statistics and manufacture signal. Since both regimes occur in this data, λ is recorded
per axis and never applied (project decision 2026-07-03, `run_lfmm_nogif.R:3-8`).

## 7. Genomic inflation — must be reported alongside any count

λ is a property of the axis, and is near-identical across the four variant classes on
every axis, which indicates residual population structure not absorbed by the 16 latent
factors rather than anything marker-class-specific.

| | λ (snp) | note |
|---|---|---|
| median over the 20 original axes | **1.72** | 19/20 axes exceed 1 |
| pc3 | **3.13** | most inflated axis in the study |
| bio15 | 2.67 | most inflated bio* axis |
| bio10 / bio1 | 2.52 / 2.44 | |
| bio19 | 0.92 | mildly deflated |
| **pc2** | **0.76** | deflated; see below |

**pc2 is genuinely deflated, not a failure.** λ = 0.740–0.762 consistently across all
four classes; the p-value distribution is directly conservative (2.56% of SNPs below
p<0.05 and 0.36% below p<0.01, against 5% and 1% expected under the null); the smallest
p across 2,016,071 records is 3.45e-05, three orders of magnitude short of the
Bonferroni bar; and the env vector was verified to match `pools.pc2` exactly. The
interpretation is that the continentality axis is largely spanned by the K=16 latent
factors, so the model absorbs the signal being tested. pc2 returns **zero** significant
records in every class.

**Because inflation is not corrected, raw p-values are not calibrated and any count
derived from them is an upper bound.** The ranking of axes by hit count closely tracks
their ranking by λ, which is what one expects if inflation rather than climate drives
the count. Report λ next to every count.

## 8. What is *not* done

Stated explicitly so the Methods section does not over-claim:

- **No GIF / genomic-control correction** (§6).
- **No coverage or call-rate filter** — the term exists but is inert (§3).
- **No NaN imputation in effect** — 0 cells imputed (§4).
- **No WZA or other block aggregation of p-values** — retired (`../archive_gea/`).
- **No multiple-testing correction across the 22 axes or the 4 classes.** Bonferroni is
  applied *within* a (class, axis) scan only.
- **No LD pruning of records** before testing.

## 9. Downstream (for orientation only)

Records are assigned to **clq0.9 BigLD tiling blocks** — LD islands defined at r² ≥ 0.9
on a MAF>0.05 / call-rate≥90% common-variant set, converted to a gap-free partition by
HapFM's rule so 100% of records keep a block (`blocks/blocks_tiling.py`,
`multiaxis/reblock_blockdef.py --how tiling`). 58,376 blocks tile 119,146,348 bp.
`merge_small_blocks` folds sub-2-record blocks into the previous block, **per class**.

A block is called a hit if ≥1 of its records clears **Bonferroni 0.05/n** within that
(class, axis) scan (`multiaxis/raw_block_significance.py`). BH-FDR is computed but not
used for hit calling: on raw p it calls 60–85% of the genome.

Two properties to carry into any block-level claim:

1. Block sizes are severely right-skewed — median 756 bp, mean 2.04 kb, p99 18.1 kb,
   max 1.22 Mb. Bonferroni-hit blocks average 10–100× the median, so "% of genome
   significant" is substantially size-driven and should be reported as an upper bound
   next to the block count.
2. Per-axis rows cannot be summed — the bioclim axes are strongly correlated and hit
   the same blocks repeatedly. Use the union table
   (`results/multiaxis/raw_block_significance_lfmm_union.csv`).

---

## Reproduce

```bash
PY=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python
PR=analysis/grenenet_selection/r2_gea_nonsnp/phase1_replication

# 1. class matrices + climate PCs (gen 9)
$PY $PR/build_class_matrices.py --classes snp sv smallindel nonsnp --gens 9

# 2. LFMM, one array task per (axis, class)
sbatch $PR/multiaxis/ma_lfmm.sbatch          # bio1-19 + pc1, 3-class split
sbatch $PR/multiaxis/ma_nonsnp_lfmm.sbatch   # pooled nonsnp
sbatch $PR/multiaxis/ma_lfmm_pc23.sbatch     # pc2/pc3, all 4 classes

# 3. block assignment
$PY $PR/multiaxis/reblock_blockdef.py --axis <axis> --cls <cls> \
      --blockdef clq0.9 --how tiling --models lfmm

# 4. hit counts + genome fraction
$PY $PR/multiaxis/raw_block_significance.py --model lfmm \
      --axes bio1 ... bio19 pc1 pc2 pc3
```

Key outputs: `results/multiaxis/lfmm/lfmm_{cls}_gen9_{axis}.csv` (raw p),
`..._gif.txt` (λ), `results/multiaxis/wza_in_clq09_tile/` (block-assigned; the name is
historical — it is the **raw per-record** table, not WZA output),
`results/multiaxis/raw_block_significance_lfmm{,_union}.csv`.
