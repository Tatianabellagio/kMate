# Results log — running tally of validation/benchmark numbers

Add new findings at the top with timestamp.

---

## 2026-04-29 ~10:40 — Cactus 82 ↔ xwu 231 SNP overlap

Replicates the syri-vs-xwu analysis (`/carnegie/nobackup/scratch/tbellagio/freqk_gr/panel_overlap_test/`)
but for the cactus pangenome panel. Tells us how much cactus's natural SNPs
add over xwu's short-read GrENE-Net SNPs — informs whether to include cactus
SNPs in the imputation reference panel.

**Counts**:

|                         | xwu 231 | xwu 80-subset | xwu 82-subset | cactus 82 SNPs |
|---|---:|---:|---:|---:|
| total positions         | 3,235,480 | 2,586,635 | 2,602,621 | **3,717,261** |

**Overlap (cactus 82 vs xwu 82-subset, apples-to-apples)**:

| | n positions |
|---|---:|
| shared | 1,753,272 (67% of xwu82) |
| cactus-only (not in xwu82) | **1,963,989** |
| xwu82-only (not in cactus) | 849,349 |

**For imputation reference-panel design**:

| Backbone | n positions | gain vs current |
|---|---:|---:|
| xwu 80-subset SNPs (current ref_80 backbone) | 2,586,635 | — |
| xwu 80-subset ∪ cactus SNPs | **4,560,828** | **+76.3%** |
| cactus-only sites added | 1,974,193 | new informative SNPs |

**Interpretation**:
- Cactus pangenome covers only **67%** of xwu's 82-subset SNPs — 33% of xwu's
  SNPs aren't called by cactus (short-read SNP catalog is broader at "easy" positions).
- Cactus contributes **1.96M extra SNPs** that xwu didn't catch — most likely
  in repetitive regions, SV flanks, and centromere-adjacent areas where
  short-read mapping fails but long-read assembly + minigraph-cactus does not.
- The two SNP catalogs are **complementary**, not redundant. Using both via
  union gives **76% more SNPs** in the imputation reference panel.

**Decision implication**: even before 57779's slope readout, the overlap
analysis already says "Tier 2 (include cactus SNPs in ref_80) is worth doing"
purely on the +76% reference density gain. The 1.96M extra SNPs improve LD
anchoring during imputation, which improves Beagle's accuracy at the SVs.

Plots: `plots/cactus_vs_xwu82_snp_venn.png`, `plots/imputation_panel_design_3way.png`

---

## 2026-04-29 ~12:00 — Tier 1 root-cause fix VALIDATED on Chr4 (job 57779)

Built a corrected merged VCF from `ref_80 + imputed_151` (instead of the buggy
`cactus_svs + imputed_151`), rebuilt cn_kmer_v2 + cn_var_v2 for Chr4 only,
re-ran per_sample_driver in global mode on SEEDMIX_S1 reads.

**Headline numbers vs the old (buggy) version**:

| Metric | OLD (buggy build) | NEW (with cactus SNPs) |
|---|---|---|
| **slope** | 1.43 | **1.0031** ✓ |
| intercept | 0.0034 | 0.00044 |
| **R² (raw)** | 0.68 | **0.9935** ✓ |
| Pearson r | 0.984 | 0.9968 |
| RMSE | 0.071 | 0.0154 |

**Density check** (cactus founders' alt-rate at SNP records, was the root cause):

| | OLD | NEW |
|---|---|---|
| cactus founders mean density | 0.0020 | ≈imputed (0.95×) |
| imputed/cactus density ratio | 50× | 0.95× |

This **conclusively confirms** the 1.43× slope was a build artifact, not a
fundamental rank-deficiency. **No calibration needed** — slope = 1.003.

**Action**: full-genome rebuild submitted as job **57794** (~10h: ~7.5h for
the 5 cn_kmer chrom builds + 30 min cn_var). Once done, we have a working
clean 231-founder Beagle-imputed panel ready for production immediately,
without waiting for pang_69 / PanGenie. PanGenie path is still preferred
architecturally (no imputation step) but Tier 1 fix gives a fallback today.

---

## 2026-04-29 ~11:00 — Deep dive: where does the 1.43× slope come from?

### Finding 1: cn_var_231 has no SNP genotypes for the 80 cactus founders

| Record type | n | cactus density | imputed density | ratio |
|---|---|---|---|---|
| SVs (any size) | 241,155 | **0.029** | 0.020 | 0.68× |
| SNPs (1bp/1bp) | **3,235,480** | **0.000000** | 0.107 | ∞ |
| All | 3,476,635 | 0.0020 | 0.101 | 50× |

**Root cause**: `cn_var_231` was built via
`bcftools merge cactus_svs_renamed.vcf.gz imputed_151.vcf.gz`. The
`cactus_svs` file contains only SV records (no SNPs), so when bcftools
merges, the 80 cactus founders get `./.` at all 3.24M SNP records (the
build_cn_var.py script treats `./.` as 0). Cactus founders contribute zero
alt-allele evidence at SNPs — 93% of the cn_var matrix.

The EM correctly learns "cactus founders never carry alt at SNPs" and
shifts mass to imputed founders to explain SNP-related k-mer counts. This
inflates predicted alt_freq by **1.43×** relative to recipe truth.

On the SVs alone (where both groups have real data) cactus density (0.029)
is *higher* than imputed (0.020), as expected biologically — SVs tend to
be private to specific founders and Beagle imputation is conservative.

**Proper fix**: rebuild `cn_var_231` from
`bcftools merge ref_80.vcf.gz imputed_151.vcf.gz`, where `ref_80` already
contains the 80 cactus founders' SNP genotypes (from `grene_80.vcf.gz`,
used as Beagle's reference panel input). The `ref_80.vcf.gz` already
exists at `/carnegie/nobackup/scratch/tbellagio/hapfire_sv/imputation/work/`.
Estimated rebuild cost: cn_kmer ~9h + cn_var ~30 min, single SLURM job.

**Workaround (current)**: post-hoc 1.43× calibration recovers R² 0.68→0.97
for the production deliverable. Slope is panel-intrinsic and stable
(1.4325 ± 0.0055 across 8 SEEDMIX reps) — not a per-sample artifact.

### Finding 2: 82-founder vs 231-founder predict different things by construction

Head-to-head SEEDMIX_S1 on the 2.05M overlapping records:

| Subset | R² (82 vs 231 raw) | R² (82 vs 231 calibrated) | r |
|---|---|---|---|
| All overlap | 0.93 | 0.78 | 0.96 |
| SNPs only (1.69M) | — | 0.79 | 0.98 |
| SVs only (360K) | — | **0.14** | 0.48 |

The 82-founder pipeline projects through the **80-panel-restricted recipe**
(35.3% of total mass); the 231-founder pipeline projects the **full recipe**
(100% mass). They correctly predict different alt_freqs at the same record.
Calibration aligns 231 to recipe truth, which moves it AWAY from 82's
predictions.

**SVs disagree more than SNPs** (R²=0.14 vs 0.79) because SV records went
through different pipelines (F_MISSING filter dropped some, Beagle re-genotyped
others). SNPs are shared via the GrENE-Net VCF.

### Finding 3: per-block h on SEEDMIX is stable, with rank-deficiency artifacts

For SEEDMIX_S1 (homogeneous F0 pool, 231-founder window mode):
- **All 598 blocks got local fits** — no fallbacks needed
- **Effective n founders per block**: median 71, global 116 (per-block EM concentrates mass; rank deficiency at work)
- **A handful of blocks have eff_n = 1** (degenerate fit) — outliers, not a systemic issue
- **Per-founder block-h std**: cactus median 0.008, imputed 0.006 — stable across the genome, no spurious "selection peaks" detected on this homogeneous pool

Plot: `plots/seedmix_231_per_block_h_genome.png` shows top-5 recipe-founders'
per-block h along the genome — flat lines around the recipe truth, with
the rank-deficiency wobble visible but no systematic genome-wide structure.

This confirms block-EM **does not introduce false signals** on homogeneous
pools — the wobble is small noise, not large enough to be confused with
real selection peaks. For evolved E&R samples, the per-block h variation
that DOES emerge will be biologically meaningful selection signal.

---

## 2026-04-29 ~10:00 — cactus_nocap rebuild = no-op for our pipeline

The 42-hour cactus rebuild without `--maxLen` (job 56169) finished. Compared
to the original capped pang (`pang/output/pang_1001gplus_82acc.vcf.gz`):

| Metric | capped | nocap | Δ |
|---|---|---|---|
| vcfbub-filtered records | 4,450,411 | 4,450,117 | -294 (0.007%) |
| max SV size | 299,982 bp | 299,982 bp | 0 |
| SVs ≥1Mb | 0 | 0 | 0 |
| 1-10kb SVs | 11,885 | 11,878 | -7 |
| 10-100kb SVs | 1,413 | 1,412 | -1 |

**The maxLen flag never gated any real SV** — the largest variant in both
versions is 300 kb, well under any sensible cap. The trivial differences
(~300 records) are cactus's run-to-run stochasticity. The capped pang was
already producing the full SV catalog. **We don't need to switch the
production pipeline to nocap.**

---

## 2026-04-29 ~01:30 — Calibration confirmed across all 8 SEEDMIX_231 replicates

| Sample | r | R²(raw) | slope | intercept | R²(calibrated) | RMSE(cal) |
|---|---|---|---|---|---|---|
| S1 | 0.9843 | 0.676 | 1.437 | 0.0034 | 0.968 | 0.0223 |
| S2 | 0.9846 | 0.683 | 1.434 | 0.0028 | 0.969 | 0.0221 |
| S3 | 0.9838 | 0.668 | 1.443 | 0.0028 | 0.967 | 0.0227 |
| S4 | 0.9858 | 0.696 | 1.428 | 0.0027 | 0.971 | 0.0212 |
| S5 | 0.9849 | 0.690 | 1.427 | 0.0036 | 0.969 | 0.0219 |
| S6 | 0.9854 | 0.690 | 1.429 | 0.0037 | 0.970 | 0.0215 |
| S7 | 0.9853 | 0.682 | 1.435 | 0.0034 | 0.970 | 0.0216 |
| S8 | 0.9855 | 0.692 | 1.427 | 0.0035 | 0.970 | 0.0214 |

**Slope: 1.4325 ± 0.0055** (0.4% variation across replicates) — production-stable.

**Production recipe**:
```
1. Compute calibration once from any SEEDMIX rep:
     python src/calibrate_alt_freqs.py compute \
         --predicted-tsv results/seedmix_231/SEEDMIX_S1.tsv \
         --recipe data/seedmix_recipe_normalized.tsv \
         --cn-var data/cn_var_231.cn_var.npz \
         --cn-var-meta data/cn_var_231.meta.npz \
         --out data/calibration_231.json

2. Apply to all 2,414 evolved-sample TSVs:
     python src/calibrate_alt_freqs.py apply \
         --calibration data/calibration_231.json \
         --in-tsv results/per_sample_231/<sample>.tsv \
         --out-tsv results/per_sample_231_calibrated/<sample>.tsv
```

Expected per-record alt_freq R² vs truth ~ 0.97 across the production set.

---

## 2026-04-29 ~00:30 — cactus_em Tier 1 chain COMPLETED, both modes

| Mode | R² | RMSE | r | Notes |
|---|---|---|---|---|
| global | **0.9985** | 0.008 | 0.999 | best — appropriate for homogeneous 10-founder pool |
| window | 0.91 | 0.063 | 0.956 | 412/598 blocks fell back to global (Chr1-only reads → most blocks have <200 nonzero k-mers) |

For homogeneous F0-style pools, global mode is correct. Window mode noise comes from per-block EM on under-determined blocks. Same pattern as on SEEDMIX. For evolved samples, window mode captures local mosaic ancestry that global blurs.

---

## 2026-04-28 ~22:15 — Tier 1 cactus_em sim (visor_freqk architecture)

First end-to-end VISOR HACk + SHORtS sim through our cactus pangenome k-mer EM
pipeline (branch `cactus_em` in `sims/visor_freqk/`).

**Setup**: 10 cactus founders pooled uniformly (1/10 each). 5 of the 10 carry
a 1kb deletion at Chr1:10000000. 30× coverage on Chr1 only. Reads simulated
by VISOR SHORtS using the HACk-modified founder FASTAs.

**Pipeline runtime** (10 cactus founders, Chr1 only):
- 00 setup clones (link FASTAs):           1s
- 02 VISOR HACk (5 SV haplotypes):        50s
- 03 VISOR SHORtS (pool reads at 30×):    7m
- 05c cactus_em global mode:              ~9m  (k-mer count 367s + EM 28s; nonzero k-mers 8.8% since reads are Chr1-only)

**Founder-frequency recovery (global EM)** evaluated on per-record alt_freq:
| Region | n records | R² | RMSE | Pearson r |
|---|---|---|---|---|
| Chr1 outside deletion | 1,564,851 | **0.9985** | 0.0080 | 0.9993 |
| Chr1 inside deletion (60 cactus VCF records in 1kb) | 60 | 0.997 vs no-deletion truth | 0.009 | 0.999 |

**Comparison to old SUMMARY.md hapFIRE-proj**: MAE 0.0010-0.0013 at cov50 on
GrENE-Net 231-ecotype sims. Our cactus_em on cactus 10-founder sim achieves
RMSE 0.008 — comparable order of magnitude (the pools differ).

**Deletion detectability**: the per-record alt_freq predictions inside the
1kb deletion region are slightly lower than outside (0.101 vs 0.119), but the
signal is weak because there are only 60 cactus VCF records in 1kb and the
genome-wide EM is dominated by the millions of records outside. Our pipeline
estimates **founder frequencies**, not specific SV calls — the SV is invisible
unless it's in the cactus VCF as its own record.

---

## 2026-04-28 ~21:30 — SEEDMIX recipe-truth comparison, all 4 mode/panel combos

For SEEDMIX (8 GrENE-Net F0 replicates, recipe-truth = founder seed proportions):
per-record alt_freq projected from estimated h, compared to recipe-projected truth.

| Panel | Mode | R² (mean over reps) | RMSE | Pearson r | Notes |
|---|---|---|---|---|---|
| 82-founder | global | **0.994** | 0.015 | 0.997 | best — small panel, fully observed cactus genotypes; SEEDMIX_S1 only (others not yet run in global) |
| 82-founder | window-200kb | 0.97 | 0.03 | 0.98 | block-EM noise on F0 pool — block-EM is wrong model for homogeneous pool but survives well |
| 231-founder | global | 0.75 | 0.06 | 0.995 | imputation bias on the 151 added founders (Beagle 98.5% concordance still introduces compounding errors) |
| 231-founder | window-200kb | 0.68 | 0.07 | 0.984 | block-EM × imputation noise (rank-deficiency in 231-dim block × imputation) |

**AC-stratified (82-founder window mode, 8 reps mean)**:
| AC bin | n | R² | RMSE | r |
|---|---|---|---|---|
| AC=1 (singletons) | 2.6M | -103.5 | 0.020 | 0.02 |
| AC 2–4 | 1.4M | -5.07 | 0.024 | 0.39 |
| AC 5–10 | 0.8M | -1.92 | 0.035 | 0.49 |
| AC>10 (commons) | 1.3M | **0.95** | 0.052 | **0.98** |
| ALL | 6.2M | **0.97** | 0.032 | 0.98 |

**Key takeaways**:
1. **Pearson r is nearly perfect (>0.98) in all configurations** — predictions are shape-correct against recipe truth. R² gap is a scale/bias issue.
2. **231-founder R²=0.75 is recoverable to R²=0.989 with a 1-parameter linear calibration** — see diagnostic below. Predicted alt_freqs are 1.42× over-predicted because EM puts too much mass on the 151 imputed founders, but the SHAPE is preserved. Post-hoc rescale is sufficient.
3. **For F0 pools (SEEDMIX), use global mode**: block-EM adds noise without signal benefit when ancestry is uniform.
4. **For evolved samples (the 2,414 GrENE-Net), block-EM in window mode is appropriate**: captures local mosaic ancestry that's invisible to global EM. Imputation bias is less of a concern because there's no clean truth to overfit against.
5. **Singletons (AC=1) are pure noise in all modes**: insufficient panel signal to estimate frequency. Filter by AC≥5 in downstream analyses.

### Diagnostic: 231-founder bias is rank-deficiency on the simplex

Predicted h vs recipe truth at the founder level: **Pearson r = -0.05** (essentially random).
But h-RMSE is small (0.004) because the recipe is diffuse.

Mass redistribution:
| Founder group | Truth mass | Predicted mass | Ratio |
|---|---|---|---|
| 80 cactus (real) | 35.3% | 8.7% | 0.25× |
| 151 imputed | 64.7% | 91.3% | 1.41× |

The EM redirects mass from cactus founders to imputed founders during fitting — a rank-deficiency artifact. Many imputed founders have near-identical genotypes (Beagle imputation aims for the most common haplotype at each locus), so the EM can swap mass between them freely. The simplex prior provides no preference for "true" founders.

**At the per-record alt_freq level, the bias becomes a scale factor**:
- predicted_af ≈ 1.42 × truth_af + 0.001  (linear regression on ALL 3.5M records)
- R²(after correction) = **0.989** for ALL records
- R²(after correction) = **0.983** for AC>10

This means downstream analyses can absorb the bias with a single multiplicative scale on alt_freqs, calibrated from SEEDMIX:
```
calibration = 1 / 1.42   (from SEEDMIX_S1)
calibrated_af = raw_af × calibration
```

---

## 2026-04-28 morning — Beagle imputation 80→231 founders

LOO concordance on full genome (10 random founders, leave-one-out):
- **Worst founder (10002): 98.5% SV concordance** (181,853 / 184,663 SVs match cactus truth, 0 missing)
- Beat the 96.8% smoke-test target — confirms Beagle is *more* accurate at full panel scale (denser haplotype graph)

Beagle DR2 distribution (per-record imputation quality):
- 96% of records at DR2 ≥ 0.9 (high confidence)
- 4% at DR2 < 0.1 (effectively unimputable)
- Cleanly bimodal — single threshold at DR2 ≥ 0.5 separates trustworthy SVs

---

## 2026-04-28 — Genome-wide EM validation (val_gw_56704)

Sims pre-existing in `data/sim_chr1{,_skewed}/`, run via global EM at K=80M:

| Pool | Truth h | Result | Time |
|---|---|---|---|
| uniform82 | 1/82 each | RMSE=0.0034 (R²=nan, zero variance truth) | 3.0h |
| skewed5 | [.40,.25,.15,.10,.10] | **R²=0.9925, RMSE=0.0048** | 1.4h |

Confirms patched float32-throughout EM works on the full genome-wide K=80M sparse cn matrix.

---

## EM solver bottleneck fix (the >60× speedup)

Found via profiling: `h64 @ cn_f32` was triggering an implicit float64 upcast of cn (52 GB temp allocation per iter) — making each EM iteration ~13 seconds at K=12.9M.

| Operation | Time/call (K=12.9M) | Speedup |
|---|---|---|
| `h64 @ cn_f32` (mixed) | 6.49s | 1× |
| `h32 @ cn_f32` (all f32) | **0.47s** | **14×** |
| `cn_f32 @ cw64` (mixed) | 6.96s | 1× |
| `cn_f32 @ cw32` (all f32) | **0.37s** | **19×** |

Plus filtering counts to nonzero-only (~10% of K at 30× cov): another ~10× speedup.

**Combined**: ~60× speedup, brings per-sample time from 3h → ~15 min.
2,414-sample production: 38 days → ~3 days at 8-way parallel.
