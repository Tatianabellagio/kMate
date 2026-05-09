# cactus_em — algorithm explanation

End-to-end description of how cactus_em estimates SNP and SV allele frequencies in pool-seq samples. Companion docs: `MODEL_SPEC.md` (math + early design notes), `BACKGROUND.md` (why we built it), `SUMMARY.md` (real-data results), `SIMULATION_RESULTS.ipynb` (full benchmark).

## The problem

Given pool-seq reads from a population that's a mixture of N known founder genomes (the "panel"), recover the allele frequency at every variant in the panel.

For GrENE-Net specifically:
- 231 *A. thaliana* founders in the panel — 80 with long-read assemblies (cactus pangenome), 151 imputed via Beagle from short-read SNPs
- ~3.4 M variants per chrom panel (SNPs + INS + DEL + SVs)
- Pool-seq short reads from one sample (e.g. SEEDMIX_S1, or one of the 2,415 evolved GrENE-Net samples)

Output: per-VCF-record alt allele frequency in the pool, one number per (chrom, pos, ref, alt) per sample.

## Two key matrices, computed once per panel

### `cn_full` — founder × k-mer (the EM evidence)

For each pangenome bubble, **PanGenie** produces a list of unique-to-this-bubble k-mers (post var-dedup + ref-dedup, capped at 16/allele). For each k-mer `k` in bubble `v`, we check which of the N founders' haplotypes through bubble `v` contain `k`.

```
cn_full[f, k] = 1   if founder f carries k-mer k
                0   otherwise
```

Genome-wide: ~80 M unique k-mers × 231 founders, sparse int8 CSR. ~30-50 GB in RAM.

Built by `poolfreq/src/build_kmer_cn.py` from PanGenie's `kmers.tsv.gz` + the panel VCF + reference FASTA. One-time cost: ~5-8 h SLURM.

### `cn_var` — founder × variant (the projection target)

One column per VCF record (per ALT allele after `bcftools norm -m -`), N founder rows.

```
cn_var[f, r] = 1   if founder f's GT at variant r is the ALT (any GT > 0)
               0   otherwise
```

Genome-wide: ~3.4 M variants × 231 founders.

Built by `poolfreq/src/build_cn_var.py` from the panel VCF GT field. ~30 min.

For the 80 cactus founders, GTs come from the cactus pangenome VCF directly. For the 151 imputed founders, GTs come from Beagle imputation against the GrENE-Net SNP haplotypes — see `imputation/IMPUTATION_PLAN.md`.

## Per-sample run (the heavy step)

**Input:** paired FASTQs from a pool-seq sample.

**Step 1 — count k-mers in reads.** Use jellyfish with the genome-wide k-mer set as the query. Output: `c[k]` = number of read-hits to each of the 80M k-mers in cn_full's index. Discard k-mers not in the index — irrelevant for downstream.

**Step 2 — estimate sample coverage `λ`.** From total k-mer hits, derive an expected per-k-mer rate (≈ read coverage / 1 since each k-mer position is hit once per copy).

**Step 3 — EM solve for founder frequencies `h`.** This is the core innovation. The generative model is

```
c[k] ~ Poisson( λ · h^T · cn_full[:, k] )
```

for each k-mer `k`, with `h` on the founder simplex (`h ≥ 0`, `Σh = 1`). Equivalently: treat each unit of count at k-mer `k` as a latent draw from one of the founders carrying that k-mer, with mixture proportion `h`.

The simplex-constrained MLE is found by the standard multiplicative EM update (Lee–Seung-style; identical in form to multinomial-mixture EM and to NMF with simplex normalization):

```
μ_k(h)    =  h^T · cn_full[:, k]                          (E-step denominator)
h_new[f]  ∝  h[f] · Σ_k  cn_full[f, k] · c[k] / μ_k(h)    (M-step, then renormalize Σh = 1)
```

`λ` cancels in the M-step ratio (it's constant across `k` under the mixture interpretation), so the solver doesn't actually consume the coverage estimate from Step 2 — coverage is computed for sanity-checking only.

The iteration monotonically increases the likelihood and converges in tens of iterations from a uniform start. The float32-throughout implementation runs in ~5–15 min per sample on the 80M-k-mer matrix. See `poolfreq/src/em_solver.py`.

**Why pool k-mers across the whole genome?** Each individual k-mer is sparse (most founders don't carry it), but jointly the 80M k-mers over-determine the 231-vector h. This is why we get clean h estimates even at low pool-seq coverage — the joint signal is huge.

**Step 4 — project h to per-record AFs.** Once h is known, the alt allele frequency at variant `r` is just:

```
pred_af[r] = h · cn_var[:, r]
```

Vector-matrix multiply. Linear, exact, no estimation noise added at this step. Output: per-VCF-record alt freq. See `poolfreq/src/per_sample_driver.py:run_one_sample()`.

## Optional: per-block (window) EM

For evolved samples with mosaic ancestry (1-3 generations of recombination + selection), the assumption "every haplotype is a clean founder" is violated. The mitigation in `poolfreq/src/block_em.py`:

- Partition the genome into ~200 kb windows (`--window-bp 200000`, default)
- Solve the EM per-window → per-window `h_w`
- For each variant `r` in window `w`, project `pred_af = h_w · cn_var[:, r]`

Window mode trades global accuracy for local resolution — handles recombinants who switch founders across windows. This is the production setting for evolved samples (`--block-mode window`). For homogeneous F0 pools (SEEDMIX), `--block-mode global` is correct.

## What's structurally different from freqk

| | cactus_em | freqk |
|---|---|---|
| **Inference** | joint EM across all 80M k-mers under simplex constraint | per-bubble independent counts |
| **Founder modelling** | explicit (h is per-founder) | none (no founder concept) |
| **Cross-bubble pooling** | yes — bubbles with sparse k-mers borrow strength from bubbles with rich k-mers via h | no |
| **Per-record output** | h × cn_var — exact projection | per-bubble alt-count / total-count |

The joint EM is the structural differentiator. freqk treats each bubble independently; cactus_em pools k-mer evidence across all bubbles via the founder-frequency vector h. On panels where individual bubbles have sparse k-mer evidence (after PanGenie's cap + dedup), cactus_em's pooling is what buys the extra accuracy.

## What's structurally similar to hapFIRE

| | cactus_em | hapFIRE |
|---|---|---|
| **Joint inference** | EM on all k-mers | HARP per LD-block, then CVXPY founder projection |
| **Founder modelling** | h on simplex | haplotype freqs → founder freqs via CVXPY |
| **Per-record output** | h @ cn_var | haplotype_freq @ haplotype_to_SNP_matrix |

That's why cactus_em and hapFIRE agree at R² = 0.996 on real SEEDMIX_S1 SNPs — they're solving the same inverse problem with the same panel info, just via different inference routes (k-mer EM vs HARP+CVXPY).

cactus_em adds SVs (which hapFIRE can't do, because HARP is per-base SNP-only), and uses the cactus pangenome graph k-mers (richer than HARP's per-SNP windows).

## Files in the pipeline

```
poolfreq/
├── src/
│   ├── build_kmer_cn.py       # builds cn_full from PanGenie kmers.tsv + VCF
│   ├── build_cn_var.py        # builds cn_var from VCF GT field
│   ├── em_solver.py           # core EM (float32-throughout)
│   ├── kmer_count.py          # jellyfish wrapper for read k-mer counts
│   ├── block_em.py            # window-mode partition + per-block EM
│   ├── per_sample_driver.py   # one sample → per-record AF TSV
│   └── batch_runner.py        # multi-sample parallel
└── data/
    ├── cn_full_231_v2/cn_Chr{1..5}.cn.npz   # 231 × 80M matrix
    └── cn_var_231_v2.cn_var.npz             # 231 × 3.4M matrix
```

## Tradeoffs

**Strengths:**
- Joint k-mer EM gives ~1% per-record MAE (matches hapFIRE's accuracy on SNPs)
- Single pipeline handles SNPs + INS + DEL + SVs (no separate runs)
- Coverage-saturated below 1× — the joint-evidence is over-determined

**Costs:**
- 100-200 GB per-sample memory (vs freqk's 16-32 GB) — see `HANDOFF.md` for the production-scale memory note
- 30-45 min per-sample runtime (vs freqk's 10-20 min)
- Requires pre-built cn_full and cn_var (one-time, ~10h SLURM)
- Founder genotypes must be available (from long-read assemblies + Beagle imputation for missing founders)

## Empirical performance

See `SIMULATION_RESULTS.ipynb` for the full benchmark suite. Headlines:

**Simulation (cov 10×, 82-founder skewed pool, n=2.69M polymorphic records):**
- cactus_em (231-cn): R² = 0.998, MAE = 0.008, slope = 1.02, intercept = -0.001
- freqk on cactus pangenome: R² ≈ 0.43 on biallelic single-record windows
- Coverage saturation: same MAE at 1× as at 30×

**Real SEEDMIX_S1 SNPs vs hapFIRE (n=3.24M):**
- cactus_em: R² = 0.996, RMSE = 0.013, MAE = 0.008 (slope 1.013, basically reproduces hapFIRE)
- freqk: R² = 0.94, RMSE = 0.056, MAE = 0.034 (per `freqk_gr/results/`)

**Real evolved samples (site04, 57 samples)** — running now via SLURM array `58866`. Will give cactus_em-vs-hapFIRE-vs-freqk on real recombinant pools.

## Open questions / loose ends

1. **Recombinant-pool simulation** (HANDOFF.md item #2) — all sims so far pool intact founder FASTAs. Real evolved samples have recombination. The window-mode EM should handle this but is only validated on real SEEDMIX (no truth) — a controlled mosaic-FASTA sim is the right way to verify.
2. **Memory mitigation** (HANDOFF.md memory note) — cn_full mmap from disk would cut RAM 5-10× and unlock more concurrent SLURM tasks. Worth doing before scaling to >2,500 samples.
3. **SV recovery on noisy multi-allelic records** (Tier 6 of the simulation notebook) — clean per-allele 1:1 alignment in the v2 aggregator showed SVs at R² ≥ 0.988 once the aggregation bug was fixed; the open question is whether real evolved samples preserve this.
4. **Replicate consistency** — only SEEDMIX_S1 has been run through the v2 cactus_em pipeline. Running S2-S8 would tighten the hapFIRE comparison from 1 rep to 8.
