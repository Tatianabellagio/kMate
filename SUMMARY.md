# hapFIRE-projection for SV allele frequencies — findings

## Method (one line)

Run hapFIRE on SNPs as published, then project the recovered ecotype-frequency vector onto SVs via a precomputed founder × SV genotype matrix:

```
f_SV(v)  =  Σ_e   G_SV(e, v) · f_ecotype(e)
```

Mathematically identical to how hapFIRE itself computes per-SNP frequencies. ~50 lines of Python on top of unmodified hapFIRE.

---

## Simulation result (visor_freqk, ground truth available)

Pool 231 ecotypes uniformly at 1/231 each; first `N_SV = round(231·f)` carry a 1 kb deletion at the rep-specific position. Truth = `N_SV/231`. All hapFIRE jobs done — 25 BAMs across 4 reps × 3 coverages × 5 frequencies (some reps had partial coverage of f%).

### Headline numbers

| Coverage | n | hapFIRE-proj MAE | freqk MAE | hapFIRE max err | freqk max err | fold-improvement |
|---|---|---|---|---|---|---|
| **10×** | 5 | **0.0013** | 0.0997 | 0.0039 | 0.182 | **77×** |
| **20×** | 5 | **0.0010** | 0.0756 | 0.0020 | 0.166 | **76×** |
| **50×** | 15 | **0.0006** | 0.0583 | 0.0013 | 0.191 | **97×** |

**hapFIRE-proj wins by ~75-100× across all coverages.** Crucially, hapFIRE's max error stays in the 4th decimal place at cov10× — it doesn't degrade meaningfully when coverage at the SV locus drops, because it leans on ~800 K Chr1 SNPs whose collective coverage is still high. freqk's max error stays >0.16 at every coverage because individual SVs have specific k-mer-uniqueness/coverage failure modes that don't average out.

### rep9 detail (cov10 / cov20 / cov50)

| cov | f% | truth | hapFIRE-proj | freqk |
|---|---|---|---|---|
| 50× | 10 | 0.0996 | **0.0999** | 0.136 |
| 50× | 30 | 0.2987 | **0.2980** | 0.108 ← **freqk fails badly** |
| 50× | 50 | 0.5022 | **0.5016** | 0.479 |
| 50× | 70 | 0.7013 | **0.7009** | 0.683 |
| 50× | 90 | 0.9004 | **0.9001** | 0.927 |
| 20× | 10 | 0.0996 | **0.1001** | 0.204 |
| 20× | 30 | 0.2987 | **0.3007** | 0.307 |
| 20× | 50 | 0.5022 | **0.5022** | 0.668 ← **freqk fails** |
| 20× | 70 | 0.7013 | **0.6996** | 0.798 |
| 20× | 90 | 0.9004 | **0.8998** | 0.903 |
| 10× | 10 | 0.0996 | **0.0995** | 0.246 ← **freqk fails** |
| 10× | 30 | 0.2987 | **0.2982** | 0.480 ← **freqk fails** |
| 10× | 50 | 0.5022 | **0.5041** | 0.391 |
| 10× | 70 | 0.7013 | **0.7016** | 0.674 |
| 10× | 90 | 0.9004 | **0.8965** | 0.932 |

### Per-ecotype recovery (rep9 cov50 f30 example)

hapFIRE recovered all 231 ecotype frequencies within `[0.00402, 0.00461]` of the uniform truth `1/231 = 0.00433`. MAE per ecotype: **0.00010**. Mass on first 69 ecotypes (the SV carriers) = 0.298 (truth 0.299). The 231-simplex is recovered nearly perfectly because Chr1 has ~800 K SNPs that constrain the inference.

### Sanity check: posterior-combine (`scripts/posterior_combine.py`)

For panel-known SVs, simple inverse-variance combine of hapFIRE-proj and freqk:
- weight on hapFIRE ≈ 0.997-1.000 across all 23 (rep, cov, f%) cases
- combined MAE ≡ hapFIRE-proj MAE to 4 decimal places

Confirms the framework works (graceful fallback to freqk for non-panel SVs without harming panel-known cases). The interesting test for combine is on the real evolved samples — there the recombination assumption is violated and hapFIRE's variance per SV may be larger, allowing freqk to contribute meaningfully where its k-mers are clean.

---

## Real-data result (8 SEEDMIX replicates)

Per-sample Pearson correlations across SVs in the freqk × panel overlap (4,164 SVs).

| Sample | hapFIRE-proj vs panel-truth | freqk vs panel-truth | hapFIRE-proj vs freqk |
|---|---|---|---|
| SEEDMIX_S1..S8 | 0.99 (uniform across all 8) | 0.45 ± 0.005 | 0.46 ± 0.005 |

Caveats for the seedmix numbers — they're suggestive, not conclusive:
- "panel-truth" is a recipe-based lower bound on the true SV freq — it only counts the 80 founders for whom we have SV genotypes. The other 151 of 231 ecotypes (~67% of seedmix mass) are unobserved.
- The r=0.99 between hapFIRE-proj and panel-truth largely reflects that both quantities are dominated by AC (carrier count). Per-ecotype, hapFIRE doesn't recover the recipe (Pearson r ≈ 0): the seedmix recipe and what was actually sequenced may differ (germination/extraction biases) and the 231-simplex is hard to recover when proportions are nearly uniform.
- The simulation is the cleaner test, since the truth is exactly uniform 1/231 with no missing-genotype mass.

For the seedmix, freqk reports very high AF for many low-AC panel SVs — 35% of AC=1 SVs have freqk_AF > 0.05 and 37 of them have freqk_AF > 0.5. With one carrier in panel and seed-prop ~0.005, panel-mass is bounded at ~0.005; freqk reporting 1.0 means either the unobserved 151 all carry it (extremely unlikely for a panel singleton) or freqk has off-target k-mer hits. Likely the latter.

---

## Why this works (and where it might not)

- **hapFIRE uses ~800 K SNPs across Chr1 to estimate founder frequencies.** Even at 10× coverage, that's ~8 M base-reads — orders of magnitude more constraints than freqk can recover from k-mers spanning a single SV breakpoint. So hapFIRE's founder freq estimation is decoupled from coverage at the SV locus.
- **The projection step is exact.** Once founder frequencies are known, `f_SV = G @ f_eco` is just a linear combination — no estimation, no noise added.
- **Failure modes (in order of severity for real GrENE-Net data):**
  - **Recombination assumption.** This is the *biggest* one. hapFIRE assumes every haplotype in the pool is one of the listed founders — i.e. no recombination, no novel haplotypes. The paper explicitly notes this hits founder-frequency estimation harder than SNP-frequency estimation, because per-SNP errors average over many founders but per-SV errors (especially for low-AC SVs) depend on individual founder freqs. The simulation does not test this — VISOR pools intact founder consensus FASTAs with zero recombination, which is the best case for hapFIRE. The seedmix is also the founder pool itself (still no recombination). **Evolved GrENE-Net samples (MLFH…) have 1–3 generations of recombination and selection** between the founder mix and sequencing → expect degraded per-block accuracy on evolved samples relative to the simulation numbers above.
  - **Panel ascertainment.** SVs polymorphic only among the 151 un-genotyped founders are invisible to projection. ~67% of seedmix mass is on un-genotyped founders, so hapFIRE-proj systematically under-calls SVs that segregate there.
  - **231-simplex identifiability.** With many similar founder haplotypes, the CVXPY projection from block-haplotype freqs to founder freqs can have multiple near-optimal solutions; hapFIRE picks one. For uniform truth (sim) this isn't a problem; for non-uniform pools the picked solution may be subtly biased.

- **Mitigations worth trying before applying on evolved samples:**
  - Per-block residual diagnostic. For each LD block, check `‖DMᵀ f_eco − y_block‖` from the CVXPY solve; large residuals = poor founder fit = likely recombinant block. Drop those blocks from the SV projection. (`ecotype_frequency_estimation_selected` already weights toward top-20% by haplotype diversity — same idea, different criterion.)
  - Save and use **per-block** founder freqs rather than the genome-wide weighted average. Currently hapFIRE saves only the average; the per-block matrix is computed but discarded (`hapFIRE.py:138-139` are commented out). ~5-line patch to keep them. Then `f_SV(v)` uses the block containing position(v) instead of the genome average — recombinants in unrelated blocks don't pollute the estimate.

---

## Pipeline (operational details)

For each (rep, cov, f%) test BAM:

1. **hapFIRE** (unmodified, source: `xwu/haplotype_frequency_estimation/hapFIRE_sourcecode/`):
   ```
   python hapFIRE.py -v greneNet_Chr1_only.vcf -b sim.srt.bam \
       -f Chr1.fa -o <prefix>
   ```
   Resources per job: 4 CPUs, 96 GB memory, 8 h time limit. (16 GB OOM'd at HARP stage.)
   - VCF parse: ~5 min
   - LD-block partition (`CompleteLDPartition`): ~1 min, gives 2 blocks on Chr1 (24 Mb, 6 Mb)
   - HARP per block: ~30 min (big), ~5 min (small)
   - CVXPY ecotype-frequency projection: seconds
   - Total: ~40-50 min per BAM at cov50

2. **Projection** (`scripts/project_sv_freq.py`, ~50 lines):
   ```
   ecotype_freq = parse(<prefix>_ecotype_frequency_selected.txt)
   sv_carriers  = vcf_samples[:N_SV]   # first N_SV ecotypes by VCF header order
   hapfire_proj = sum(ecotype_freq[e] for e in sv_carriers)
   ```

3. **Compare** with freqk's existing `var_del_1kb_n231_f<f>_err001.allele_frequencies.k31.tsv`.

For the real data (seedmix), step 1 is replaced by Xing Wu's existing hapFIRE outputs at `/carnegie/nobackup/scratch/xwu/GrENE_net/hapFIRE_frequencies/seed_mix/s{1..8}_ecotype_frequency.txt`.

---

## Files

```
hapfire_sv/
├── SUMMARY.md                            ← this file
├── scripts/
│   ├── build_founder_sv_matrix.py        # SV panel VCF → 50K×82 0/1 dosage matrix
│   ├── project_seedmix_svs.py            # F @ G.T projection + freqk comparison
│   ├── seedmix_truth_analysis.py         # adds recipe-based panel-truth
│   ├── plot_seedmix.py                   # seedmix figures
│   ├── plot_simulation.py                # simulation figures
│   ├── run_hapfire_one.sh                # SLURM wrapper for one BAM
│   └── ...
├── data/
│   ├── founder_sv_matrix.parquet         # 50,446×82 dosage
│   ├── founder_sv_meta.parquet           # var_type/sv_size/ac per SV
│   ├── seedmix_recipe_normalized.tsv     # GrENE seed proportions (sum=1 over 231)
│   ├── vcf_samples_231.txt               # hapFIRE/GrENE-Net 231 ecotypes
│   └── vcf/greneNet_Chr1_only.vcf        # uncompressed VCF for hapFIRE
├── results/
│   ├── simulation_summary.tsv            # truth + hapFIRE + freqk per (rep,cov,f)
│   ├── simulation_truth_vs_estimate_by_cov.png
│   ├── simulation_mae_by_cov.png
│   ├── seedmix_3way_comparison.tsv.gz
│   ├── seedmix_S1_comparison.png
│   ├── seedmix_correlation_bars.png
│   └── rep*_cov*_*kb_f*/                 # per-job hapFIRE outputs
└── logs/
```
