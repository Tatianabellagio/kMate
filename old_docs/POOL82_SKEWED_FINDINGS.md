# pool_sweep_82_skewed — findings (2026-05-01)

## TL;DR

On the realistic skewed pool (multinomial 200-individual draw from 82 cactus
founders, seed 42, 9/82 with count=0), running the same reads through both
pipelines:

| panel | n_poly | MAE | r | slope | intercept |
|---|---|---|---|---|---|
| 82cn | 5.81M | 0.0498 | 0.785 | 0.594 | +0.127 |
| **231cn** | **2.69M** | **0.0118** | **0.985** | **1.016** | **+0.004** |

**231-cn (post Tier-1 fix) is essentially unbiased on simulation.** The
1.43× slope artifact reported on SEEDMIX_S1 is gone in `cn_var_231_v2`.
4× lower MAE than 82-cn at every coverage tested.

Coverage doesn't matter — 5x ≡ 30x to 4 decimal places (same flat behaviour
the uniform pool showed). cactus_em is k-mer-saturated on Chr1+ at all
panel-covered loci.

## Why the headline 82-cn looks bad: it's almost entirely COMPLEX-record drag

Stratifying by var_type at cov10:

| panel | var_type | n | MAE | r | slope | intercept |
|---|---|---|---|---|---|---|
| 82cn | SNP | 3.13M | 0.011 | 0.973 | 0.996 | +0.007 |
| 82cn | INS | 752K | 0.067 | 0.698 | 0.920 | +0.072 |
| 82cn | DEL | 818K | 0.041 | 0.801 | 0.947 | +0.043 |
| 82cn | **COMPLEX** | **1.11M** | **0.155** | **0.526** | 0.802 | **+0.174** |
| 231cn | SNP | 2.52M | 0.009 | 0.998 | 1.024 | 0.000 |
| 231cn | INS | 112K | 0.065 | 0.412 | 1.247 | +0.057 |
| 231cn | DEL | 54K | 0.039 | 0.743 | 0.986 | +0.037 |

`cn_var_82` includes 1.1M COMPLEX records (ref_len ≠ alt_len, both > 1)
that `cn_var_231_v2` filters out during the Beagle/F_MISSING pipeline. Those
records have MAE 0.155 and pull the 82-cn headline down. **On SNPs alone the
two pipelines are essentially equivalent (slope ~1.0 both panels).** 231-cn's
edge is on the ~150K SVs it does keep, plus generally cleaner SNPs.

## AF-bin breakdown (231cn at cov10)

| AF bin | n | MAE | r |
|---|---|---|---|
| (0, 0.05) | 1.42M | 0.009 | 0.220 |
| [.05, .10) | 439K | 0.010 | 0.429 |
| [.10, .25) | 405K | 0.014 | 0.855 |
| [.25, .50) | 216K | 0.019 | 0.944 |
| [.50, .75) | 107K | 0.021 | 0.956 |
| [.75, 1) | 100K | 0.022 | 0.982 |

MAE is bounded above by 0.022 across the entire AF spectrum on 231cn. Pearson
r climbs as AF range opens up (low-AF bins compress truth variance, depressing
r mechanically). Production-grade at every AF bucket.

## Loose ends worth probing

1. **231-cn INS slope = 1.247** at cov10 (over-call by 25%). Smaller n
   (112K) than DEL or SNP. Could be a real Beagle-imputation effect for
   insertions specifically (insertions are harder to impute than deletions
   because they introduce sequence not present in the SNP haplotype). Not
   in 82-cn (slope 0.92) — so it's a 231-cn-only artifact.
2. **Recipe-weighted (231-founder) pool** — the production target is
   2,415 evolved samples drawn from a 231-founder seed mix, not an
   82-founder one. The 82-founder pool was used here because all 82
   have long-read assemblies (clean truth). Once we trust the 231 panel,
   doing the same sim with a 231-founder skewed pool tests the
   production-realistic regime.

## Amplification sims in flight (launched 2026-05-01 ~11:20)

| Job | Coverage | Seed | Purpose |
|---|---|---|---|
| 58457 | 10x | 43 | Monte Carlo replicate |
| 58458 | 10x | 44 | Monte Carlo replicate |
| 58459 | 10x | 45 | Monte Carlo replicate |
| 58460 | 1x  | 42 | Saturation floor |
| 58461 | 2x  | 42 | Saturation floor |
| 58462 | 3x  | 42 | Saturation floor |

After these complete, run `scripts/refresh_pool_sweep_82_skewed.sh` to
update `summary.tsv` and re-execute the notebook. The new plots will
populate:
- `pool82s_saturation_curve.png` — full 1-30x range
- `pool82s_cov10_seed_variance.png` — Monte Carlo variance across 4 seeds

## Next batch — recipe-weighted pool (queued, not launched)

`scripts/run_pool_sweep_82_recipe.sh COVERAGE` is ready. Pool composition
matches the SEEDMIX recipe restricted to the 80 cactus founders (the 2
non-overlap cactus founders contribute 0 since they aren't in GrENE-Net
recipe). Total cactus recipe mass = 35.3% of full SEEDMIX; renormalized to
sum=1 over 80 founders.

Launch when the in-flight 6 jobs have cleared:
```
sbatch scripts/run_pool_sweep_82_recipe.sh 5
sbatch scripts/run_pool_sweep_82_recipe.sh 10
sbatch scripts/run_pool_sweep_82_recipe.sh 20
sbatch scripts/run_pool_sweep_82_recipe.sh 30
```

This is the most production-realistic test of the cactus_em pipeline on
sim data (no recombination). Compare side-by-side to the multinomial pool:
recipe is more skewed (max-frac 1.4% vs multinomial-uniform ~1%) and
includes 2 zero-mass founders, so recovery may degrade slightly.

## Files

```
pool_sweep_82_skewed/
├── summary.tsv                                    # all strata, both panels, all 4 covs
├── per_record_cov{5,10,20,30}_{82,231}cn.tsv.gz  # per-record (truth, pred)
└── cov{N}_n200_s42/cactus_em_result_{82,231}cn.tsv  # raw EM output

scripts/aggregate_pool_sweep_82_skewed.py
summaries/pool_sweep_82_skewed.ipynb              # 7 plots, all rendered

plots/pool82s_metrics_vs_coverage.png             # MAE/RMSE/r vs cov, by panel
plots/pool82s_slope_vs_coverage.png               # slope=1 reference line
plots/pool82s_cov10_scatter_{82cn,231cn}.png      # 3-panel hexbin per var_type
plots/pool82s_cov10_mae_by_vartype.png            # bar
plots/pool82s_cov10_mae_by_af.png                 # bar
plots/pool82s_cov10_mae_by_svsize.png             # INS+DEL by size bin
plots/pool82s_cov10_head_to_head.png              # 82-cn predicted vs 231-cn predicted
```
