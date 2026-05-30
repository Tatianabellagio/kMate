# benchmarks/

End-to-end accuracy benchmarks for the kMate estimator: simulate pool-seq reads from
a known founder mixture, run the full pipeline, and compare estimated allele
frequencies against the simulation truth.

Both benchmarks drive the shared pool-seq simulation framework in
[`../sims/`](../sims/) (mosaic construction → VISOR pooled read simulation →
truth tables; see `../sims/README.md` and `../docs/SIMULATIONS_METHODS.md`) and
the estimator in `../src/`.

| dir | panel | what it isolates |
|---|---|---|
| `p80/` | **homogeneous** 80-cactus-founder panel (all long-read assemblies, no PanGenie founders) | Control: with no cactus-vs-PanGenie k-mer asymmetry, does the +41% h-bias and the off-diagonal streak vanish? Established the **panel-conditional caveat** that ω_k=1/m_b slightly under-performs uniform EM on a balanced panel. |
| `p231/` | **full** production 231-founder panel (78 cactus + 153 PG) | Headline benchmark on the real heterogeneous panel kMate ships against, across the regime matrix (coverage, generations, recombination, skew). |

Each benchmark dir follows the same numbered-stage layout:

```
NN_*.sh        # build panel inputs → simulate → run estimator → evaluate
scripts/       # stage scripts + sim-mosaic builders
sims/          # per-panel sim outputs (gitignored; regeneratable)
results/       # FINAL_RESULTS_*.ipynb + summary TSVs + plots
```

See each subdir's `README.md` for the stage-by-stage walk-through. Results notebooks
(`results/FINAL_RESULTS_*.ipynb`) hold the regime-by-method MAE tables and figures.

> Historical note: `p80/` succeeds the earlier `control_p82/` control (dropped after the
> 2026-05-16 panel-QC decisions removed 2 cactus assemblies).
