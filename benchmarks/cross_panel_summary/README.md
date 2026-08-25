# benchmarks/cross_panel_summary/

The **one table that spans every panel benchmark**, plus the figures built from it and
the tooling that produces both.

Individual benchmarks (`p231/`, `p80/`, …) each score themselves in their own
`results/`. This dir is where those are pooled into a single comparable table, so you
can ask "how does kMate do across panel × regime × estimator mode" without opening
five directories.

Was called `results/` until 2026-08-25 — a name that said nothing at the top level.

## Contents

    benchmark_table.tsv    the master table (544 rows)
    plots/                 figures built from it
    scripts/               the tooling that writes both

`benchmark_table.tsv` columns:

    tool panel mode coverage n_founders generation mating selection seed
    var_class n MAE RMSE R2 pearson_r bias outlier truth_col est_file

Spans **panels** {p80, p231} × **modes** {global, chrom, block} × the regime matrix
(coverage, generations, mating, selection, seed) × var_class {SNP, indel, SV}.

⚠️ `mode` mixes estimator eras: `global` and `block` rows predate the 2026-07-07
estimator refresh, `chrom` rows are post-refresh. **Do not compare a `global` row
against a `chrom` row and read the difference as a regime effect** — it is confounded
with the estimator change. `est_file` names the run each row came from; check it before
drawing a cross-mode conclusion.

## Scripts

| script | what it does |
|---|---|
| `build_benchmark_table.py` | append scored rows to `benchmark_table.tsv` (one invocation per run) |
| `build_benchmark_figure.py` | RMSE barplots from that table → `plots/benchmark_panel_rmse*.png` |
| `plot_scenario_grids.py` | the per-panel scenario grids, written into each panel's own `<panel>/results/plots/` |
| `plot_scenario_grids_ldr01.py` | same, for the `--unit ld --ld-r2 0.1` control arm |
| `_accuracy_panel.py` | shared density-scatter aesthetic imported by both grid scripts |

The grid scripts write *into the panel dirs*, not here — they are cross-panel **tooling**
whose outputs belong with the panel they describe. They live here because they are
shared by p231 and p80 and would otherwise be duplicated or stranded at the root.

Run in the `basic` env (matplotlib hangs in `plotting`):

```bash
PY=/global/home/users/tbellg/miniforge3/envs/basic/bin/python
$PY scripts/build_benchmark_figure.py               # defaults resolve to ../benchmark_table.tsv
$PY scripts/plot_scenario_grids.py
```
