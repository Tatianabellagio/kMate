# benchmarks/

End-to-end accuracy benchmarks for the kMate estimator: simulate pool-seq reads from
a known founder mixture, run the full pipeline, and compare estimated allele
frequencies against the simulation truth.

**Every simulation and kMate fit behind the paper's kMate figures runs from [`simulations/`](simulations/)**: one
simulator, one fit runner, one run table (`simulations/runs.py`, grouped by figure) and one scorer, whose tables the
benchmark notebook (`analysis/grenenet_selection/notebooks/_build_kmate_benchmarks_nb.py`) draws. The per-benchmark
scripts that ran them before are in `archive/scripts_before_simulations_2026-10-17/`. The directories below keep the
inputs those runs read (founder genomes, k-mer and variant panels under `p231/` and `p80/`), hapFIRE's arm
(`vs_hapfire/`) and the exploratory benchmarks that are not in the paper.

---

## Directory index

Every dir below has its own `README.md` with the question it answers.

**Panel benchmarks (kMate vs simulation truth)**

| dir | what it isolates |
|---|---|
| `simulations/` | **every pool and fit of the paper's kMate figures** (run table, simulator, fit runner, scorer) |
| `p231/` | inputs of the production 231-founder panel benchmarks (founder genomes, k-mer panels) |
| `p80/` | inputs of the long-read (78-founder) panel benchmarks (founder genomes and chains, 78-founder panels) |
| `founder_h_accuracy/` | founder-mixture `h` recovery (most sensitive to estimator internals) |
| `founder_h_imbalance/` | is the private-k-mer filter still needed post-normalization? |
| `founder_h_uncertainty/` | closed-out: `h` error is identifiability bias, not sampling variance |
| `unit_localonly_p231/` | local-only window fitting + the block floor derivation |
| `unit_ldblock_window/` | ⚠️ window-mode block-size sweep — pre-refresh, see its README |
| `realdata_outcross/` | outcrossed real-data application |

**Tool comparisons** — ⚠️ read `SCORING_RULES.md` first

| dir | comparator | status |
|---|---|---|
| `vs_hapfire/` | hapFIRE | ✅ hapFIRE's inputs (GrENE-Net genomes, VCFs, genotype matrix, bwa index) and `ecotype_resolution/`; its pools and runs are in `simulations/` (matched founders, **tool-native panels**) |
| `vs_vg_giraffe/` | vg giraffe | ✅ vg arm live; its hapFIRE arm is retired to `archive/` |

**Infrastructure**

| dir | what it is |
|---|---|
| `kmate_output_store/` | shared per-pool kMate output store (~11 GB) — ⚠️ pre-refresh, see its README |
| `archive/` | superseded and invalidated work, each with a README saying why (`scripts_before_simulations_2026-10-17/`: what `simulations/` replaced, with `MANIFEST.csv`) |

### The one rule worth repeating

**Never compare kMate and hapFIRE on a shared panel.** hapFIRE errors on missing GTs,
so any panel built for it imputes `MISSING -> 0|0`, handing it false REF homozygotes
where kMate's `var_called` mask simply excludes the site. A shared panel is not a
shared *input*. Two benchmarks were built that way and both are archived — see
`archive/shared_panel_hapfire_comparison_invalid/`. Give each tool its native panel;
a shared *truth* is still fine when the quantity is panel-independent (e.g. founder
pool weights).
