# benchmarks/archive/

Superseded benchmark artifacts, kept for provenance (do not treat as current).

- `benchmark_summary_global_vs_block_PREFIX_STALE.csv` — pre-fix (2026-06-20)
  global-vs-block summary. Its g0 global-vs-block ordering is now REVERSED by the
  refreshed `../benchmark_table.tsv` (2026-07-07), because the corrected estimator
  (per-founder M-step normalization + `--kmer-weight uniform`, superseding `inv_mb`)
  changes the global result. Numbers are stale; use `../benchmark_table.tsv`.

- `speed_benchmark_abandoned/` (2026-07-10) — an earlier, incomplete multi-tool
  (kMate/hapFIRE/vg) runtime harness. Only kMate+hapFIRE were ever run (p80,
  2 coverage points); vg was never executed and no results table was produced.
  Superseded by `../speed_vs_hapfire/`, the maintained, documented speed
  comparison. Unreferenced elsewhere in the repo — safe to have moved.

- `run_benchmark_pool.sbatch`, `submit_benchmark_matrix.sh` (2026-07-10) — the
  pre-`--unit`-unification pool-benchmark harness (comments still say "count
  k-mers ONCE, run GLOBAL + BLOCK(dynld_K500)"). Superseded by the `07f`/`07g`
  runners under `p231/scripts/` and `p80/scripts/`. Unreferenced elsewhere.

**Not archived despite looking stale — still load-bearing, needs a follow-up
decision before moving:**
- `../benchmark_runs/` (~11 GB, pre-unification per-pool AF tables) — still read
  by `accuracy_vs_competitors/scripts/{run_hapfire_pool.sbatch,
  plot_block_nofallback_af.py, stratify_vg_multiplicity.py,
  score_all_competitors.sh}`.
- `../ldblock_window_test/` (~4 GB, pre-unification LD-clique coarseness sweep)
  — still read by `analysis/grenenet_gea/{block_kmer_coverage.py,
  _build_blocks_units_nb.py}` and `accuracy_vs_competitors/scripts/
  {run_ldblock_window.sbatch, plot_blockpanel.py}`.
- `../p231/results/kmate_global_*` / `../p80/results/cactus_em_*` (old
  weighting-arm result dirs) — intentionally retained per
  `ROADMAP_GLOBAL_REFRESH.md` ("keep old dirs for the old↔new delta column")
  and actively read by `plot_scenario_grids.py` and
  `p80/results/{FINAL_RESULTS_cov10_p80.ipynb,FILT2_MB_VS_UNIFORM_p80.ipynb}`.
- `../localonly_p231/` — cited as evidence in
  `analysis/grenenet_gea/GLOBAL_MODE_DECISION.md` and
  `_build_floor_derivation_nb.py`.
- `../h_uncertainty/` — closed-out investigation, but not yet confirmed fully
  orphaned; its module concept is referenced from core docs (`ALGORITHM.md`,
  `src/README.md`), so double-check before moving the benchmark dir itself.
