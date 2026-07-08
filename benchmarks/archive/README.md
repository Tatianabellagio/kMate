# benchmarks/archive/

Superseded benchmark artifacts, kept for provenance (do not treat as current).

- `benchmark_summary_global_vs_block_PREFIX_STALE.csv` — pre-fix (2026-06-20)
  global-vs-block summary. Its g0 global-vs-block ordering is now REVERSED by the
  refreshed `../benchmark_table.tsv` (2026-07-07), because the corrected estimator
  (per-founder M-step normalization + `--kmer-weight uniform`, superseding `inv_mb`)
  changes the global result. Numbers are stale; use `../benchmark_table.tsv`.
