# orphaned_aggregators/ (archived 2026-08-25)

The former `benchmarks/scripts/` — four aggregation scripts with **zero references**
anywhere in the repo (each matched only itself on a `git grep`).

They are orphaned because the drivers whose output they aggregate were themselves
archived to `src/archive/`:

| script | aggregated output of | driver now at |
|---|---|---|
| `aggregate_results.py` | `batch_runner` per-sample TSVs → wide matrix | `src/archive/batch_runner.py` |
| `aggregate_seedmix_validation.py` | `validate_seedmix_recipe` metrics over the 8 SEEDMIX reps | `src/archive/validate_seedmix_recipe.py` |
| `analyze_seedmix_v3.py` | the v3-era SEEDMIX analysis | superseded panel (v3qc lineage, `PIPELINE_STATE.md` §4) |
| `build_cactus_em_vs_hapfire_site04.py` | a single-site `cactus_em` vs hapFIRE comparison | pre-kMate naming; shared-panel design |

Note the last one is also a shared-panel kMate-vs-hapFIRE comparison, i.e. the same
design retired in `../shared_panel_hapfire_comparison_invalid/`.

Current SEEDMIX validation lives at `scripts/validate_seedmix_vs_hapfire.py` (repo
root `scripts/`, not `benchmarks/`), writing to
`analysis/panel_qc/seedmix_validation_vs_hapfire.tsv`.
