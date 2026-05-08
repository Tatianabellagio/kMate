# Archived: superseded cn_var matrices

Moved here on 2026-05-08 from `poolfreq/data/`. Production cn_var is
`poolfreq/data/cn_var_231_v2.{cn_var,meta}.npz` and is unaffected by this move.

**Canonical version-history table lives at the project level:**
- `METHODS_TRIED.md` § 5b — "Panel / cn_var version history" — lists each
  superseded build, the era it served, and why it was replaced
- Project-memory note `project_cn_var_231_v2_is_beagle_imputed.md` clarifies
  the v2 imputation status (intentional, sim-validation only — does NOT
  conflict with the unimputed production VCF decision)

This directory is a candidate for deletion once disk is needed. The
superseded matrices are not referenced by any production code path; a few
historical/smoke-test scripts still mention them by name (`run_plots.py`,
`poolfreq/tests/run_seedmix_*.sh`, archived hapfire-projection notebook), but
those are not on the production critical path.

## What's in this directory

| file | era |
|---|---|
| `cn_var_82.{cn_var,meta}.npz` | hapFIRE-projection era, 82-founder cactus panel |
| `cn_var_82_renamed_to_1001g.*` | same era, ID-renamed |
| `cn_var_231.{cn_var,meta}.npz` (v1) | first 231-founder build (had a carrier-status bug, fixed in v2) |
| `cn_var_decomposed_chr1.*` | Chr1-only smoke test |

## How to restore

```bash
mv archive/cn_var_old/cn_var_*.npz poolfreq/data/
```
