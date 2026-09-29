# common — shared inputs

Data preparation used by more than one results section. Nothing here answers a
question; everything downstream reads it.

## The chain

```
kMate per-sample AF TSVs  (rerun_kfw_hb/)
      |  build_af_store.py            -> results/af_store/      (per-sample NPY store)
      +--build_gen_matrices.py        -> results/gen_matrices/  (per-GENERATION per-SAMPLE, float16)
      +--build_gen_matrix.py          -> results/gen_matrices/  (per-GENERATION per-TIMEPOINT stacked, uint16)
      +--build_pool_matrix.py         -> results/pool_matrices/ (flower-weighted site_gen_plot pools)
      +--build_p0.py                  -> founding gen-0 p0 from the 8 SEEDMIX reps
      +--build_sample_h_cache.py      -> results/fitness/       (per-sample founder h, 231-vector)
```

> ⚠ **`build_gen_matrices.py` and `build_gen_matrix.py` are different scripts**,
> not a typo pair. Plural = per-generation per-**sample** AF matrices (float16).
> Singular = per-generation per-**timepoint** stacked matrices (uint16). Check
> which one you mean.

| script | writes | driver |
|---|---|---|
| `build_af_store.py` | `results/af_store/` — the foundation store `lib.py` reads | `build_af_store_array.sh` |
| `build_gen_matrices.py` | `results/gen_matrices/` per-sample AF by class | `build_gen_matrices_array.sh` |
| `build_gen_matrix.py` | per-timepoint stacked matrices | — |
| `build_pool_matrix.py` | `results/pool_matrices/` — the most-read derived product (~25 consumers across r1/r2) | — |
| `build_p0.py` | founding p0, aligned to the AF store | — |
| `build_sample_h_cache.py` | `results/fitness/` founder-h cache | `run_sample_h_cache.sbatch` |
| — | group_means cache | `rebuild_group_means.sbatch` |

`rerun_kfw_hb/` (1.6 TB) is the production AF input: the `--unit chrom` +
full-panel Kf_w rerun. See its own `RERUN_MAP.md`.
