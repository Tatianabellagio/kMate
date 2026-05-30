# sims — pool-seq simulation framework

Self-contained home for the simulations that benchmark `kMate`. The full
methods write-up (regime matrix, parameters, citations, validation) is
**`docs/SIMULATIONS_METHODS.md`** — this README is the code/layout map.

> Decoupled from the former external `visor_freqk` repo (2026-05-30): the
> framework code that the benchmark actually depends on now lives here and is
> tracked in the kmate repo. There is no longer a nested git repo under `sims/`.

## Layout

```
sims/
  scripts/
    make_recomb_mosaics.py     canonical mosaic-FASTA generator (Stage 1: founder draw,
                               optional G-generation recombination at hotspot or random
                               crossovers, stitch consensus FASTAs into per-individual
                               haploid mosaics). Supports --source-weights and
                               --gen0-no-replace (SEEDMIX-mimicry balanced allocation).
    compute_recomb_truth.py    per-record realized-pool AF truth from ancestry tracks
                               (MAR projection through var_pa / var_called; §4 of the
                               methods doc).
    build_g0_uniform_sim.py    g0 uniform-pool sim generator helper.
  data/
    hapfire_block_index_chr1.npz   LD-block boundaries (BigLD on the GrENE-Net 231-panel
                               SNP set; Kim et al. 2018) used as hotspot crossover
                               positions for recombinant regimes. Fixed input, gitignored
                               (binary data, 7.5M); regenerate from the 231-panel SNP set
                               via BigLD if absent (see docs/SIMULATIONS_METHODS.md §9).
  data/ logs/ results/         sim outputs (gitignored)
```

## Who calls this

The per-benchmark drivers consume the framework above; they hold panel-specific
config + the VISOR read-sim stage, and call the **canonical** scripts here (no
per-benchmark copies of the mosaic builder anymore):

- `benchmarks/p80/scripts/06_run_sim_p80.sh`, `06b_run_sim_p80_skewed.sh`
  (80-cactus control; uses the LD-block hotspot crossovers from `data/`)
- `benchmarks/p231/scripts/06_run_sim_p231.sh`, `06b_run_sim_p231_skewed.sh`
  (231-founder headline; random crossovers, does not use the block index)

Benchmark-local scoring/validation (`compare_af_vs_truth.py`,
`sim_from_raw_assemblies.py`) stays in `benchmarks/p80/scripts/` — it is
evaluation tooling, not the simulation framework.

## Reproducibility

Python ≥3.10 in conda env `hapfm`; VISOR + samtools/bcftools in env `pang`. All
randomness is seed-controlled — same seed + params reproduce bit-identically.
See `docs/SIMULATIONS_METHODS.md` §10.
