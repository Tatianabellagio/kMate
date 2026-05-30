# kMate source inventory

What's actually used to estimate founder frequencies (`h`) and project to
per-record allele frequencies. Last cleaned: 2026-05-26.

There are exactly **two estimators**: `global` and `window`. Nothing else.

## Active — the estimator (h estimation + AF projection)

| file | role |
|---|---|
| `em_solver.py` | EM core. `solve_em()` — the founder-mixture estimator (E/M multiplicative update; supports the Dirichlet anchor used by window mode). |
| `kmer_count.py` | Jellyfish wrapper. Counts canonical k=31 k-mers in the reads → the count vector fed to the EM. |
| `block_em.py` | Window-mode pieces: `define_windows`, `assign_kmers_to_blocks`, `assign_records_to_blocks`, `solve_em_per_block` (per-window EM + global anchor), `project_blocks_to_records` (per-window → per-record AF, missing-aware). |
| `block_haplotype_em.py` | `smooth_h_across_blocks` — Li–Stephens-style smoothing of per-window `h` (window mode). |
| `per_sample_per_chrom.py` | Production driver. FASTQ/BAM → counts → EM (per chromosome) → AF projection → output TSV. Dispatches `--block-mode {global, window}`. |

### Production recipes

```bash
# global — selfing / inbred / F0 pools (e.g. SEEDMIX)
python per_sample_per_chrom.py \
  --kmer-pa-prefix <kmer_pa>/kmer_pa \
  --var-pa <panel>.var_pa.npz --var-called <panel>.var_called.npz \
  --var-meta <panel>.meta.npz \
  --reads R1.fq R2.fq --sample <name> --out <name>.tsv \
  --threads 8 --chroms Chr1 --block-mode global

# window ("star2") — recombinant pools. The window defaults ARE this recipe,
# so `--block-mode window` alone reproduces it:
#   --window-bp 10000 --global-anchor-weight 0.3 --hmm-smooth-passes 5 --hmm-smooth-alpha 0.5
python per_sample_per_chrom.py [same inputs] --block-mode window
```
(Recipe source: `benchmarks/p80/scripts/07_run_cactus_em_p80.sh`, method `star2`.)

## Active — panel / matrix / sim prep (inputs to the estimator, not h itself)

| file | role |
|---|---|
| `build_kmer_pa.py` | builds `kmer_pa` (founder × k-mer membership matrix) |
| `build_var_pa.py` | builds `var_pa` / `var_called` (founder × variant carrier + called mask) |

`src/` now holds only the kMate program (estimator + kmer_pa builders). Tooling that
consumes the program's *outputs* or makes *simulation inputs* lives elsewhere:
- results aggregation / scoring → `benchmarks/scripts/` (`aggregate_results.py`, `aggregate_seedmix_validation.py`)
- pool simulation → `sims/scripts/` (`build_g0_uniform_sim.py`)

## Archived (moved here, not deleted) — `archive/`

Superseded or experimental code, kept for provenance. Not part of the
production estimator.

| archived file | why |
|---|---|
| `pre_cleanup_2026-05-26/{em_solver,block_em,per_sample_per_chrom}.py` | full pre-cleanup snapshots of the three edited files |
| `ld_blocks.py` | LD-block partitioning — only used by the removed `ld_gabriel`/`ld_complete` modes |
| `per_sample_driver.py` | older genome-wide (non-per-chrom) driver, superseded by `per_sample_per_chrom.py`; archived together with its legacy validation scripts (now in `tests/archive/`). The current harness uses the `*_per_chrom*` tests. Its window path was also broken (stale `project_blocks_to_records` signature). |
| `batch_runner.py`, `calibrate_alt_freqs.py` | helpers that depended on `per_sample_driver.py` |
| `per_sample_bigld_haplotype.py` | BigLD per-haplotype estimator (hapFIRE comparison) |
| `block_haplotype_bigld.py` | BigLD per-block-haplotype EM driver, moved out of `block_haplotype_em.py` (which now keeps only `smooth_h_across_blocks` for window mode) |
| `build_kmer_pa_from_fastas.py` | alt kmer_pa builder from founder FASTAs; no production caller (production kmer_pa uses `build_kmer_pa.py` from the PanGenie index) |
| `em_fixes/` | experimental H-estimation variant sweep (combined/balanced/balancedbubble/perbubble/whitening/invac/correlation/normalizer_factorial + score_*/run_* drivers); tried, documented in `docs/METHODS_TRIED_AND_RESULTS.md` |
| `build_subsampled_cn.py` | kmer_pa k-mer-balancing transform (per-stratum subsample / cactus-vs-PG class-match). **Tested and rejected** — all subsample variants lost to plain `filt2` on AF MAE (`docs/METHODS_TRIED_AND_RESULTS.md` §1–2). Production filter is **filt2 only** (drop ac=1). Kept for provenance. |
| `validate_seedmix_recipe.py` | recipe-validation helper |
| `block_solver.py`, `joint_solver.py`, `hapfire_solver.py` | earlier solver prototypes |
| `sweep_shape_norm_h_only.py` | a one-off k-mer-rebalancing sweep |

### Code removed from the active files (preserved in the snapshots)

- `em_solver.py`: `solve_em_with_omega` (contamination-ω variant; never used).
- `block_em.py`: `solve_em_clustered_per_block` + `find_haplotype_clusters` (earlier per-block EM, superseded by `solve_em_per_block`); `assign_kmers_to_blocks_multi` + `project_blocks_to_records_overlap` (overlapping-windows "Route 2").
- `per_sample_per_chrom.py`: the `ld_gabriel`/`ld_complete`/`bigld_panel` block modes; overlapping windows (`--window-step`); the k-mer-budget balancing flags (`--row-normalize-kmer_pa`, `--kf-correction-alpha`) and `_apply_kf_correction`; and `--ac-weight-counts`. All were off-by-default experiments (see `ALGORITHM.md` §8).
