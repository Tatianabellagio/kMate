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
  --cn-kmer-prefix <cn_full>/cn \
  --cn-var <panel>.cn_var.npz --cn-var-called <panel>.cn_var_called.npz \
  --cn-var-meta <panel>.meta.npz \
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
| `build_kmer_cn.py` | builds `cn` (founder × k-mer membership matrix) |
| `build_cn_var.py` | builds `cn_var` / `cn_var_called` (founder × variant carrier + called mask) |
| `build_kmer_cn_from_fastas.py` | alt cn builder from founder FASTAs |
| `build_g0_uniform_sim.py`, `build_subsampled_cn.py` | simulation / subsampling helpers |
| `aggregate_results.py`, `aggregate_seedmix_validation.py`, `validate_seedmix_recipe.py` | downstream aggregation / validation |

## Archived (moved here, not deleted) — `archive/`

Superseded or experimental code, kept for provenance. Not part of the
production estimator.

| archived file | why |
|---|---|
| `pre_cleanup_2026-05-26/{em_solver,block_em,per_sample_per_chrom}.py` | full pre-cleanup snapshots of the three edited files |
| `ld_blocks.py` | LD-block partitioning — only used by the removed `ld_gabriel`/`ld_complete` modes |
| `per_sample_driver.py` | the older genome-wide (non-per-chrom) driver; superseded by `per_sample_per_chrom.py` |
| `batch_runner.py`, `calibrate_alt_freqs.py` | helpers that depended on `per_sample_driver.py` |
| `per_sample_bigld_haplotype.py` | BigLD per-haplotype estimator (hapFIRE comparison) |
| `block_solver.py`, `joint_solver.py`, `hapfire_solver.py` | earlier solver prototypes |
| `sweep_shape_norm_h_only.py` | a one-off k-mer-rebalancing sweep |

### Code removed from the active files (preserved in the snapshots)

- `em_solver.py`: `solve_em_with_omega` (contamination-ω variant; never used).
- `block_em.py`: `solve_em_clustered_per_block` + `find_haplotype_clusters` (earlier per-block EM, superseded by `solve_em_per_block`); `assign_kmers_to_blocks_multi` + `project_blocks_to_records_overlap` (overlapping-windows "Route 2").
- `per_sample_per_chrom.py`: the `ld_gabriel`/`ld_complete`/`bigld_panel` block modes; overlapping windows (`--window-step`); the k-mer-budget balancing flags (`--row-normalize-cn`, `--kf-correction-alpha`) and `_apply_kf_correction`; and `--ac-weight-counts`. All were off-by-default experiments (see `ALGORITHM.md` §8).
