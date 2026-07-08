# kMate source inventory

What's actually used to estimate founder frequencies (`h`) and project to
per-record allele frequencies. Last cleaned: 2026-05-26.

The code is now an installable package under [`kmate/`](kmate/) (`pip install -e .`
from the repo root), exposing a `kmate` command — `kmate run` is the estimator,
`kmate build-kmer-pa` / `build-var-pa` / `filter-pa` build the panel, and
`kmate selftest` verifies an install on a bundled fixture. The thin `src/*.py`
files are back-compat shims so existing `python src/<script>.py ...` callers
still work. Module roles below are unchanged by the packaging.

There is **one estimator**, selected by the unit it fits: `--unit {chrom,ld,bp,tsv}`.
Each unit is fit locally (haploblock collapse → EM → project). The default is
**`--unit chrom`** — one founder mixture per chromosome — which is the production
estimator for selfing / inbred / F0 pools (e.g. GrENE-Net). `--block-mode {global,window}`
remain as deprecated aliases (`global`→`--unit chrom`, `window`→`--unit bp`).

## Active — the estimator (h estimation + AF projection)

| file | role |
|---|---|
| `em_solver.py` | EM core. `solve_em()` — the founder-mixture estimator (E/M multiplicative update; supports the Dirichlet anchor used by the legacy window recipe). M-step normalization defaults to `normalize="per_founder"` (each founder divided by its own full-panel k-mer content `Kf_w`). |
| `kmer_count.py` | Jellyfish wrapper. Counts canonical k=31 k-mers in the reads → the count vector fed to the EM. |
| `ld_partition.py` | `CompleteLDPartition` — r²-LD blocks from the panel's own `var_pa` (the `--unit ld` partition; blocks are a panel property, computed once and cached). |
| `h_uncertainty.py` | Collapse-aware uncertainty on the `K_b` haplotype classes (Fisher info + bootstrap) for `--emit-af-se`. |
| `block_em.py` | Per-unit pieces (bp/tsv/ld units): `define_windows`, `assign_kmers_to_blocks`, `assign_records_to_blocks`, `solve_em_per_block` (per-unit EM), `project_blocks_to_records` (per-unit → per-record AF, missing-aware). |
| `block_haplotype_em.py` | `smooth_h_across_blocks` — Li–Stephens-style smoothing of per-window `h` (legacy `--no-local-only` window recipe only). |
| `per_sample_per_chrom.py` | Production driver. FASTQ/BAM → counts → EM (per unit) → AF projection → output TSV. Dispatches `--unit {chrom,ld,bp,tsv}` (default `chrom`); `--block-mode {global,window}` are deprecated aliases. |

### Production recipes

```bash
# --unit chrom (DEFAULT) — selfing / inbred / F0 pools (e.g. SEEDMIX), the
# GrENE-Net production estimator; one founder mixture per chromosome. Uses the
# defaults --normalize per_founder + --kmer-weight uniform (2026-07-06;
# see docs/FOUNDER_NORMALIZATION_FIX.md), so no weight/unit flags needed.
kmate run \
  --kmer-pa-prefix <kmer_pa>/kmer_pa \
  --var-pa <panel>.var_pa.npz --var-called <panel>.var_called.npz \
  --var-meta <panel>.meta.npz \
  --reads R1.fq R2.fq --sample <name> --out <name>.tsv \
  --threads 8 --chroms Chr1        # --unit chrom is the default

# --unit bp / legacy "star2" window — recombinant pools. Fixed-bp windows, fit
# local-only by default. The anchored+smoothed star2 recipe is now opt-in behind
# --no-local-only (--window-bp 10000 --global-anchor-weight 0.3
#  --hmm-smooth-passes 5 --hmm-smooth-alpha 0.5), and keeps the inv_mb weighting:
kmate run [same inputs] --unit bp --window-bp 10000 --no-local-only --kmer-weight inv_mb
```
(`--unit ld --ld-r2 0.1` fits r²-LD blocks; note it collapses in low-diversity
blocks (centromere) so it is not appropriate for selfing pools.)
(`python src/per_sample_per_chrom.py ...` still works via the shim.)
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
| `build_kmer_pa_from_fastas.py` | alt kmer_pa builder from founder FASTAs; no production caller (production kmer_pa uses `build_kmer_pa.py` from the in-house `ours_Chr{N}` index; the PanGenie index is a comparator only — see `scripts/build_kmer_pa_arch3.sh`) |
| `em_fixes/` | experimental H-estimation variant sweep (combined/balanced/balancedbubble/perbubble/whitening/invac/correlation/normalizer_factorial + score_*/run_* drivers); tried, documented in `docs/METHODS_TRIED_AND_RESULTS.md` |
| `build_subsampled_cn.py` | kmer_pa k-mer-balancing transform (per-stratum subsample / cactus-vs-PG class-match). **Tested and rejected** — all subsample variants lost to plain `filt2` on AF MAE (`docs/METHODS_TRIED_AND_RESULTS.md` §1–2). Production filter is **filt2inv** (drop ac=1 singletons AND ac=F invariants; `--filter-production`, see `docs/PIPELINE_STATE.md` §0). Kept for provenance. |
| `validate_seedmix_recipe.py` | recipe-validation helper |
| `block_solver.py`, `joint_solver.py`, `hapfire_solver.py` | earlier solver prototypes |
| `sweep_shape_norm_h_only.py` | a one-off k-mer-rebalancing sweep |

### Code removed from the active files (preserved in the snapshots)

- `em_solver.py`: `solve_em_with_omega` (contamination-ω variant; never used).
- `block_em.py`: `solve_em_clustered_per_block` + `find_haplotype_clusters` (earlier per-block EM, superseded by `solve_em_per_block`); `assign_kmers_to_blocks_multi` + `project_blocks_to_records_overlap` (overlapping-windows "Route 2").
- `per_sample_per_chrom.py`: the `ld_gabriel`/`ld_complete`/`bigld_panel` block modes; overlapping windows (`--window-step`); the k-mer-budget balancing flags (`--row-normalize-kmer_pa`, `--kf-correction-alpha`) and `_apply_kf_correction`; and `--ac-weight-counts`. All were off-by-default experiments (see `ALGORITHM.md` §8).
