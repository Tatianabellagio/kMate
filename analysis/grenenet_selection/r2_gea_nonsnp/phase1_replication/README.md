# GEA: genome-environment association of allele-frequency change, SNPs, short indels and SV events

Every GEA number and figure in the paper comes from the chain below. `bash _submit_gea.sh` submits it to the cluster
as one chain of Slurm jobs, in this order. (Rewritten 2026-10-20 at the GEA audit; the earlier README and `STATUS.md`
described the retired WZA / bio1-only phase-1 replication and are in `archive_gea/`.)

**Panel.** Production is **arch3v2** (since 2026-10-22, `panel/arch3v2/README.md`). `_submit_gea.sh` reads
`KMATE_PANEL` (default arch3v2), exports it to every job, and every output folder below carries the panel suffix
`_arch3v2` (`lib.VSUF`); logs go to `logs/gea_arch3v2/` at the project root, job ids to `results/_submit_gea_arch3v2_jobids.txt`. The
unsuffixed folders (`class_matrices_cr50/`, `lfmm_cr50/`, `tile_cr50/`, `event_gea/cr50/`, `genes/convergence/results/`)
are the previous panel (arch3, 352 pools in 31 gardens including site 33), kept only to reproduce old outputs
(`KMATE_PANEL=arch3`).

**Design.** Response = change in allele frequency of each variant from the founders to generation 9 (delta p, 346
pools in 30 gardens; site 33 is excluded, `lib.EXCLUDED_SITES`), regressed on each of 22 climate axes (BIOCLIM bio1-19 and climate PCs 1-3 of the gardens,
z-scored over pools) with LFMM ridge, K = 16 latent factors. **p-values are raw** (`lfmm_test` `$pvalue`, no
genomic-control or WZA correction: project decision); lambda is computed and stated, never applied. Significance is
**Bonferroni 0.05 / n tested units, per class x axis** (never divided across axes). Units: a SNP or short indel is one
record; an SV (|length change| > 50 bp) is one **Truvari event** (all panel alleles of a bubble collapsed). Filters,
both arms: pool-mean minor-allele frequency >= 0.05 and called in >= 116 of 231 founders (50%; an event is called at a
founder if any member record is). LD blocks: the clq r2 = 0.9 tiling (BigLD on dense common variants), assigned to
every class with `blocks/blocks_tiling.assign_tiling`, **not merged** (so a block id names the same interval for every
class).

| step | script (driver) | does | writes |
|---|---|---|---|
| 1 | `build_class_matrices.py` (`run_class_matrices_cr50.sbatch`) | record arm (snp, smallindel): the pools x records frequency matrix after the MAF and call-rate filters, and the pools' climate (bio1-19, pc1-3) | `results/class_matrices_cr50_arch3v2/` |
| 2 | `build_lfmm_input.py` + `multiaxis/run_lfmm_nogif.R` (`multiaxis/ma_lfmm_cr50.sbatch`, 44 tasks) | **record arm** LFMM: snp and smallindel x 22 axes | `results/multiaxis/lfmm_cr50_arch3v2/` |
| 3 | `genes/convergence/event_gea/build_event_matrix.py` (`run_event_matrix_cr50.sbatch`) | **event arm** inputs: delta p per Truvari event (sum of member frequencies), from `r1_sv_negative_selection/results/sv_adaptive_arch3v2/nonsnp_event_units.npz` | `genes/convergence/event_gea/cr50_arch3v2/{nonsnp,sv}_gen9_{units.csv,dp.npy}` |
| 4 | `build_lfmm_input.py --dp` + `multiaxis/run_lfmm_nogif.R` (`genes/convergence/event_gea/run_event_lfmm_cr50.sbatch`, 44 tasks) | event arm LFMM: the SV-event scan and the pooled non-SNP scan (short indels + SV events, one Bonferroni) | `event_gea/cr50_arch3v2/lfmm_event_{sv,nonsnp}_gen9_<axis>.csv.gz` |
| 5 | `multiaxis/reblock_cr50.py` (`multiaxis/run_reblock_cr50.sbatch`) | joins both arms into one schema and assigns the tiling block; short indels get their integer `unit` id | **`results/multiaxis/tile_cr50_arch3v2/lfmm_{snp,smallindel,sv,nonsnp}_gen9_<axis>.csv`** (the production GEA table) |
| 6 | `genes/convergence/class_signal_counts.py` | significant blocks per class: a block is significant for a class when its lead unit of that class clears that class's Bonferroni on >= 1 axis | `genes/convergence/results_arch3v2/class_signal_{counts,cells,blocks}.csv` |
| 7 | `genes/convergence/annotate_gea_universe.py`, `gea_hit_landscape_data.py`, `gea_sv_cargo_data.py` | every tested non-SNP unit: region (project location rule), hit axes (any of the sv / smallindel / pooled scans), SV content class | `genes/convergence/results_arch3v2/gea_{universe_regions,hit_landscape}.csv.gz`, `gea_sv_unit_cargo.csv` |
| 8 | `genes/convergence/gea_results_numbers.py` | **every number of the GEA Results section**, from the tables of steps 6-7 | `genes/convergence/results_arch3v2/gea_results_numbers.csv` |
| 9 | `multiaxis/raw_block_significance.py` (diagnostic) | BH and GIF counts, genome fraction, lambda per class x axis | `results/multiaxis/raw_block_significance_lfmm_tile_cr50_arch3v2{,_union}.csv` |

Figures: `notebooks/gea_draft_arch3v2.ipynb` (main GEA figure, `paper_figures/panel_figure_gea/`; kernel `basic`) and
`notebooks/lfmm/gea_hit_landscape_arch3v2.ipynb` (built by `KMATE_PANEL=arch3v2 python
genes/convergence/_build_gea_hit_landscape_nb.py`). The candidate list
(`genes/convergence/results_arch3v2/candidate_list_v2_events_cr50/`) is built from `tile_cr50_arch3v2` by
`genes/convergence/v2/run_v2.sh`. `notebooks/gea_draft.ipynb` and `notebooks/lfmm/gea_hit_landscape.ipynb` are the
previous-panel (arch3) versions, kept as they are.

Inputs not made here: the kMate pool frequency store (`common/results/af_store_arch3v2`,
`common/results/pool_matrices_arch3v2`; its generation-9 matrices, one pool per plot at the plot's last generation, are
made by `assemble_lastgen.py --kind snp` and `--kind nonsnp` from phase 1's `final_gen.csv`), the arch3v2 founder panel (`panel/arch3v2/`), the event units (`r1_sv_negative_selection/_build_event_units.py`), the tiling
block map (`blocks/`), garden climate (`common/build_garden_climate.py`).

`build_event_matrix.py` asserts that its short-indel delta p equals the step-1 matrices of the same panel
(`class_matrices_cr50_arch3v2/smallindel_gen9.*`), and reads its pools from `common/results/pool_matrices_arch3v2`
(the same 346 pools, in the same order, as `class_matrices_cr50_arch3v2/gen9.pools.csv`).

Environments: `kmate` (Python), `lfmm_env` (R, `lfmm`). Compute nodes only; steps 2 and 4 need ~180 GB / 48 GB per task.

## The model

Per class x axis, `run_lfmm_nogif.R` fits `lfmm_ridge(Y, X, K = 16)` and tests with `lfmm_test(..., calibrate = "gif")`.
Y = delta p (pools x tested units, from `build_lfmm_input.py`); X = one climate axis, z-scored over pools (`ddof = 0`).
`calibrate = "gif"` is passed only so `lfmm_test` returns lambda (`$gif`): the p-value used is `$pvalue`, which
`lfmm_test` sets before its calibration branch and never overwrites (`$calibrated.pvalue` is a separate field, unused).
Lambda recomputed from the written p-values matches `$gif`, which would be ~1 had the p-values been calibrated. K = 16
is chosen in `notebooks/10_lfmm_k_selection.ipynb`. The climate PCs are the PCs of the 19 standardized BIOCLIM
variables over the pools, z-scored (`build_class_matrices.add_climate_pcs`); pc1 is oriented to +bio1, pc2 and pc3 to
their largest loading (bio7, bio10). The kMate pool frequencies have no missing values (`build_lfmm_input.py` asserts
it), so nothing is imputed.

## Checks

arch3v2 (`tile_cr50_arch3v2`, counted from the tables 2026-10-23):

- Each class has the same number of rows on all 22 axes; tested units per axis: SNP 1,923,468, short indel 649,971,
  SV events 24,308, pooled non-SNP 674,279 (as in the executed `gea_draft_arch3v2.ipynb`). On bio1, every table has
  MAF >= 0.05 and no NaN p, and the SV and non-SNP tables give every SV event the same block. 30 gardens, 346 pools; the
  record arm (`class_matrices_cr50_arch3v2/gen9.pools.csv`) and the event arm (`pool_matrices_arch3v2`) read the same
  pools in the same order.

Previous panel (arch3, 2026-10-20 audit; not repeated on arch3v2):

- All 88 `tile_cr50` tables: MAF > 0.05, no NaN p, the same rows on all 22 axes; tested units per axis: SNP 1,790,024,
  short indel 603,630, SV events 20,134, pooled non-SNP 623,764; the SV and non-SNP tables give every SV event the same
  block. 31 gardens, 352 pools; the record and event arms read identical pool tables.
- `run_lfmm_nogif.R` is the one LFMM script of both arms: on the production SV inputs (sv x bio1, sv x pc3) it gives the
  same p, B and z as the archived copy the event arm used to call, and reproduces the production tables to <= 2.5e-12.
- An independent recomputation from `tile_cr50` (no project code) reproduces every block and SV count in
  `gea_results_numbers.csv`.

## Caveats that travel with these numbers

- Raw p is uncalibrated (arch3v2 lambda per class x axis 0.78-3.21, median 1.78, highest sv x pc3; from
  `raw_block_significance_lfmm_tile_cr50_arch3v2.csv`); the hit list is a candidate net.
- Two hit rules exist: block counts use the class's own scan; the candidate pool also takes the pooled non-SNP scan,
  which adds one SV event (`gea_results_numbers.csv` reports both).
- The garden climate PCs used by the LFMM (`class_matrices_cr50_arch3v2/gen9.pools.csv`) differ slightly from
  `lib.garden_climate()` (`common/results/garden_climate.csv`, 30 gardens): pc1 by up to 0.036, pc2 by up to 0.071,
  pc3 by up to 0.050; bio1-19 identical.

## Other scripts in this folder

- `annotate_genes_tair_uniprot.py` (gene function; imports `annotate_gene_function.py` for its `mygene_batch` helper).
- `build_class_matrices_seedmix.py` (`build_class_matrices_seedmix.sbatch`, generation 0 = the 8 SEEDMIX reps) and
  `rebuild_class_matrices_gen123.sbatch` (generations 1-3, via `build_class_matrices.py`): build
  `results/class_matrices_cr50_gen0123_arch3v2/` (MAF >= 0.05 and call rate >= 50%, the GEA floor), the gen0-3 tables
  of the shared-Drive allele-frequency export (`push_af_matrices_to_drive.sh`). Not a GEA input.
- `build_lfmm_input_hap.py`: phase-1 replication experiment (haploblock LFMM input); not used by the paper. The LFMM
  drivers `run_lfmm_hap.sbatch` and `run_lfmm_lastgen.R` were archived 2026-10-23 (`archive_unused_2026-10-23/MANIFEST.csv`).
- `multiaxis/callrate_calibration.py`: evidence behind the call-rate floor (previous panel). `multiaxis/axis_inflation.py`,
  `multiaxis/callrate_floor_effect.py` and the exploratory `multiaxis/_build_*_nb.py` builders were archived 2026-10-23
  (`archive_unused_2026-10-23/MANIFEST.csv`).
