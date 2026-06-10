# kMate SV-GEA export → MOI-LAB shared Drive

Pushed to `gdrive:PROJECTS/grenenet/GrENE-net_PHASE1SV/kmate_gea_export/` by
`push_gea_to_drive.sh`, so the GEA work survives the cluster outage (started 2026-06-08).

## What's in the export
| Drive subfolder | Source on cluster | Size |
|---|---|---|
| `code/` | `kmate/analysis/grenenet_gea/` | ~8 MB |
| `external/` | scattered absolute paths (below) | ~146 MB |
| `results_grenenet_gea/` | `kmate/results/grenenet_gea/` (af_store 52G + matrices + group_means + gea/lfmm) | ~67 GB |

## External inputs — original absolute paths (hardcoded in `lib.py`)
To run the notebooks off the Drive copy, the absolute paths in `lib.py` must be
re-pointed at `external/`:

| `lib.py` const | Original path | In export |
|---|---|---|
| `T5` | `/global/scratch/users/tbellg/pang/grenenet_reads/Table_S5_sample_collection_sequencing_library.csv` | `external/Table_S5_sample_collection_sequencing_library.csv` |
| `SAMPLES_DATA` | `…/grenenet-phase1/frequency/hapFIRE_frequencies/samples_data_fix57.csv` | `external/samples_data_fix57.csv` |
| `BIOCLIM` | `…/grenenet-phase1/drive_zenodo/data-intermediate/bioclimvars_experimental_sites_era5.csv` | `external/bioclimvars_experimental_sites_era5.csv` |
| `LD_BLOCKS` | `…/gea_grene-net/ARCHIVE/linages_wza_picmin/kendall_0_w_id_n_blocks.csv` | `external/kendall_0_w_id_n_blocks.csv` |
| `ARA_KEYS` | `/global/home/users/tbellg/ara_key_files` | `external/ara_key_files/` |
| `PROJ` / `GEA` | `/global/scratch/users/tbellg/kmate` → `results/grenenet_gea` | `results_grenenet_gea/` |

## NOT exported (stays on cluster; rebuildable when it's back)
- `results/grenenet_kmate_arch3/` — 1.6 TB raw per-sample AF TSVs (af_store is the
  distilled version of these; rebuild af_store with `build_af_store.py` if needed).
- The panel / benchmarks / archive (the other ~2.3 TB).
