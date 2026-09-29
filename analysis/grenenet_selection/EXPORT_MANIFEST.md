# kMate phase-1 GEA results → MOI-LAB shared Drive

Pushed to `gdrive:PROJECTS/grenenet/GrENE-net_PHASE1SV/data/grenenet_gea/phase1_replication/`
by `push_gea_to_drive.sh`. **Results-only**, so the advisor can plot the phase-1
replication himself for talks. The Drive layout mirrors the cluster folder names
(`grenenet_gea/phase1_replication/wza`), rooted under `data/`.

> **⚠️ Stale exports (2026-07-08).** The WZA CSVs below were built on **PRE-FIX AF**
> (`results/grenenet_kmate_arch3/`, now retired) and **predate the full-panel Kf_w / `--unit chrom`
> regeneration**. They should be **re-exported** after the downstream regen (production AF is now
> `analysis/grenenet_selection/common/rerun_kfw_hb/{seedmix,evolved}/`) before being cited or presented.

## What's in the export (~7 MB)
| Drive path (under `data/grenenet_gea/phase1_replication/`) | Source on cluster | What it is |
|---|---|---|
| `wza/wza_*_gen9_bio1_deg7cap2000.csv` | same path under `results/` | Final per-block WZA p-values: last-gen (gen9), bio1, primary deg7-cap2000 correction, 3 models (kendall / lfmm / binomial) × 3 classes (snp / smallindel / sv) = 9 CSVs |
| `manhattan_3models_deg7cap2000.png` | same path under `results/` | Reference 3-model Manhattan figure |

Each WZA CSV is one row per LD block; the per-block significance is the `Z_pVal`
column (Manhattan y = −log10 `Z_pVal`).

## Deliberately NOT exported
- **Code** — lives on GitHub (`Tatianabellagio/kMate`); no need to duplicate on Drive.
- **Raw allele frequencies** — `analysis/grenenet_selection/common/results/af_store/` (52 GB), `gen_matrices/`
  (43 GB), `pool_matrices/` (35 GB), `class_matrices/`, LFMM `*.f64` memmaps. All large
  and rebuildable from the pipeline (`build_af_store.py` → `build_gen_matrices.py` →
  `build_pool_matrix.py`); they are not plotting inputs for the advisor.
- **Raw per-sample AF TSVs** — `results/grenenet_kmate_arch3/` (1.6 TB), the panel,
  benchmarks, archive.

## History
A previous push (2026-06-08, `kmate_gea_export/`) parked the full ~67 GB working set
(code + external inputs + all of `results/grenenet_gea/`) on Drive ahead of a cluster
outage. That bulk copy is superseded by this lightweight results-only export.
