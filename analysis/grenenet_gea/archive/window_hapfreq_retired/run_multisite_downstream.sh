#!/bin/bash
# Downstream of the multisite founder GWAS, parametrized by membership tag + output suffix:
#   clade-level winners (kmate) -> climate permutation null (kmate) -> notebook (basic).
# Usage: run_multisite_downstream.sh <MEMB_TAG> <OUT_SUFFIX>
#   e.g.  run_multisite_downstream.sh clq90 _clq90
#         run_multisite_downstream.sh clq50 _clq50
#         run_multisite_downstream.sh K500 ""        # rebuild the original
set -euo pipefail
cd /global/scratch/users/tbellg/kmate
export MEMB_TAG="${1:?membership tag}" OUT_SUFFIX="${2-}"
KPY=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python
BPY=/global/home/users/tbellg/miniforge3/envs/basic/bin/python
echo "[$(date)] downstream MEMB_TAG=$MEMB_TAG OUT_SUFFIX='$OUT_SUFFIX'"
echo "[$(date)] 1/3 clade-level winners (kmate)"
$KPY analysis/grenenet_gea/archive/window_hapfreq_retired/cross_site_winners_multisite.py
echo "[$(date)] 2/3 climate permutation null (kmate)"
$KPY analysis/grenenet_gea/archive/window_hapfreq_retired/multisite_climate_perm.py
echo "[$(date)] 3/3 notebook (basic)"
$BPY analysis/grenenet_gea/archive/window_hapfreq_retired/_build_multisite_gwas_nb.py
echo "[$(date)] downstream done -> notebooks/multisite_founder_gwas${OUT_SUFFIX}.ipynb"
