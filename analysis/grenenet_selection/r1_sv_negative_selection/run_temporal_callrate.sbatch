#!/bin/bash
#SBATCH --job-name=sv_temporal_cr
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G
#SBATCH --time=6:00:00
#SBATCH --output=logs/temporal_callrate_%j.out
#SBATCH --error=logs/temporal_callrate_%j.out
# Temporal (Fig. 2d) chain on Truvari events with a relaxed founder call-rate floor.
#   sbatch --export=ALL,SV_CALL_MIN=0.7 run_temporal_callrate.sbatch
set -eo pipefail
cd /global/home/users/tbellg/kmate/analysis/grenenet_selection/r1_sv_negative_selection
PY=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python
export SV_CALL_MIN=${SV_CALL_MIN:-0.7} SV_UNIT=event SV_MIN_MAC=2
CR=$($PY -c "import os; print(int(round(float(os.environ['SV_CALL_MIN'])*100)))")
echo "[$(date)] host=$(hostname) SV_CALL_MIN=$SV_CALL_MIN"
$PY -u _sv_temporal_climate_matched.py
$PY -u _sv_purging_classify.py
SV_TAG=_event2_cr${CR} $PY -u _sv_purging_annotate.py
echo "[$(date)] done"
