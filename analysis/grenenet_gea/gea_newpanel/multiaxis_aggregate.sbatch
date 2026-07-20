#!/bin/bash
#SBATCH --job-name=ma_aggregate
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G
#SBATCH --time=2:00:00
#SBATCH --output=/global/scratch/users/tbellg/kmate/analysis/grenenet_gea/gea_newpanel/_regen2_logs/ma_aggregate_%j.log
#SBATCH --error=/global/scratch/users/tbellg/kmate/analysis/grenenet_gea/gea_newpanel/_regen2_logs/ma_aggregate_%j.log
set -euo pipefail
cd /global/scratch/users/tbellg/kmate/analysis/grenenet_gea/gea_newpanel
BPY=/global/home/users/tbellg/miniforge3/envs/basic/bin/python
$BPY -u _build_gea_class_gain_multiaxis_nb.py
echo DONE_MA_AGGREGATE
