#!/bin/bash
#SBATCH --job-name=ksweep_hits_re
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=2:00:00
#SBATCH --output=/global/scratch/users/tbellg/kmate/analysis/grenenet_gea/gea_newpanel/_regen2_logs/ksweep_hits_%j.log
#SBATCH --error=/global/scratch/users/tbellg/kmate/analysis/grenenet_gea/gea_newpanel/_regen2_logs/ksweep_hits_%j.log
# K* selection + RAW-p hit calling on the fresh ksweep_calib lambda tables + raw_p dumps.
set -euo pipefail
cd /global/scratch/users/tbellg/kmate
PY=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python
$PY -u analysis/grenenet_gea/gea_newpanel/ksweep_hits.py
echo DONE_KSWEEP_HITS
