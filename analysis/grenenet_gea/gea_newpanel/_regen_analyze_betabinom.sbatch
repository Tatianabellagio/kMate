#!/bin/bash
#SBATCH --job-name=an_beta
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=3:00:00
#SBATCH --output=/global/scratch/users/tbellg/kmate/results/grenenet_gea/gea_newpanel/_regen_logs/analyze_betabinom_%j.log
set -euo pipefail
cd /global/scratch/users/tbellg/kmate
PY=/global/home/users/tbellg/miniforge3/envs/basic/bin/python
$PY -u analysis/grenenet_gea/gea_newpanel/analyze_betabinom.py
echo DONE_ANALYZE_BETABINOM
