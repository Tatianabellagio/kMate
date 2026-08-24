#!/bin/bash
#SBATCH --job-name=an_quasi
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=2:00:00
#SBATCH --output=/global/scratch/users/tbellg/kmate/analysis/grenenet_gea/gea_newpanel/results/_regen_logs/analyze_quasibinom_%j.log
set -euo pipefail
cd /global/scratch/users/tbellg/kmate
PY=/global/home/users/tbellg/miniforge3/envs/basic/bin/python
$PY -u analysis/grenenet_gea/gea_newpanel/analyze_quasibinom.py
echo DONE_ANALYZE_QUASIBINOM
