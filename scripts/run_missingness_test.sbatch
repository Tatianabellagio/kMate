#!/bin/bash
#SBATCH --job-name=miss_test
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=2
#SBATCH --mem=48G
#SBATCH --time=4:00:00
#SBATCH --output=logs/miss_test_%j.out
#SBATCH --error=logs/miss_test_%j.err
# Headless run of the missingness-causation test (notebook -> script).
# Uses the 'basic' env (has numpy/scipy/matplotlib/pysam; old 'hapfm' is gone).
mkdir -p logs
set -euo pipefail
BASE=/global/scratch/users/tbellg/kmate
PY=/global/home/users/tbellg/miniforge3/envs/basic/bin/python
CHR=${1:-Chr1}
echo "[$(date)] missingness test, $CHR, env=basic"
$PY -u $BASE/scripts/run_missingness_test.py "$CHR"
echo "[$(date)] DONE"
