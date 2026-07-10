#!/bin/bash
#SBATCH --job-name=af_calib
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=160G
#SBATCH --time=2:00:00
#SBATCH --output=logs/af_calib_%j.out
#SBATCH --error=logs/af_calib_%j.out
# Calibrate the AF identifiability floor c for a panel.  Usage: sbatch run_af_calibrate.sbatch {p80|231}
set -euo pipefail
cd /global/scratch/users/tbellg/kmate
mkdir -p logs benchmarks/h_uncertainty/results
PANEL="${1:?pass panel: p80 or 231}"
PY=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python
$PY benchmarks/h_uncertainty/af_calibrate_floor.py --panel "$PANEL"
echo "=== done $PANEL $(date) ==="
