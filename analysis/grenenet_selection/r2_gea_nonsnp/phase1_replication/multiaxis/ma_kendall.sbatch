#!/bin/bash
#SBATCH --job-name=ma_kendall
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=2:00:00
#SBATCH --array=0-59
#SBATCH --output=/global/scratch/users/tbellg/kmate/logs/ma_kendall_%A_%a.out
#SBATCH --error=/global/scratch/users/tbellg/kmate/logs/ma_kendall_%A_%a.out
# Multi-axis Kendall-tau: per-record tau(AF, axis) across 355 pools, gen9.
# One array task per (axis, class). Output -> multiaxis/kendall/kendall_{cls}_gen9_{axis}.csv
set -euo pipefail
cd /global/scratch/users/tbellg/kmate
source analysis/grenenet_selection/r2_gea_nonsnp/phase1_replication/multiaxis/axes.sh
PY=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python
i=${SLURM_ARRAY_TASK_ID}; AXIS=$(axis_for $i); C=$(cls_for $i)
OUT=analysis/grenenet_selection/r2_gea_nonsnp/phase1_replication/results/multiaxis/kendall
echo "== host $(hostname) | kendall $C gen9 $AXIS =="
$PY -u analysis/grenenet_selection/r2_gea_nonsnp/phase1_replication/run_kendall.py \
    --class "$C" --gen 9 --climate "$AXIS" --threads "${SLURM_CPUS_PER_TASK:-8}" --out "$OUT"
echo "== ALLDONE_MA_KENDALL $C $AXIS =="
