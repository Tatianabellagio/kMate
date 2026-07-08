#!/bin/bash
#SBATCH --job-name=g9win_extract
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --cpus-per-task=2
#SBATCH --mem=6G
#SBATCH --time=0:20:00
#SBATCH --requeue
#SBATCH --output=logs/g9win_%A_%a.out
#SBATCH --error=logs/g9win_%A_%a.err
set -euo pipefail
cd /global/scratch/users/tbellg/kmate
PY=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python
LIST=results/grenenet_gea/gen9_window/gen9_members.txt
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$LIST")
[ -z "$SAMPLE" ] && { echo "no sample at ${SLURM_ARRAY_TASK_ID}"; exit 0; }
$PY analysis/grenenet_gea/gen9_window/extract_sample.py "$SAMPLE"
