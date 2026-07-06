#!/bin/bash
#SBATCH --job-name=hap_extract
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --cpus-per-task=2
#SBATCH --mem=6G
#SBATCH --time=0:15:00
#SBATCH --requeue
#SBATCH --output=logs/hap_extract_%A_%a.out
#SBATCH --error=logs/hap_extract_%A_%a.err
set -euo pipefail
cd /global/scratch/users/tbellg/kmate
PY=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python
LIST=results/grenenet_gea/window_vs_global/samples.txt
OFFSET=${OFFSET:-0}
ROW=$(( OFFSET + SLURM_ARRAY_TASK_ID ))
SAMPLE=$(sed -n "${ROW}p" "$LIST")
[ -z "$SAMPLE" ] && { echo "no sample at $ROW"; exit 0; }
$PY analysis/grenenet_gea/gen9_window/extract_hap_sample.py "$SAMPLE"
