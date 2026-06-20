#!/bin/bash
#SBATCH --job-name=fix_oldpanel
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=0:30:00
#SBATCH --requeue
#SBATCH --output=logs/fix_oldpanel_%A_%a.out
#SBATCH --error=logs/fix_oldpanel_%A_%a.err
set -euo pipefail
cd /global/scratch/users/tbellg/kmate
PY=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python
DIR=${DIR:-results/grenenet_kmate_arch3}
MASK=results/grenenet_gea/af_store/old2new_mask.npy
ARCH=${ARCH:-results/grenenet_kmate_arch3_oldpanel_archive}
LIST=${LIST:-results/grenenet_gea/window_vs_global/oldpanel_175.txt}

SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$LIST")
[ -z "$SAMPLE" ] && { echo "no sample at row $SLURM_ARRAY_TASK_ID"; exit 0; }
$PY analysis/grenenet_gea/correct_oldpanel_sample.py "$SAMPLE" "$DIR" "$MASK" "$ARCH"
