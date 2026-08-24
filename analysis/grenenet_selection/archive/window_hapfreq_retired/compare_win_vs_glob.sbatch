#!/bin/bash
#SBATCH --job-name=af_wvg
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --cpus-per-task=2
#SBATCH --mem=6G
#SBATCH --time=0:30:00
#SBATCH --requeue
#SBATCH --output=logs/af_wvg_%A_%a.out
#SBATCH --error=logs/af_wvg_%A_%a.err
set -euo pipefail
cd /global/scratch/users/tbellg/kmate

PY=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python
GDIR=results/grenenet_kmate_arch3
WDIR=results/grenenet_kmate_window
OUTDIR=results/grenenet_gea/window_vs_global
LIST=${LIST:-$OUTDIR/samples.txt}
OFFSET=${OFFSET:-0}

ROW=$(( OFFSET + SLURM_ARRAY_TASK_ID ))
SAMPLE=$(sed -n "${ROW}p" "$LIST")
[ -z "$SAMPLE" ] && { echo "no sample at row $ROW"; exit 0; }

OUT=$OUTDIR/per_sample/${SAMPLE}.csv
mkdir -p $OUTDIR/per_sample
$PY analysis/grenenet_gea/compare_window_vs_global_af.py "$SAMPLE" "$GDIR" "$WDIR" "$OUT"
