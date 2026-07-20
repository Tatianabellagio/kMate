#!/bin/bash
#SBATCH --job-name=tagging_sweep
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH --time=00:45:00
#SBATCH --array=0-35
#SBATCH --output=logs/tagging_sweep_%A_%a.out
#SBATCH --error=logs/tagging_sweep_%A_%a.out

set -uo pipefail
cd /global/scratch/users/tbellg/kmate/analysis/grenenet_gea
PY=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python3
mkdir -p logs

LINE=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" tagging_sweep_jobs.txt)
read -r MAC MINPAIR OUTDIR CLS MODE CHROM <<< "$LINE"
mkdir -p "$OUTDIR"

echo "[$(date)] task $SLURM_ARRAY_TASK_ID: $OUTDIR $CLS $MODE $CHROM (mac=$MAC min_pair_frac=$MINPAIR)"
"$PY" build_tagging_masked.py --chrom "$CHROM" --anchor-class "$CLS" --panel-mode "$MODE" \
    --mac-floor "$MAC" --min-pair-frac "$MINPAIR" --out-dir "$OUTDIR"
echo "[$(date)] done"
