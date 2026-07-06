#!/bin/bash
#SBATCH --job-name=hapfreq
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH --time=2:00:00
#SBATCH --array=0-15
#SBATCH --output=logs/hapfreq_%A_%a.out
#SBATCH --error=logs/hapfreq_%A_%a.out
# STEP 2: project per-block founder h -> per-sample haplotype-frequency matrix (chunked).
set -eo pipefail
cd /global/scratch/users/tbellg/kmate
mkdir -p logs results/grenenet_gea/hapfreq/chunks
PY=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python
N=2168; NCHUNK=16
SZ=$(( (N + NCHUNK - 1) / NCHUNK ))         # 136
START=$(( SLURM_ARRAY_TASK_ID * SZ ))
END=$(( START + SZ )); if [ $END -gt $N ]; then END=$N; fi
echo "[$(date)] chunk task $SLURM_ARRAY_TASK_ID : samples [$START,$END)"
$PY -u analysis/grenenet_gea/build_hapfreq_matrix.py \
    --start $START --end $END \
    --out-dir results/grenenet_gea/hapfreq/chunks 2>&1 | grep -vE "variants kept"
echo "[$(date)] done"
