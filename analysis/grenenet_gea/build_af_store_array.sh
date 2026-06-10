#!/bin/bash
#SBATCH --job-name=af_store
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=2:00:00
#SBATCH --requeue
#SBATCH --output=logs/af_store_%A_%a.out
#SBATCH --error=logs/af_store_%A_%a.err

# Collapse the 2,168 per-sample AF TSVs into the compact per-sample NPY store
# (uint16 4-decimal AF + uint8 n_called, split SNP / non-SNP; shared index once).
# One array task converts CHUNK consecutive samples; per-sample vectors are
# skipped if already present, so a requeued/preempted task resumes cleanly.
#
# Usage (run `init` ONCE first, then submit the array):
#   PY=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python
#   $PY analysis/grenenet_gea/build_af_store.py init --out-dir results/grenenet_gea/af_store
#   N=$($PY -c "import json;print(json.load(open('results/grenenet_gea/af_store/meta.json'))['n_samples'])")
#   CHUNK=70; NT=$(( (N + CHUNK - 1) / CHUNK ))
#   sbatch --array=1-${NT}%32 --export=ALL,OUT_DIR=results/grenenet_gea/af_store,CHUNK=$CHUNK \
#     analysis/grenenet_gea/build_af_store_array.sh

set -uo pipefail
cd /global/scratch/users/tbellg/kmate
mkdir -p logs
PY=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python

: ${OUT_DIR:?Set OUT_DIR to the af_store directory (after running init)}
CHUNK=${CHUNK:-70}
START=$(( (SLURM_ARRAY_TASK_ID - 1) * CHUNK ))

echo "[$(date)] task=${SLURM_ARRAY_TASK_ID} start=${START} count=${CHUNK} out=${OUT_DIR}"
$PY -u analysis/grenenet_gea/build_af_store.py convert \
    --out-dir "$OUT_DIR" --start "$START" --count "$CHUNK"
echo "[$(date)] task=${SLURM_ARRAY_TASK_ID} done"
