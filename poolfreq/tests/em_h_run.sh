#!/bin/bash
#SBATCH --job-name=em_h
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=3:00:00
#SBATCH --output=logs/em_h_%j.out
#SBATCH --error=logs/em_h_%j.err
mkdir -p logs
set -euo pipefail
# Plain-EM h estimate (alpha=0, no shape norm) for one cn_full on one sim.
# Usage: sbatch em_h_run.sh <CN_PREFIX> <R1> <R2> <SAMPLE> <OUT_PREFIX>
CN_PREFIX=$1; R1=$2; R2=$3; SAMPLE=$4; OUT_PREFIX=$5
ROOT=/global/scratch/users/tbellg/hapfire_sv
PY=/global/home/users/tbellg/miniforge3/envs/hapfm/bin/python
# script lives in src/archive/ but imports em_solver/kmer_count from src/
export PYTHONPATH=$ROOT/poolfreq/src:${PYTHONPATH:-}
mkdir -p "$(dirname "$OUT_PREFIX")"
echo "[$(date)] EM h: cn=$CN_PREFIX sample=$SAMPLE"
$PY -u $ROOT/poolfreq/src/archive/sweep_shape_norm_h_only.py \
    --cn-prefix "$CN_PREFIX" \
    --reads "$R1" "$R2" \
    --sample "$SAMPLE" \
    --out-prefix "$OUT_PREFIX" \
    --alphas 0 --threads 8 \
    --counts-cache "${OUT_PREFIX}.counts.npy"
echo "[$(date)] DONE -> ${OUT_PREFIX}.h_sweep.npz"
