#!/bin/bash
#SBATCH --job-name=predict_ecoperf
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=4:00:00
#SBATCH --output=logs/predict_ecoperf_%j.out
#SBATCH --error=logs/predict_ecoperf_%j.out
# Build cached per-founder-per-site selection matrix S (heavy Lustre I/O, once) then fit the
# kinship-vs-climate-match variance-partition + leave-one-site-out CV model.
set -euo pipefail
cd /global/scratch/users/tbellg/kmate
mkdir -p logs
PY=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python
echo "[$(date)] host=$(hostname)"
$PY analysis/grenenet_gea/predict_ecotype_performance.py
echo "[$(date)] done"
