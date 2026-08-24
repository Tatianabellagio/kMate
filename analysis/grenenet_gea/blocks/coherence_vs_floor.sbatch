#!/bin/bash
#SBATCH --job-name=coh_floor
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=4:00:00
#SBATCH --output=logs/coh_floor_%j.out
#SBATCH --error=logs/coh_floor_%j.out
# Coherence cost of coarsening: per-window PC1-VE on gen9 pools, base + floor 8/15/25/40.
set -euo pipefail
cd /global/scratch/users/tbellg/kmate
mkdir -p logs
PY=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python
echo "[$(date)] coherence vs floor"
$PY -u analysis/grenenet_gea/blocks/coherence_vs_floor.py
echo "[$(date)] done"
