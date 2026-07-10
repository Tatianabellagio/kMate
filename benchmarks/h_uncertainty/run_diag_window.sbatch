#!/bin/bash
#SBATCH --job-name=diag_emwin
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=2:00:00
#SBATCH --output=logs/diag_emwin_%j.out
#SBATCH --error=logs/diag_emwin_%j.out
# Window-mode: per-window ĥ error vs coverage + local k-mer support.
set -euo pipefail
cd /global/scratch/users/tbellg/kmate
mkdir -p logs benchmarks/h_uncertainty/results
PY=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python
$PY benchmarks/h_uncertainty/diag_window.py \
    --pool cov10_n50_g0_s42_hotspots_p80_chr1 \
    --window-bp 10000 --covs 3,10,30 --reps 3 --min-kmers 200
echo "=== done $(date) ==="
