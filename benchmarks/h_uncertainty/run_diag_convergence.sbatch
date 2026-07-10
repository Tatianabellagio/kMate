#!/bin/bash
#SBATCH --job-name=diag_emconv
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=3:00:00
#SBATCH --output=logs/diag_emconv_%j.out
#SBATCH --error=logs/diag_emconv_%j.out
# Decompose EM error: coverage(variance) vs identifiability(bias); f32 vs f64.
set -euo pipefail
cd /global/scratch/users/tbellg/kmate
mkdir -p logs benchmarks/h_uncertainty/results
PY=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python
$PY benchmarks/h_uncertainty/diag_convergence.py \
    --pool cov10_n50_g0_s42_hotspots_p80_chr1 \
    --covs 3,10,30,100 --f64-covs 10 --max-iter 2000
echo "=== done $(date) ==="
