#!/bin/bash
#SBATCH --job-name=bench_hunc
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=6:00:00
#SBATCH --output=logs/bench_hunc_%j.out
#SBATCH --error=logs/bench_hunc_%j.out
# Benchmark analytic Fisher SE on ĥ vs parametric bootstrap (gold) on p80.
set -euo pipefail
cd /global/scratch/users/tbellg/kmate
mkdir -p logs benchmarks/h_uncertainty/results
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8
PY=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python
# Lighter than the 4h-timeout first attempt: B=50, faster EM refits (max_iter
# 120, tol 1e-5 in the script). 50×2×3=300 refits instead of 720.
$PY benchmarks/h_uncertainty/bench_h_uncertainty.py \
    --pool cov10_n50_g0_s42_hotspots_p80_chr1 \
    --covs 5,10,30 --B 50 \
    --rconds 1e-10,1e-4,1e-2,5e-2 --support-eps 1e-3
echo "=== done $(date) ==="
