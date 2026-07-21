#!/bin/bash
#SBATCH --job-name=k_sweep_re
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=6:00:00
#SBATCH --output=/global/scratch/users/tbellg/kmate/analysis/grenenet_gea/gea_newpanel/_regen2_logs/k_sweep_%j.log
#SBATCH --error=/global/scratch/users/tbellg/kmate/analysis/grenenet_gea/gea_newpanel/_regen2_logs/k_sweep_%j.log
# LF-adjusted partial-rank K sweep (pure python, fresh class matrices) + figure (basic env).
set -euo pipefail
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8
cd /global/scratch/users/tbellg/kmate
PYK=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python
PYB=/global/home/users/tbellg/miniforge3/envs/basic/bin/python
G=analysis/grenenet_gea/gea_newpanel
echo "== run_k_sweep (kmate) =="
$PYK -u $G/run_k_sweep.py
echo "== plot_k_sweep (basic) =="
$PYB -u $G/plot_k_sweep.py || echo "WARN plot_k_sweep failed (figure only)"
echo DONE_K_SWEEP
