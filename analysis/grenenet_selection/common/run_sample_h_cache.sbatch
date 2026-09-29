#!/bin/bash
#SBATCH --job-name=samp_h_cache
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=3:00:00
#SBATCH --output=analysis/grenenet_selection/logs/samp_h_cache_%j.out
#SBATCH --error=analysis/grenenet_selection/logs/samp_h_cache_%j.err
# One-time per-sample founder-h cache for the Price/fitness analyses (Lustre-heavy).
set -euo pipefail
hostname
cd /global/scratch/users/tbellg/kmate
source ~/miniforge3/etc/profile.d/conda.sh && conda activate kmate
python analysis/grenenet_selection/common/build_sample_h_cache.py
