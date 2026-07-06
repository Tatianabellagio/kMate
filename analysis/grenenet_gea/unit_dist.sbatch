#!/bin/bash
#SBATCH --job-name=unit_dist
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=4
#SBATCH --mem=120G
#SBATCH --time=4:00:00
#SBATCH --output=logs/unit_dist_%j.out
#SBATCH --error=logs/unit_dist_%j.out
set -eo pipefail
cd /global/scratch/users/tbellg/kmate
mkdir -p logs
/global/home/users/tbellg/miniforge3/envs/kmate/bin/python -u \
  analysis/grenenet_gea/unit_distributions.py 2>&1 | grep -vE "QUALITY CHECK|window_size|variants kept"
