#!/bin/bash
#SBATCH --job-name=sm_p0
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=0:20:00
#SBATCH --output=logs/sm_p0_%j.out
#SBATCH --error=logs/sm_p0_%j.err
# Project the 8 SEEDMIX window reps -> founding haplotype-frequency p0 aligned to
# the cohort hapfreq matrix. Chain after run_seedmix_window.sbatch (afterok).
set -euo pipefail
hostname
cd /global/scratch/users/tbellg/kmate
source ~/miniforge3/etc/profile.d/conda.sh && conda activate kmate
python analysis/grenenet_gea/build_hapfreq_p0_seedmix.py
