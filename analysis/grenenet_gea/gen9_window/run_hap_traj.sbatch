#!/bin/bash
#SBATCH --job-name=hap_traj
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=0:20:00
#SBATCH --output=logs/hap_traj_%j.out
#SBATCH --error=logs/hap_traj_%j.err
# Pipeline B v2 Step 1+2: per-site haplotype trajectories p̄_{t,g} + point variance V_{t,g}.
set -euo pipefail
hostname
cd /global/scratch/users/tbellg/kmate
source ~/miniforge3/etc/profile.d/conda.sh && conda activate kmate
python analysis/grenenet_gea/build_hap_trajectories.py
