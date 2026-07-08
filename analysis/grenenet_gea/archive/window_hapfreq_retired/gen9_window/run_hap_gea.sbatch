#!/bin/bash
#SBATCH --job-name=hap_gea
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=0:40:00
#SBATCH --output=logs/hap_gea_%j.out
#SBATCH --error=logs/hap_gea_%j.err
# Pipeline B v2 Step 3+4: per-site s_g, climate β1 + site-permutation null, block FDR.
set -euo pipefail
hostname
cd /global/scratch/users/tbellg/kmate
source ~/miniforge3/etc/profile.d/conda.sh && conda activate kmate
export N_PERM=${N_PERM:-5000}
echo "N_PERM=$N_PERM"
python analysis/grenenet_gea/build_hap_gea.py
