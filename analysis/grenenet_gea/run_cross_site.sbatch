#!/bin/bash
#SBATCH --job-name=cross_site_win
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=2:30:00
#SBATCH --output=analysis/grenenet_gea/logs/cross_site_%j.out
#SBATCH --error=analysis/grenenet_gea/logs/cross_site_%j.err
# Cross-site winners / clade x climate structure. Heavy Lustre I/O (per-founder global_h across
# all sites' evolved samples) -> sbatch with time. See cross_site_winners.py.
set -euo pipefail
hostname
cd /global/scratch/users/tbellg/kmate
source ~/miniforge3/etc/profile.d/conda.sh && conda activate kmate
python analysis/grenenet_gea/cross_site_winners.py
echo DONE
