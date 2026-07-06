#!/bin/bash
#SBATCH --job-name=rdavp
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=2:00:00
#SBATCH --output=/global/scratch/users/tbellg/kmate/analysis/grenenet_gea/gea_newpanel/rdavp_%j.out
#SBATCH --error=/global/scratch/users/tbellg/kmate/analysis/grenenet_gea/gea_newpanel/rdavp_%j.out
set -euo pipefail
cd /global/scratch/users/tbellg/kmate
/global/home/users/tbellg/miniforge3/envs/kmate/bin/python -u analysis/grenenet_gea/gea_newpanel/rdavp_varpart.py
echo DONE_RDAVP
