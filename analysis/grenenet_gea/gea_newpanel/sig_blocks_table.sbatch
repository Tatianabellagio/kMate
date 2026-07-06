#!/bin/bash
#SBATCH --job-name=sig_blocks
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G
#SBATCH --time=1:00:00
#SBATCH --output=/global/scratch/users/tbellg/kmate/analysis/grenenet_gea/gea_newpanel/sig_blocks_%j.out
#SBATCH --error=/global/scratch/users/tbellg/kmate/analysis/grenenet_gea/gea_newpanel/sig_blocks_%j.out
set -euo pipefail
cd /global/scratch/users/tbellg/kmate
/global/home/users/tbellg/miniforge3/envs/kmate/bin/python -u analysis/grenenet_gea/gea_newpanel/sig_blocks_table.py
echo DONE_SIG_BLOCKS
