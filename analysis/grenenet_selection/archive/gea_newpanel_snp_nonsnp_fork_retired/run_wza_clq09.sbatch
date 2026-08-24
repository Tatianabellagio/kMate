#!/bin/bash
#SBATCH --job-name=wza_clq09
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G
#SBATCH --time=3:00:00
#SBATCH --output=/global/scratch/users/tbellg/kmate/analysis/grenenet_gea/gea_newpanel/wza_clq09_%j.out
#SBATCH --error=/global/scratch/users/tbellg/kmate/analysis/grenenet_gea/gea_newpanel/wza_clq09_%j.out
# Reassign clq0.9 blocks to the raw LFMM p-values and run canonical WZA per class.
set -euo pipefail
cd /global/scratch/users/tbellg/kmate
PY=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python
$PY -u analysis/grenenet_gea/gea_newpanel/run_wza_clq09.py
echo DONE_WZA_CLQ09
