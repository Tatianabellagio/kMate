#!/bin/bash
#SBATCH --job-name=wza_3model
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=2:00:00
#SBATCH --output=/global/scratch/users/tbellg/kmate/analysis/grenenet_gea/gea_newpanel/_regen2_logs/wza_3model_%j.log
#SBATCH --error=/global/scratch/users/tbellg/kmate/analysis/grenenet_gea/gea_newpanel/_regen2_logs/wza_3model_%j.log
# Fresh 3-model (binomial/kendall/lfmm) x 2-class WZA on the rerun_kfw_hb cohort,
# clq0.9 blocks, deg7-cap2000. Feeds the unified SNP-vs-nonSNP 3x2 raw+WZA grid.
set -euo pipefail
cd /global/scratch/users/tbellg/kmate/analysis/grenenet_gea/gea_newpanel
PYK=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python
echo "== host $(hostname) =="
$PYK -u run_wza_3model_clq09.py
echo "== ALLDONE_WZA_3MODEL =="
