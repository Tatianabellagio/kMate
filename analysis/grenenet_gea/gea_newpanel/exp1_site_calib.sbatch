#!/bin/bash
#SBATCH --job-name=exp1_site_calib
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --cpus-per-task=8
#SBATCH --mem=192G
#SBATCH --time=6:00:00
#SBATCH --dependency=afterok:35492923
#SBATCH --output=/global/scratch/users/tbellg/kmate/analysis/grenenet_gea/gea_newpanel/exp1_site_calib_%j.out
#SBATCH --error=/global/scratch/users/tbellg/kmate/analysis/grenenet_gea/gea_newpanel/exp1_site_calib_%j.out
set -euo pipefail
cd /global/scratch/users/tbellg/kmate
R=/global/home/users/tbellg/miniforge3/envs/lfmm_env/bin/Rscript
PY=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python
NP=analysis/grenenet_gea/gea_newpanel
SITE=analysis/grenenet_gea/gea_newpanel/results/lfmm_site
ES=analysis/grenenet_gea/gea_newpanel/results/env_site
EP=analysis/grenenet_gea/gea_newpanel/results/env_plot
PLOT=analysis/grenenet_gea/phase1_replication/results/lfmm
POW=analysis/grenenet_gea/gea_newpanel/results/power
mkdir -p "$POW"
echo "=== EXP1 site-level bio1 K sweep (raw vs GIF) ==="
for CLS in snp nonsnp; do for K in 1 3 5; do
  $R $NP/run_lfmm_both_p.R "$SITE/lfmm_${CLS}_site" "$ES/env_site_bio1.csv" $K \
     "$POW/exp1_site_${CLS}_bio1_K${K}_bothp.csv"
done; done
$PY -u $NP/report_power.py $POW/exp1_*_bothp.csv
echo DONE_EXP1
