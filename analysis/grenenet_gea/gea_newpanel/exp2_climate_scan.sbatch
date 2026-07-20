#!/bin/bash
#SBATCH --job-name=exp2_climate_scan
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --cpus-per-task=8
#SBATCH --mem=192G
#SBATCH --time=6:00:00
#SBATCH --output=/global/scratch/users/tbellg/kmate/analysis/grenenet_gea/gea_newpanel/exp2_climate_scan_%j.out
#SBATCH --error=/global/scratch/users/tbellg/kmate/analysis/grenenet_gea/gea_newpanel/exp2_climate_scan_%j.out
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
echo "=== EXP2 site-level climate-axis scan bio1..19 + pc1 (K=3) ==="
for CLS in snp nonsnp; do for AX in $(seq -f 'bio%g' 1 19) pc1; do
  E="$ES/env_site_${AX}.csv"; [ -f "$E" ] || continue
  $R $NP/run_lfmm_both_p.R "$SITE/lfmm_${CLS}_site" "$E" 3 \
     "$POW/exp2_scan_${CLS}_${AX}_bothp.csv"
done; done
$PY -u $NP/report_power.py $POW/exp2_*_bothp.csv
echo DONE_EXP2
