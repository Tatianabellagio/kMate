#!/bin/bash
#SBATCH --job-name=lfmm_power
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --cpus-per-task=8
#SBATCH --mem=192G
#SBATCH --time=6:00:00
#SBATCH --output=/global/scratch/users/tbellg/kmate/analysis/grenenet_gea/gea_newpanel/lfmm_power_%j.out
#SBATCH --error=/global/scratch/users/tbellg/kmate/analysis/grenenet_gea/gea_newpanel/lfmm_power_%j.out
# Power experiments on the RAW LFMM: (1) site-level (31, honest unit) vs plot-level
# (355, pseudoreplicated); (2) GIF vs no-GIF; (3) climate axis scan bio1..19.
set -euo pipefail
cd /global/scratch/users/tbellg/kmate
PY=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python
R=/global/home/users/tbellg/miniforge3/envs/lfmm_env/bin/Rscript
NP=analysis/grenenet_gea/gea_newpanel
SITE=results/grenenet_gea/gea_newpanel/lfmm_site
POW=results/grenenet_gea/gea_newpanel/power
PLOT=results/grenenet_gea/phase1_replication/lfmm
mkdir -p "$POW"

echo "===== build site-level (flower-weighted per-site mean) matrices ====="
$PY -u $NP/build_site_matrices.py

echo "===== SITE-LEVEL bio1, K sweep (n=31) ====="
for CLS in snp nonsnp; do for K in 1 3 5; do
  $R $NP/run_lfmm_both_p.R "$SITE/lfmm_${CLS}_site" "$SITE/env_site_bio1.csv" $K \
     "$POW/site_${CLS}_bio1_K${K}_bothp.csv"
done; done

echo "===== SITE-LEVEL climate-axis scan bio1..19 (K=3) ====="
for CLS in snp nonsnp; do for i in $(seq 1 19); do
  E="$SITE/env_site_bio${i}.csv"; [ -f "$E" ] || continue
  $R $NP/run_lfmm_both_p.R "$SITE/lfmm_${CLS}_site" "$E" 3 \
     "$POW/site_${CLS}_bio${i}_K3_bothp.csv"
done; done

echo "===== PLOT-LEVEL bio1, K=16, raw vs GIF (pseudoreplicated baseline) ====="
for CLS in snp nonsnp; do
  $R $NP/run_lfmm_both_p.R "$PLOT/lfmm_${CLS}_gen9" "$PLOT/lfmm_${CLS}_gen9_env.csv" 16 \
     "$POW/plot_${CLS}_bio1_K16_bothp.csv"
done

echo "===== SUMMARY ====="
$PY -u $NP/report_power.py $POW/*_bothp.csv
echo DONE_LFMM_POWER
