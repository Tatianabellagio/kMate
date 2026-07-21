#!/bin/bash
#SBATCH --job-name=exp3_plot_K_pc1
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --cpus-per-task=8
#SBATCH --mem=192G
#SBATCH --time=6:00:00
#SBATCH --dependency=afterok:35492923
#SBATCH --output=/global/scratch/users/tbellg/kmate/analysis/grenenet_gea/phase1_replication/axis_scan/exp3_plot_K_pc1_%j.out
#SBATCH --error=/global/scratch/users/tbellg/kmate/analysis/grenenet_gea/phase1_replication/axis_scan/exp3_plot_K_pc1_%j.out
set -euo pipefail
cd /global/scratch/users/tbellg/kmate
R=/global/home/users/tbellg/miniforge3/envs/lfmm_env/bin/Rscript
PY=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python
NP=analysis/grenenet_gea/phase1_replication/axis_scan
SITE=analysis/grenenet_gea/phase1_replication/results/lfmm_site
ES=analysis/grenenet_gea/phase1_replication/results/env_site
EP=analysis/grenenet_gea/phase1_replication/results/env_plot
PLOT=analysis/grenenet_gea/phase1_replication/results/lfmm
POW=analysis/grenenet_gea/phase1_replication/results/power
mkdir -p "$POW"
echo "=== EXP3 plot-level bio1 K16/K20 + pc1 (raw vs GIF) ==="
for CLS in snp nonsnp; do
  for K in 16 20; do
    $R $NP/run_lfmm_both_p.R "$PLOT/lfmm_${CLS}_gen9" "$PLOT/lfmm_${CLS}_gen9_env.csv" $K \
       "$POW/exp3_plot_${CLS}_bio1_K${K}_bothp.csv"
  done
  $R $NP/run_lfmm_both_p.R "$PLOT/lfmm_${CLS}_gen9" "$EP/env_plot_pc1.csv" 16 \
     "$POW/exp3_plot_${CLS}_pc1_K16_bothp.csv"
done
$PY -u $NP/report_power.py $POW/exp3_*_bothp.csv
echo DONE_EXP3
