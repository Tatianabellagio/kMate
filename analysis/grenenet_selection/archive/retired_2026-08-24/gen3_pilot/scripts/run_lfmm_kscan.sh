#!/bin/bash
#SBATCH --job-name=lfmm_kscan
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --cpus-per-task=8
#SBATCH --mem=192G
#SBATCH --time=8:00:00
#SBATCH --output=/global/scratch/users/tbellg/kmate/analysis/grenenet_selection/kscan_%j.out
#SBATCH --error=/global/scratch/users/tbellg/kmate/analysis/grenenet_selection/kscan_%j.out
set -euo pipefail
cd /global/scratch/users/tbellg/kmate/analysis/grenenet_selection
R=/global/home/users/tbellg/miniforge3/envs/lfmm_env/bin/Rscript
L=/global/scratch/users/tbellg/kmate/analysis/grenenet_selection/r2_gea_nonsnp/results/lfmm
DP=$L/delta_p_gen3_sv.csv
ENV=$L/env_gen3_bio1.csv
# per-K LFMM on the fresh gen3-SV Δp/env; writes kscan_bio1_k{K}.{calibrated_pval,pval,beta}.csv + .gif.txt
for K in 1 2 3 4 5 6 7 8; do
  echo "=== K=$K ==="
  $R r2_gea_nonsnp/run_lfmm.R "$DP" "$ENV" "$K" "$L/kscan_bio1_k${K}"
done
echo DONE_KSCAN
