#!/bin/bash
#SBATCH --job-name=ma_fresh
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --cpus-per-task=8
#SBATCH --mem=192G
#SBATCH --time=10:00:00
#SBATCH --array=0-17
#SBATCH --output=/global/scratch/users/tbellg/kmate/analysis/grenenet_gea/gea_newpanel/_regen2_logs/ma_fresh_%A_%a.log
#SBATCH --error=/global/scratch/users/tbellg/kmate/analysis/grenenet_gea/gea_newpanel/_regen2_logs/ma_fresh_%A_%a.log
# Multi-axis fresh GEA regen: one array task per bioclim axis bio2..bio19
# (18 axes; bio1 already done). Each task: binomial+kendall+lfmm x snp/nonsnp/sv + WZA.
set -uo pipefail
AXES=(bio2 bio3 bio4 bio5 bio6 bio7 bio8 bio9 bio10 bio11 bio12 bio13 bio14 bio15 bio16 bio17 bio18 bio19)
AXIS=${AXES[$SLURM_ARRAY_TASK_ID]}
bash /global/scratch/users/tbellg/kmate/analysis/grenenet_gea/gea_newpanel/multiaxis_fresh_run_axis.sh "$AXIS"
