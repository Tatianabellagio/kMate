#!/bin/bash
# Size-stratified CDS burden on the main model: pooled SVs (ins/del is reference-dependent and therefore
# not a biological split), no origin filter, reference TE kept, local-adaptation control.
# CDS split into small (50-99 bp) / medium (100-499) / large (500+), the most balanced 3-bin split
# available (444/481/482 CDS events). MAC >= 2 throughout, matching the main model's `all` class.
#SBATCH --job-name=burden_size
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1 --mem=16G
#SBATCH --time=03:00:00
#SBATCH --array=0-7
#SBATCH --output=logs/burden_size_%A_%a.out
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
export POOL_SET=sv6 BURDEN_ORIGIN=all ONESTAGE_DIR=onestage_all FINAL_OUT=final_models_pooled_all
export BURDEN_STRATIFY=size STRATIFY_CTX=CDS FREQ_BINS="small,medium,large"
export EXTRA_FIXED="obio1c + obio1c:bio1c"
export TASK=$SLURM_ARRAY_TASK_ID
echo "[$(date)] host=$(hostname) STRATIFY=size TASK=$TASK"
/global/home/users/tbellg/miniforge3/envs/r_env/bin/Rscript _sv_burden_pooled_models.R
