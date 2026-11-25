#!/bin/bash
# The missing technology-panel fit: the pooled no-filter no-TE model on the SHORT-READ founders only
# (151 of 231 -- the larger panel). A SUBSET=cactus fit (80 founders) already exists; this completes the
# pair so the interaction can be checked within each sequencing technology rather than across them.
#SBATCH --job-name=burden_shortread
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1 --mem=16G
#SBATCH --time=03:00:00
#SBATCH --array=0-4
#SBATCH --output=logs/pooled_shortread_%A_%a.out
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
export SUBSET=shortread
export POOL_SET=sv4
export BURDEN_ORIGIN=all
export ONESTAGE_DIR=onestage_all
export FINAL_OUT=final_models_pooled_all_noTE
export TASK=$SLURM_ARRAY_TASK_ID
echo "[$(date)] host=$(hostname) task=$TASK  SUBSET=$SUBSET POOL_SET=$POOL_SET BURDEN_ORIGIN=$BURDEN_ORIGIN"
/global/home/users/tbellg/miniforge3/envs/r_env/bin/Rscript _sv_burden_pooled_models.R
