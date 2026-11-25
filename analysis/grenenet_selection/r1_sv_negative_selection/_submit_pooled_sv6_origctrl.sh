#!/bin/bash
# The local-adaptation control on the POOLED, NO-ORIGIN-FILTER, TE-INCLUSIVE arm
# (POOL_SET=sv6 + BURDEN_ORIGIN=all -> final_models_pooled_all, the fits behind
# sv_burden_pooled_nofilter.ipynb). Its no-TE twin already has this control; this completes the pair.
#SBATCH --job-name=pooled_sv6_origctrl
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1 --mem=16G
#SBATCH --time=03:00:00
#SBATCH --array=0-5
#SBATCH --output=logs/pooled_sv6_origctrl_%A_%a.out
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
export POOL_SET=sv6
export BURDEN_ORIGIN=all
export ONESTAGE_DIR=onestage_all
export FINAL_OUT=final_models_pooled_all
export EXTRA_FIXED="obio1c + obio1c:bio1c"
export TASK=$SLURM_ARRAY_TASK_ID
echo "[$(date)] host=$(hostname) task=$TASK POOL_SET=$POOL_SET ORIGIN=$BURDEN_ORIGIN EXTRA='$EXTRA_FIXED'"
/global/home/users/tbellg/miniforge3/envs/r_env/bin/Rscript _sv_burden_pooled_models.R
