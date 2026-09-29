#!/bin/bash
# Marginality check: the FULL model with the origin-climate INTERACTION but WITHOUT its main effect
# (EXTRA_FIXED="obio1c:bio1c" instead of "obio1c + obio1c:bio1c"). Tests whether keeping obio1c matters
# for the SV conclusions. Outputs carry the _x_obio1cbio1c suffix -- distinct from the FULL fit's
# _x_obio1cobio1cbio1c -- so nothing is overwritten.
#SBATCH --job-name=pooled_sv6_interonly
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1 --mem=16G
#SBATCH --time=03:00:00
#SBATCH --array=0-5
#SBATCH --output=logs/pooled_sv6_interonly_%A_%a.out
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
export POOL_SET=sv6
export BURDEN_ORIGIN=all
export ONESTAGE_DIR=onestage_all
export FINAL_OUT=final_models_pooled_all
export EXTRA_FIXED="obio1c:bio1c"
export TASK=$SLURM_ARRAY_TASK_ID
echo "[$(date)] host=$(hostname) task=$TASK EXTRA='$EXTRA_FIXED'"
/global/home/users/tbellg/miniforge3/envs/r_env/bin/Rscript _sv_burden_pooled_models.R
