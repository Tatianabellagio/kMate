#!/bin/bash
# Add the bio1c MAIN effect. The model fits sv_ctx:bio1c and obio1c:bio1c but has no bio1c main term --
# the garden random intercept absorbs any garden-level mean, including a temperature trend. That is
# defensible but it violates marginality (interactions fitted without their main effect), so it is worth
# measuring rather than assuming.
#SBATCH --job-name=bio1c_main
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1 --mem=16G
#SBATCH --time=03:00:00
#SBATCH --array=0-8
#SBATCH --output=logs/bio1c_main_%A_%a.out
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
export POOL_SET=all6 INTER_SET=all BURDEN_ORIGIN=all ONESTAGE_DIR=onestage_all
export FINAL_OUT=final_models_pooled_all BURDEN_MODE=rate CRTAG=""
export EXTRA_FIXED="obio1c + obio1c:bio1c + bio1c"
export EXTRA_RE="(0 + obio1c | garden)" RE_DROP=fid
export TASK=$SLURM_ARRAY_TASK_ID
/global/home/users/tbellg/miniforge3/envs/r_env/bin/Rscript _sv_burden_pooled_models.R
