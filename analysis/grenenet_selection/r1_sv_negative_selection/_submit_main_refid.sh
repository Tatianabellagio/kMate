#!/bin/bash
# The main model with the independent founder intercept pruned (RE_DROP=fid).
# In the standing fit that variance component sits at 4.9e-07 -- it is estimating nothing, so all
# founder-level variation is already carried by the kinship term. Removing it costs one parameter and
# gives a specification with no component at the boundary. Outputs carry the _refid suffix, so the
# standing fit is not touched and the two can be compared directly.
#SBATCH --job-name=burden_refid
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1 --mem=16G
#SBATCH --time=03:00:00
#SBATCH --array=0-5
#SBATCH --output=logs/burden_refid_%A_%a.out
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
export POOL_SET=sv6 BURDEN_ORIGIN=all ONESTAGE_DIR=onestage_all FINAL_OUT=final_models_pooled_all
export EXTRA_FIXED="obio1c + obio1c:bio1c"
export RE_DROP=fid
export TASK=$SLURM_ARRAY_TASK_ID
echo "[$(date)] host=$(hostname) RE_DROP=$RE_DROP TASK=$TASK"
/global/home/users/tbellg/miniforge3/envs/r_env/bin/Rscript _sv_burden_pooled_models.R
