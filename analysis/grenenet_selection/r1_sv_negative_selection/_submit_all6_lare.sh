#!/bin/bash
# The all6 arm under THE model: local adaptation fitted per garden, + (0 + obio1c | garden), fid pruned.
# cds2 and sv6 already have this fit (_submit_la_randomslope.sh); this closes the last arm so every
# section of both notebooks reports the same specification.
#SBATCH --job-name=burden_all6_lare
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1 --mem=16G
#SBATCH --time=03:00:00
#SBATCH --array=0-5
#SBATCH --output=logs/burden_all6_lare_%A_%a.out
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
export POOL_SET=all6 BURDEN_ORIGIN=all ONESTAGE_DIR=onestage_all FINAL_OUT=final_models_pooled_all
export EXTRA_FIXED="obio1c + obio1c:bio1c"
export EXTRA_RE="(0 + obio1c | garden)"
export RE_DROP=fid
export TASK=$SLURM_ARRAY_TASK_ID
echo "[$(date)] host=$(hostname) POOL_SET=all6 TASK=$TASK EXTRA_RE='$EXTRA_RE'"
/global/home/users/tbellg/miniforge3/envs/r_env/bin/Rscript _sv_burden_pooled_models.R
