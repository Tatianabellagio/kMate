#!/bin/bash
# Why do rate and count disagree on the INTERACTIONS at the 0.9 floor, where the two predictors
# correlate at 0.97-0.99? The part where they differ is the denominator, which correlates with founder
# SV call rate (0.79-0.84), assembly route (0.39-0.54) and founder ORIGIN CLIMATE (-0.20 to -0.33).
# The model fits obio1c:bio1c but has no cactus:bio1c and no callrate_sv:bio1c, so a burden variable
# carrying either can show a spurious burden x temperature term. Fit those interactions explicitly.
#SBATCH --job-name=burden_techint
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1 --mem=16G
#SBATCH --time=03:00:00
#SBATCH --array=0-23
#SBATCH --output=logs/burden_techint_%A_%a.out
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
i=$SLURM_ARRAY_TASK_ID; blk=$(( i / 6 )); export TASK=$(( i % 6 ))
export POOL_SET=sv6 BURDEN_ORIGIN=all ONESTAGE_DIR=onestage_all FINAL_OUT=final_models_pooled_all
export EXTRA_RE="(0 + obio1c | garden)" RE_DROP=fid CRTAG=""
case $blk in
  0) export BURDEN_MODE=rate ; export EXTRA_FIXED="obio1c + obio1c:bio1c + cactus:bio1c" ;;
  1) export BURDEN_MODE=count; export EXTRA_FIXED="obio1c + obio1c:bio1c + cactus:bio1c" ;;
  2) export BURDEN_MODE=rate ; export EXTRA_FIXED="obio1c + obio1c:bio1c + callrate_sv + callrate_sv:bio1c" ;;
  3) export BURDEN_MODE=count; export EXTRA_FIXED="obio1c + obio1c:bio1c + callrate_sv + callrate_sv:bio1c" ;;
esac
echo "[$(date)] host=$(hostname) MODE=$BURDEN_MODE EXTRA_FIXED='$EXTRA_FIXED' TASK=$TASK"
/global/home/users/tbellg/miniforge3/envs/r_env/bin/Rscript _sv_burden_pooled_models.R
