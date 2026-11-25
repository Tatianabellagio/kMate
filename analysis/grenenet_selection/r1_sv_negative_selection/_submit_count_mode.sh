#!/bin/bash
# Burden as the raw CARRIED COUNT instead of carried/called, at each call-rate floor.
# Motivation: below the 0.9 floor the denominator is founder-specific (CV 0.177, r = 0.89 with the
# numerator), so the rate may be measuring genotypability. Dropping the denominator removes that -- but
# it also stops correcting for how many events a founder is genotyped at, so the last two arms add
# callrate_sv as an explicit covariate to see whether that recovers anything.
#SBATCH --job-name=burden_count
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1 --mem=16G
#SBATCH --time=03:00:00
#SBATCH --array=0-29
#SBATCH --output=logs/burden_count_%A_%a.out
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
i=$SLURM_ARRAY_TASK_ID
blk=$(( i / 6 )); export TASK=$(( i % 6 ))
export BURDEN_MODE=count POOL_SET=sv6
export BURDEN_ORIGIN=all ONESTAGE_DIR=onestage_all FINAL_OUT=final_models_pooled_all
export EXTRA_RE="(0 + obio1c | garden)" RE_DROP=fid
case $blk in
  0) export CRTAG="";      export EXTRA_FIXED="obio1c + obio1c:bio1c" ;;
  1) export CRTAG=_cr70;   export EXTRA_FIXED="obio1c + obio1c:bio1c" ;;
  2) export CRTAG=_cr50;   export EXTRA_FIXED="obio1c + obio1c:bio1c" ;;
  3) export CRTAG=_cr70;   export EXTRA_FIXED="obio1c + obio1c:bio1c + callrate_sv" ;;
  4) export CRTAG=_cr50;   export EXTRA_FIXED="obio1c + obio1c:bio1c + callrate_sv" ;;
esac
echo "[$(date)] host=$(hostname) COUNT mode CRTAG=${CRTAG:-0.9} EXTRA_FIXED='$EXTRA_FIXED' TASK=$TASK"
/global/home/users/tbellg/miniforge3/envs/r_env/bin/Rscript _sv_burden_pooled_models.R
