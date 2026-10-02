#!/bin/bash
# THE model fitted on each sequencing-technology panel separately, as a robustness check
# alongside the frequency and size strata: cactus-assembled founders only, short-read
# (PanGenie-genotyped) founders only, and all 231 together.
# `cactus` drops out of the fixed effects inside a panel (constant), handled by the fitter.
#SBATCH --job-name=burden_tech
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1 --mem=16G
#SBATCH --time=04:00:00
#SBATCH --array=0-17
#SBATCH --output=logs/burden_tech_%A_%a.out
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
export POOL_SET=all6 INTER_SET=all BURDEN_ORIGIN=all ONESTAGE_DIR=onestage_all
export FINAL_OUT=final_models_pooled_all BURDEN_MODE=rate CRTAG=""
export EXTRA_FIXED="obio1c + obio1c:bio1c"
export EXTRA_RE="(0 + obio1c | garden)" RE_DROP=fid
i=$SLURM_ARRAY_TASK_ID
if [ "$i" -lt 9 ]; then export SUBSET=cactus;    export TASK=$i
else                   export SUBSET=shortread; export TASK=$(( i - 9 )); fi
echo "[$(date)] SUBSET=$SUBSET TASK=$TASK"
/global/home/users/tbellg/miniforge3/envs/r_env/bin/Rscript _sv_burden_pooled_models.R
