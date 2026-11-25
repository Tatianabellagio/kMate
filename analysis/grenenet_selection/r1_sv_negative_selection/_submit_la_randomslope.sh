#!/bin/bash
# THE TEST: does the burden x garden-temperature interaction survive when local adaptation is allowed its
# own strength in each garden, instead of being forced linear in garden temperature?
#   adds (0 + obio1c | garden) to the main model, keeping obio1c + obio1c:bio1c as fixed effects.
# Motivation: per-garden local-adaptation slopes correlate with garden temperature at only r = 0.39, and
# the per-garden burden slope correlates with the per-garden LA slope at r = 0.66 (nonCDS) / 0.47 (CDS).
#SBATCH --job-name=burden_lare
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1 --mem=16G
#SBATCH --time=03:00:00
#SBATCH --array=0-10
#SBATCH --output=logs/burden_lare_%A_%a.out
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
i=$SLURM_ARRAY_TASK_ID
if [ "$i" -lt 6 ]; then export POOL_SET=sv6;  export TASK=$i
else                   export POOL_SET=cds2; export TASK=$(( i - 6 )); fi
export BURDEN_ORIGIN=all ONESTAGE_DIR=onestage_all FINAL_OUT=final_models_pooled_all
export EXTRA_FIXED="obio1c + obio1c:bio1c"
export EXTRA_RE="(0 + obio1c | garden)"
export RE_DROP=fid
echo "[$(date)] host=$(hostname) POOL_SET=$POOL_SET TASK=$TASK EXTRA_RE='$EXTRA_RE'"
/global/home/users/tbellg/miniforge3/envs/r_env/bin/Rscript _sv_burden_pooled_models.R
