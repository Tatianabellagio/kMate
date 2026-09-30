#!/bin/bash
# THE model refit on burdens built at a relaxed founder call-rate floor (0.7, 0.5) against the reported
# 0.9. Everything else is held: same fitness table, same kinship GRM (built from SNPs at the 0.9 floor),
# same covariates, same random effects including local adaptation per garden.
#SBATCH --job-name=burden_cr_fit
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1 --mem=16G
#SBATCH --time=03:00:00
#SBATCH --array=0-21
#SBATCH --output=logs/burden_cr_fit_%A_%a.out
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
i=$SLURM_ARRAY_TASK_ID
if [ "$i" -lt 11 ]; then export CRTAG=_cr70; else export CRTAG=_cr50; i=$(( i - 11 )); fi
if [ "$i" -lt 6 ]; then export POOL_SET=sv6;  export TASK=$i
else                   export POOL_SET=cds2; export TASK=$(( i - 6 )); fi
export BURDEN_ORIGIN=all ONESTAGE_DIR=onestage_all FINAL_OUT=final_models_pooled_all
export EXTRA_FIXED="obio1c + obio1c:bio1c"
export EXTRA_RE="(0 + obio1c | garden)"
export RE_DROP=fid
echo "[$(date)] host=$(hostname) CRTAG=$CRTAG POOL_SET=$POOL_SET TASK=$TASK"
/global/home/users/tbellg/miniforge3/envs/r_env/bin/Rscript _sv_burden_pooled_models.R
