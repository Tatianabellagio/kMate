#!/bin/bash
# Fit THE model with the coding burden split into the two call-rate tiers, both present at once:
#   sv_CDS_highcall  the events the reported 0.9 model uses
#   sv_CDS_lowcall   the events relaxing the floor newly admits
#SBATCH --job-name=burden_crtier_fit
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1 --mem=16G
#SBATCH --time=03:00:00
#SBATCH --array=0-6
#SBATCH --output=logs/burden_crtier_fit_%A_%a.out
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
export POOL_SET=sv6 BURDEN_ORIGIN=all ONESTAGE_DIR=onestage_all FINAL_OUT=final_models_pooled_all
export BURDEN_STRATIFY=callrate STRATIFY_CTX=CDS FREQ_BINS="lowcall,highcall"
export EXTRA_FIXED="obio1c + obio1c:bio1c"
export EXTRA_RE="(0 + obio1c | garden)"
export RE_DROP=fid
export TASK=$SLURM_ARRAY_TASK_ID
echo "[$(date)] host=$(hostname) call-rate tiers, TASK=$TASK"
/global/home/users/tbellg/miniforge3/envs/r_env/bin/Rscript _sv_burden_pooled_models.R
