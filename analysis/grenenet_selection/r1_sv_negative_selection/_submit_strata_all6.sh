#!/bin/bash
# Frequency and size strata under THE model as it now stands: all six genomic contexts as their own
# burdens, each with its own garden-temperature interaction, local adaptation per garden, fid pruned,
# 0.9 call-rate floor, carried/called. Same stratification as before; only the parameterisation changes.
#SBATCH --job-name=burden_strata_all6
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1 --mem=16G
#SBATCH --time=03:00:00
#SBATCH --array=0-31
#SBATCH --output=logs/burden_strata_all6_%A_%a.out
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
export POOL_SET=all6 INTER_SET=all BURDEN_ORIGIN=all ONESTAGE_DIR=onestage_all
export FINAL_OUT=final_models_pooled_all
export EXTRA_RE="(0 + obio1c | garden)" RE_DROP=fid STRATIFY_CTX=CDS CRTAG=""
i=$SLURM_ARRAY_TASK_ID
if   [ "$i" -lt 10 ]; then   # freq tier 1: all 231 founders, rare/common -> 7 inter terms -> 10 tasks
  export BURDEN_STRATIFY=freq FREQ_BINS="rare,common" SUBSET=""
  export EXTRA_FIXED="obio1c + obio1c:bio1c + callrate_sv"; export TASK=$i
elif [ "$i" -lt 21 ]; then   # freq tier 2: 80 cactus founders, private/rare/common -> 8 -> 11 tasks
  export BURDEN_STRATIFY=freq FREQ_BINS="private,rare,common" SUBSET=cactus
  export EXTRA_FIXED="obio1c + obio1c:bio1c + callrate_sv"; export TASK=$(( i - 10 ))
else                         # size bins -> 8 -> 11 tasks
  export BURDEN_STRATIFY=size FREQ_BINS="small,medium,large" SUBSET=""
  export EXTRA_FIXED="obio1c + obio1c:bio1c"; export TASK=$(( i - 21 ))
fi
echo "[$(date)] host=$(hostname) all6 STRATIFY=$BURDEN_STRATIFY BINS=$FREQ_BINS SUBSET=${SUBSET:-all} TASK=$TASK"
/global/home/users/tbellg/miniforge3/envs/r_env/bin/Rscript _sv_burden_pooled_models.R
