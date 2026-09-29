#!/bin/bash
# Refit the stratified arms (frequency, size) under THE model -- the full one, with local adaptation
# allowed its own strength in each garden: + (0 + obio1c | garden), fid pruned.
# Same stratification as _submit_freq_tiers.sh / _submit_size_bins.sh; the only change is the random
# effect, so §6 and §7 of the notebook report the same specification as every other section.
#SBATCH --job-name=burden_strata_lare
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1 --mem=16G
#SBATCH --time=03:00:00
#SBATCH --array=0-22
#SBATCH --output=logs/burden_strata_lare_%A_%a.out
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
export POOL_SET=sv6 BURDEN_ORIGIN=all ONESTAGE_DIR=onestage_all FINAL_OUT=final_models_pooled_all
export EXTRA_RE="(0 + obio1c | garden)"
export RE_DROP=fid
export STRATIFY_CTX=CDS
i=$SLURM_ARRAY_TASK_ID
if   [ "$i" -lt 7 ];  then   # freq tier 1: all 231 founders, rare/common
  export BURDEN_STRATIFY=freq FREQ_BINS="rare,common" SUBSET=""
  export EXTRA_FIXED="obio1c + obio1c:bio1c + callrate_sv"; export TASK=$i
elif [ "$i" -lt 15 ]; then   # freq tier 2: 80 cactus founders, private/rare/common
  export BURDEN_STRATIFY=freq FREQ_BINS="private,rare,common" SUBSET=cactus
  export EXTRA_FIXED="obio1c + obio1c:bio1c + callrate_sv"; export TASK=$(( i - 7 ))
else                         # size bins
  export BURDEN_STRATIFY=size FREQ_BINS="small,medium,large" SUBSET=""
  export EXTRA_FIXED="obio1c + obio1c:bio1c"; export TASK=$(( i - 15 ))
fi
echo "[$(date)] host=$(hostname) STRATIFY=$BURDEN_STRATIFY BINS=$FREQ_BINS SUBSET=${SUBSET:-all} TASK=$TASK"
/global/home/users/tbellg/miniforge3/envs/r_env/bin/Rscript _sv_burden_pooled_models.R
