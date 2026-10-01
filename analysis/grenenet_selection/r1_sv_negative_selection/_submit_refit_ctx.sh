#!/bin/bash
# Refit every arm the two main notebooks read, on the burdens rebuilt under the 1 kb promoter rule.
#SBATCH --job-name=refit_ctx
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1 --mem=16G
#SBATCH --time=04:00:00
#SBATCH --array=0-45
#SBATCH --output=logs/refit_ctx_%A_%a.out
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
export BURDEN_ORIGIN=all ONESTAGE_DIR=onestage_all FINAL_OUT=final_models_pooled_all
export EXTRA_RE="(0 + obio1c | garden)" RE_DROP=fid CRTAG="" BURDEN_MODE=rate
i=$SLURM_ARRAY_TASK_ID
if   [ "$i" -lt 9 ]; then            # all6 + every context interacted  (THE model)
  export POOL_SET=all6 INTER_SET=all EXTRA_FIXED="obio1c + obio1c:bio1c"; export TASK=$i
elif [ "$i" -lt 14 ]; then           # cds2
  export POOL_SET=cds2 EXTRA_FIXED="obio1c + obio1c:bio1c"; export TASK=$(( i - 9 ))
else                                 # strata under all6
  export POOL_SET=all6 INTER_SET=all STRATIFY_CTX=CDS
  j=$(( i - 14 ))
  if   [ "$j" -lt 10 ]; then
    export BURDEN_STRATIFY=freq FREQ_BINS="rare,common" SUBSET=""
    export EXTRA_FIXED="obio1c + obio1c:bio1c + callrate_sv"; export TASK=$j
  elif [ "$j" -lt 21 ]; then
    export BURDEN_STRATIFY=freq FREQ_BINS="private,rare,common" SUBSET=cactus
    export EXTRA_FIXED="obio1c + obio1c:bio1c + callrate_sv"; export TASK=$(( j - 10 ))
  else
    export BURDEN_STRATIFY=size FREQ_BINS="small,medium,large" SUBSET=""
    export EXTRA_FIXED="obio1c + obio1c:bio1c"; export TASK=$(( j - 21 ))
  fi
fi
echo "[$(date)] POOL_SET=$POOL_SET STRATIFY=${BURDEN_STRATIFY:-none} TASK=$TASK"
/global/home/users/tbellg/miniforge3/envs/r_env/bin/Rscript _sv_burden_pooled_models.R
