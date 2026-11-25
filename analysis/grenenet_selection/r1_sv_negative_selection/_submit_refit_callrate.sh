#!/bin/bash
# Refit every arm sv_burden_callrate.ipynb reads, under the 1 kb promoter rule.
#SBATCH --job-name=refit_cr
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1 --mem=16G
#SBATCH --time=04:00:00
#SBATCH --array=0-72
#SBATCH --output=logs/refit_cr_%A_%a.out
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
export BURDEN_ORIGIN=all ONESTAGE_DIR=onestage_all FINAL_OUT=final_models_pooled_all
export EXTRA_RE="(0 + obio1c | garden)" RE_DROP=fid
E0="obio1c + obio1c:bio1c"; EC="$E0 + callrate_sv"
i=$SLURM_ARRAY_TASK_ID
if   [ "$i" -lt 18 ]; then           # sv6 RATE at 0.9 / 0.7 / 0.5
  export POOL_SET=sv6 BURDEN_MODE=rate EXTRA_FIXED="$E0"
  b=$(( i / 6 )); export TASK=$(( i % 6 ))
  case $b in 0) export CRTAG="";; 1) export CRTAG=_cr70;; 2) export CRTAG=_cr50;; esac
elif [ "$i" -lt 36 ]; then           # sv6 COUNT at 0.9 / 0.7 / 0.5
  j=$(( i - 18 )); export POOL_SET=sv6 BURDEN_MODE=count EXTRA_FIXED="$E0"
  b=$(( j / 6 )); export TASK=$(( j % 6 ))
  case $b in 0) export CRTAG="";; 1) export CRTAG=_cr70;; 2) export CRTAG=_cr50;; esac
elif [ "$i" -lt 48 ]; then           # sv6 COUNT + callrate_sv at 0.7 / 0.5
  j=$(( i - 36 )); export POOL_SET=sv6 BURDEN_MODE=count EXTRA_FIXED="$EC"
  b=$(( j / 6 )); export TASK=$(( j % 6 ))
  case $b in 0) export CRTAG=_cr70;; 1) export CRTAG=_cr50;; esac
elif [ "$i" -lt 55 ]; then           # call-rate tiers
  export POOL_SET=sv6 BURDEN_MODE=rate CRTAG="" EXTRA_FIXED="$E0"
  export BURDEN_STRATIFY=callrate STRATIFY_CTX=CDS FREQ_BINS="lowcall,highcall"
  export TASK=$(( i - 48 ))
elif [ "$i" -lt 64 ]; then           # tier x gene density
  export POOL_SET=sv6 BURDEN_MODE=rate CRTAG="" EXTRA_FIXED="$E0"
  export BURDEN_STRATIFY=tierdens STRATIFY_CTX=CDS FREQ_BINS="hi_rich,hi_poor,lo_rich,lo_poor"
  export TASK=$(( i - 55 ))
else                                 # tier x size
  export POOL_SET=sv6 BURDEN_MODE=rate CRTAG="" EXTRA_FIXED="$E0"
  export BURDEN_STRATIFY=tiersize STRATIFY_CTX=CDS FREQ_BINS="hi_small,hi_big,lo_small,lo_big"
  export TASK=$(( i - 64 ))
fi
echo "[$(date)] MODE=$BURDEN_MODE CRTAG=${CRTAG:-0.9} STRAT=${BURDEN_STRATIFY:-none} TASK=$TASK"
/global/home/users/tbellg/miniforge3/envs/r_env/bin/Rscript _sv_burden_pooled_models.R
