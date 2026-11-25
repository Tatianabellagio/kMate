#!/bin/bash
# Re-run the REML fit of THE model (sv6 and cds2) only to write per-observation partial residuals, so the
# per-garden burden-vs-fitness figure comes from the fitted model rather than raw per-garden regressions.
# Identical specification to the reported fits -- same files are rewritten byte-for-byte, plus resid_*.csv.
#SBATCH --job-name=burden_resid
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1 --mem=16G
#SBATCH --time=03:00:00
#SBATCH --array=0-1
#SBATCH --output=logs/burden_resid_%A_%a.out
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
if [ "$SLURM_ARRAY_TASK_ID" -eq 0 ]; then export POOL_SET=sv6; else export POOL_SET=cds2; fi
export BURDEN_ORIGIN=all ONESTAGE_DIR=onestage_all FINAL_OUT=final_models_pooled_all
export EXTRA_FIXED="obio1c + obio1c:bio1c"
export EXTRA_RE="(0 + obio1c | garden)"
export RE_DROP=fid
export TASK=1          # M2p_reml -- the reported fit
echo "[$(date)] host=$(hostname) POOL_SET=$POOL_SET TASK=1 (residual export)"
/global/home/users/tbellg/miniforge3/envs/r_env/bin/Rscript _sv_burden_pooled_models.R
