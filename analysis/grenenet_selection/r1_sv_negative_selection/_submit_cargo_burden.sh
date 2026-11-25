#!/bin/bash
# Cargo-aware burden model: the REPORTED model (sv_burden_model.ipynb; `_submit_refit_ctx.sh` tasks 0-8) with
# its six location contexts replaced by content x location -- TE / nTE x CDS, UTR, intron, promoter,
# intergenic, insertions and deletions pooled (`_sv_burden_cargo_build.py` -> `_sv_burden_cargo_export.py`).
# Every other setting is the reported one: origin all, rate burdens, local-adaptation control
# obio1c + obio1c:bio1c with its own garden slope, no independent founder intercept, every burden with its
# own temperature interaction and random garden slope.
# Tasks: 0 M1p_reml, 1 M2p_reml, 2 M2p_ml, 3-12 M2p_ml_drop_<term> (LRT per interaction).
#SBATCH --job-name=burden_cargo
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1 --mem=16G
#SBATCH --time=08:00:00
#SBATCH --array=0-12
#SBATCH --exclude=n0041.savio4
#SBATCH --output=logs/burden_cargo_%A_%a.out
set -euo pipefail
# absolute, not $SLURM_SUBMIT_DIR: from inside an salloc shell that points at wherever salloc was run
cd /global/scratch/projects/fc_moilab/tbellg/kmate/analysis/grenenet_selection/r1_sv_negative_selection
export BURDEN_ORIGIN=all ONESTAGE_DIR=onestage_all FINAL_OUT=final_models_cargo
export EXTRA_RE="(0 + obio1c | garden)" RE_DROP=fid CRTAG="" BURDEN_MODE=rate
export POOL_SET=cargo10 INTER_SET=all EXTRA_FIXED="obio1c + obio1c:bio1c"
export BURDEN_FILE=founder_burden_event_cargo.csv
export TASK=$SLURM_ARRAY_TASK_ID
echo "[$(date)] host=$(hostname) POOL_SET=$POOL_SET TASK=$TASK BURDEN_FILE=$BURDEN_FILE"
/global/home/users/tbellg/miniforge3/envs/r_env/bin/Rscript _sv_burden_pooled_models.R
