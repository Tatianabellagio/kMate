#!/bin/bash
# Stability of the main burden result against the founder call-rate floor.
# _sv_burden_build.py rebuilds the per-founder burdens at SV_CALL_MIN = 0.7 and 0.5 (0.9 is the reported
# one). Only the SV burden table changes: fitness, the SNP coding burden and the kinship GRM are all
# rebuilt from SNPs held at the 0.9 floor, so the arms differ in exactly one thing.
#SBATCH --job-name=burden_cr_build
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1 --mem=64G
#SBATCH --time=04:00:00
#SBATCH --array=0-1
#SBATCH --output=logs/burden_cr_build_%A_%a.out
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
if [ "$SLURM_ARRAY_TASK_ID" -eq 0 ]; then export SV_CALL_MIN=0.7; else export SV_CALL_MIN=0.5; fi
export SV_UNIT=event
echo "[$(date)] host=$(hostname) SV_CALL_MIN=$SV_CALL_MIN SV_UNIT=event"
/global/home/users/tbellg/miniforge3/envs/kmate/bin/python _sv_burden_build.py
