#!/bin/bash
#SBATCH --job-name=chr1_truthest
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=00:20:00
#SBATCH --output=logs/E1_truthest_%j.out
#SBATCH --error=logs/E1_truthest_%j.err
mkdir -p logs
set -euo pipefail

# Recreate AF_TRUTH_VS_ESTIMATE_v3qc_v3_mixedloose for Arch 3.
# Inputs already on disk; this is light compute (sparse @ vector + pandas joins).
cd /global/scratch/users/tbellg/kmate/panel/arch3/chr1
PY=/global/home/users/tbellg/miniforge3/envs/hapfm/bin/python
$PY -u af_truth_vs_estimate_arch3.py
echo "[$(date)] DONE E1"
