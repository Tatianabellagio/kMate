#!/bin/bash
#SBATCH --job-name=chr1_truthest
#SBATCH --partition=bse
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=00:20:00
#SBATCH --output=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/arch3/chr1/E1_truthest_%j.out
#SBATCH --error=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/arch3/chr1/E1_truthest_%j.err
set -euo pipefail

# Recreate AF_TRUTH_VS_ESTIMATE_v3qc_v3_mixedloose for Arch 3.
# Inputs already on disk; this is light compute (sparse @ vector + pandas joins).
cd /carnegie/nobackup/scratch/tbellagio/hapfire_sv/arch3/chr1
PY=/home/tbellagio/miniforge3/envs/hapfm/bin/python
$PY -u af_truth_vs_estimate_arch3.py
echo "[$(date)] DONE E1"
