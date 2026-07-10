#!/bin/bash
#SBATCH --job-name=af_cert
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=1:30:00
#SBATCH --output=logs/af_cert_%j.out
#SBATCH --error=logs/af_cert_%j.out
# Validate that the AF functional is unbiased + Fisher AF-SE is calibrated.
set -euo pipefail
cd /global/scratch/users/tbellg/kmate
mkdir -p logs benchmarks/h_uncertainty/results
PY=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python
$PY benchmarks/h_uncertainty/af_certainty.py
echo "=== done $(date) ==="
