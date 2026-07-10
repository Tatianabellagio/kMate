#!/bin/bash
#SBATCH --job-name=kb_diag
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G
#SBATCH --time=1:00:00
#SBATCH --output=/global/scratch/users/tbellg/kmate/benchmarks/p231/logs/kb_diag_%j.out
#SBATCH --error=/global/scratch/users/tbellg/kmate/benchmarks/p231/logs/kb_diag_%j.out
set -euo pipefail
source ~/miniforge3/etc/profile.d/conda.sh && conda activate kmate
cd /global/scratch/users/tbellg/kmate
python benchmarks/p231/scripts/kb_per_block.py
echo DONE-KB
