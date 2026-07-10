#!/bin/bash
#SBATCH --job-name=metric_proto
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=1:00:00
#SBATCH --output=logs/metric_proto_%j.out
#SBATCH --error=logs/metric_proto_%j.out
# Two clean per-region metrics: panel-only resolvability + analytic coverage.
set -euo pipefail
cd /global/scratch/users/tbellg/kmate
mkdir -p logs benchmarks/h_uncertainty/results
PY=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python
$PY benchmarks/h_uncertainty/metric_prototype.py
echo "=== done $(date) ==="
