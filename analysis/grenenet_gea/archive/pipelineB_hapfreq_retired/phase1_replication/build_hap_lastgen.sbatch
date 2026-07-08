#!/bin/bash
#SBATCH --job-name=hap_lastgen
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=2:00:00
#SBATCH --output=analysis/grenenet_gea/phase1_replication/logs/hap_lastgen_%j.out
#SBATCH --error=analysis/grenenet_gea/phase1_replication/logs/hap_lastgen_%j.out
# Build the clq0.9 haploblock-frequency last-gen pool matrix (global-h projection ->
# 355 gen9 pools) for phase-1-style GEA. I/O-bound on 907 member-sample window-h reads.
set -euo pipefail
cd /global/scratch/users/tbellg/kmate
mkdir -p analysis/grenenet_gea/phase1_replication/logs
PY=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python
$PY -u analysis/grenenet_gea/phase1_replication/build_hap_lastgen_matrix.py
