#!/bin/bash
#SBATCH --job-name=kendall_hap
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=1:00:00
#SBATCH --output=analysis/grenenet_gea/phase1_replication/logs/kendall_hap_%j.out
#SBATCH --error=analysis/grenenet_gea/phase1_replication/logs/kendall_hap_%j.out
# Kendall-tau GEA on clq0.9 haploblock frequencies vs bio1, across the 355 gen9 pools.
set -euo pipefail
cd /global/scratch/users/tbellg/kmate
PY=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python
$PY -u analysis/grenenet_gea/phase1_replication/run_kendall.py \
    --class hap --gen 9 --climate bio1 --threads 8
