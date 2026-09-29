#!/bin/bash
#SBATCH --job-name=classmat0
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --requeue
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=2:00:00
#SBATCH --output=/global/scratch/projects/fc_moilab/tbellg/kmate/analysis/grenenet_selection/r2_gea_nonsnp/phase1_replication/logs/classmat0_%j.out
#SBATCH --error=/global/scratch/projects/fc_moilab/tbellg/kmate/analysis/grenenet_selection/r2_gea_nonsnp/phase1_replication/logs/classmat0_%j.out
set -euo pipefail
cd /global/scratch/projects/fc_moilab/tbellg/kmate/analysis/grenenet_selection/r2_gea_nonsnp/phase1_replication
PY=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python
echo "== host $(hostname) =="
$PY -u build_class_matrices_seedmix.py --classes snp sv smallindel
echo "== ALLDONE_CLASSMAT0 =="
