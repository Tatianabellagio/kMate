#!/bin/bash
#SBATCH --job-name=classmat123
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --requeue
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=4:00:00
#SBATCH --output=/global/scratch/projects/fc_moilab/tbellg/kmate/analysis/grenenet_gea/phase1_replication/logs/classmat123_%j.out
#SBATCH --error=/global/scratch/projects/fc_moilab/tbellg/kmate/analysis/grenenet_gea/phase1_replication/logs/classmat123_%j.out
set -euo pipefail
cd /global/scratch/projects/fc_moilab/tbellg/kmate/analysis/grenenet_gea/phase1_replication
PY=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python
echo "== host $(hostname) =="
$PY -u build_class_matrices.py --classes snp sv smallindel --gens 1 2 3
echo "== ALLDONE_CLASSMAT123 =="
