#!/bin/bash
#SBATCH --job-name=cm_cr50
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --cpus-per-task=4
#SBATCH --mem=192G
#SBATCH --time=4:00:00
#SBATCH --output=/global/scratch/projects/fc_moilab/tbellg/kmate/analysis/grenenet_selection/r2_gea_nonsnp/phase1_replication/results/_build_cm_cr50_%j.out
#SBATCH --error=/global/scratch/projects/fc_moilab/tbellg/kmate/analysis/grenenet_selection/r2_gea_nonsnp/phase1_replication/results/_build_cm_cr50_%j.out
# Class matrices for gen9 WITH the founder call-rate floor of 0.50 (see build_class_matrices.py).
# Writes to results/class_matrices_cr50/ so the unfiltered production set stays intact and
# the two can be compared record for record.
set -euo pipefail
cd /global/scratch/projects/fc_moilab/tbellg/kmate
PY=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python
PR=analysis/grenenet_selection/r2_gea_nonsnp/phase1_replication
V=${KMATE_PANEL:-arch3v2}; S=""; [ "$V" = arch3 ] || S="_$V"     # panel version (lib.PANEL): outputs under *_arch3v2 dirs
$PY -u $PR/build_class_matrices.py \
    --classes snp sv smallindel nonsnp --gens 9 \
    --min-call-rate 0.50 \
    --out $PR/results/class_matrices_cr50$S
echo "CM_CR50_DONE"
