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
# GEA step 1 (_submit_gea.sh): the generation-9 class matrices of the record arm (snp, smallindel), filtered to
# MAF >= 0.05 and founder call rate >= 0.50 (build_class_matrices.py defaults) -> results/class_matrices_cr50_arch3v2/
set -euo pipefail
cd /global/scratch/projects/fc_moilab/tbellg/kmate
PR=analysis/grenenet_selection/r2_gea_nonsnp/phase1_replication
V=${KMATE_PANEL:-arch3v2}; S=""; [ "$V" = arch3 ] || S="_$V"     # panel version (lib.PANEL): outputs under *_arch3v2 dirs
/global/home/users/tbellg/miniforge3/envs/kmate/bin/python -u $PR/build_class_matrices.py \
    --classes snp smallindel --gens 9 --out $PR/results/class_matrices_cr50$S
echo "CM_CR50_DONE"
