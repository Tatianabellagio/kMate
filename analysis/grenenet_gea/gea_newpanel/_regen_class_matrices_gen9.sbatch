#!/bin/bash
#SBATCH --job-name=class_mat_g9
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --account=co_moilab
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G
#SBATCH --time=4:00:00
#SBATCH --output=/global/scratch/users/tbellg/kmate/results/grenenet_gea/gea_newpanel/_regen_logs/class_matrices_g9_%j.log

set -euo pipefail
cd /global/scratch/users/tbellg/kmate
PY=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python
# Rebuild gen9 (last-gen) class matrices from the FRESH pool_gen9 (assemble_lastgen)
# + fresh af_store index. Classes used by gea_newpanel: snp, nonsnp, sv, smallindel.
$PY -u analysis/grenenet_gea/phase1_replication/build_class_matrices.py \
    --gens 9 --classes snp nonsnp sv smallindel
echo DONE_CLASS_MATRICES_G9
