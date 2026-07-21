#!/bin/bash
#SBATCH --job-name=quasibinom
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=16
#SBATCH --mem=110G
#SBATCH --time=18:00:00
#SBATCH --array=0-3
#SBATCH --output=analysis/grenenet_gea/gea_newpanel/quasibinom_%A_%a.out
#SBATCH --error=analysis/grenenet_gea/gea_newpanel/quasibinom_%A_%a.out
# Overdispersion-corrected binomial GEA: (A) per-variant quasi-binomial + (B) effective-N,
# each WITH K latent factors. gen9 bio1. array: 0=snp/K16 1=nonsnp/K16 2=snp/K0 3=nonsnp/K0
set -euo pipefail
cd /global/scratch/users/tbellg/kmate
PY=/global/home/users/tbellg/miniforge3/envs/basic/bin/python
G=analysis/grenenet_gea/gea_newpanel
CLS=(snp nonsnp snp nonsnp); K=(16 16 0 0)
i=${SLURM_ARRAY_TASK_ID}
echo "== host $(hostname) | quasibinom class=${CLS[$i]} K=${K[$i]} gen9 bio1 =="
$PY -u $G/run_quasibinom_latent.py --class "${CLS[$i]}" --gen 9 --climate bio1 \
    --k "${K[$i]}" --threads "${SLURM_CPUS_PER_TASK:-16}"
echo "== ALLDONE_QUASIBINOM ${CLS[$i]} K${K[$i]} =="
