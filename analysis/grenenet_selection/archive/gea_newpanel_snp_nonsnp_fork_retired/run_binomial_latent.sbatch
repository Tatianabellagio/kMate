#!/bin/bash
#SBATCH --job-name=binom_lf
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --time=10:00:00
#SBATCH --output=analysis/grenenet_gea/gea_newpanel/binom_lf_%j.out
#SBATCH --error=analysis/grenenet_gea/gea_newpanel/binom_lf_%j.out
# Latent-factor (K=16) binomial GEA, bio1, gen9, SNP vs non-SNP — LFMM-style structure covariate.
set -euo pipefail
cd /global/scratch/users/tbellg/kmate
PY=/global/home/users/tbellg/miniforge3/envs/basic/bin/python   # statsmodels
G=analysis/grenenet_gea/gea_newpanel
echo "== host $(hostname) | latent-factor binomial K=16 gen9 bio1 =="
for CLS in nonsnp snp; do
    $PY -u $G/run_binomial_latent.py --class "$CLS" --gen 9 --climate bio1 --k 16 \
        --threads "${SLURM_CPUS_PER_TASK:-8}"
done
echo "== ALLDONE_BINOM_LF =="
