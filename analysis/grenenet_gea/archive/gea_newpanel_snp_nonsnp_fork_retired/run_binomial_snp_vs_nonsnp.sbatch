#!/bin/bash
#SBATCH --job-name=binom_snpnon
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --time=8:00:00
#SBATCH --output=analysis/grenenet_gea/gea_newpanel/binom_snpnon_%j.out
#SBATCH --error=analysis/grenenet_gea/gea_newpanel/binom_snpnon_%j.out
# Raw binomial-regression climate-GEA (bio1, gen9), SNP vs non-SNP, new-panel comparison.
# Per-record binomial GLM [alt,ref] ~ const + z(bio1) across the 355 gen9 pools. No WZA.
set -euo pipefail
cd /global/scratch/users/tbellg/kmate
PY=/global/home/users/tbellg/miniforge3/envs/basic/bin/python   # statsmodels lives in basic
PR=analysis/grenenet_gea/phase1_replication
OUT=analysis/grenenet_gea/gea_newpanel/results/binomial
echo "== host $(hostname) | binomial SNP vs non-SNP gen9 bio1 | cpus=${SLURM_CPUS_PER_TASK:-8} =="
for CLS in nonsnp snp; do   # nonsnp first (smaller/faster), snp second
    $PY -u $PR/run_binomial.py --class "$CLS" --gen 9 --climate bio1 \
        --out "$OUT" --threads "${SLURM_CPUS_PER_TASK:-8}"
done
echo "== ALLDONE_BINOM snp+nonsnp gen9 bio1 =="
