#!/bin/bash
#SBATCH --job-name=site_bk
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=4:00:00
#SBATCH --output=analysis/grenenet_gea/gea_newpanel/site_bk_%j.out
#SBATCH --error=analysis/grenenet_gea/gea_newpanel/site_bk_%j.out
# Site-level (honest 31-unit) binomial + Kendall, bio1, gen9, SNP vs non-SNP.
set -euo pipefail
cd /global/scratch/users/tbellg/kmate
PY=/global/home/users/tbellg/miniforge3/envs/basic/bin/python   # statsmodels
G=analysis/grenenet_gea/gea_newpanel
echo "== host $(hostname) | site-level binomial+kendall gen9 bio1 =="
for CLS in nonsnp snp; do
    $PY -u $G/run_site_binom_kendall.py --class "$CLS" --gen 9 --climate bio1 \
        --threads "${SLURM_CPUS_PER_TASK:-8}"
done
echo "== ALLDONE_SITE_BK =="
