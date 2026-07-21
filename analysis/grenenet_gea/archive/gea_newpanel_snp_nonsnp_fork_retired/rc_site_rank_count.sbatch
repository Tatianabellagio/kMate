#!/bin/bash
#SBATCH --job-name=rc_bk
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=6:00:00
#SBATCH --array=0-17
#SBATCH --output=analysis/grenenet_gea/gea_newpanel/rc_%A_%a.out
#SBATCH --error=analysis/grenenet_gea/gea_newpanel/rc_%A_%a.out
# Site-level (31-unit) Kendall + binomial GEA, all class x axis combos.
set -euo pipefail
cd /global/scratch/users/tbellg/kmate
PY=/global/home/users/tbellg/miniforge3/envs/basic/bin/python   # statsmodels
G=analysis/grenenet_gea/gea_newpanel

CLASSES=(snp nonsnp sv)
AXES=(bio5 pc1 bio12 bio13 bio16 bio19)
i=${SLURM_ARRAY_TASK_ID}
CLS=${CLASSES[$((i / 6))]}
AXIS=${AXES[$((i % 6))]}

echo "== host $(hostname) | task $i -> class=$CLS axis=$AXIS =="
$PY -u $G/rc_site_rank_count.py --class "$CLS" --axis "$AXIS" \
    --threads "${SLURM_CPUS_PER_TASK:-8}"
echo "== ALLDONE task $i $CLS $AXIS =="
