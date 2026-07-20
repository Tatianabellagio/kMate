#!/bin/bash
#SBATCH --job-name=kendall_snpnon
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=2:00:00
#SBATCH --output=analysis/grenenet_gea/gea_newpanel/kendall_snpnon_%j.out
#SBATCH --error=analysis/grenenet_gea/gea_newpanel/kendall_snpnon_%j.out
# Raw Kendall-tau climate-GEA (bio1, gen9), SNP vs non-SNP, new-panel comparison.
# Per-record tau-b(AF, bio1) across the 355 gen9 site_gen_plot pools. No WZA.
set -euo pipefail
cd /global/scratch/users/tbellg/kmate
PY=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python
PR=analysis/grenenet_gea/phase1_replication
OUT=analysis/grenenet_gea/gea_newpanel/results/kendall
echo "== host $(hostname) | kendall SNP vs non-SNP gen9 bio1 =="
for CLS in snp nonsnp; do
    $PY -u $PR/run_kendall.py --class "$CLS" --gen 9 --climate bio1 \
        --out "$OUT" --threads "${SLURM_CPUS_PER_TASK:-8}"
done
echo "== ALLDONE_KENDALL snp+nonsnp gen9 bio1 =="
