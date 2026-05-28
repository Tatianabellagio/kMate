#!/bin/bash
#SBATCH --job-name=build_dose
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=4:00:00
#SBATCH --output=logs/build_dose_%j.out
#SBATCH --error=logs/build_dose_%j.err

mkdir -p logs
set -uo pipefail
cd /global/scratch/users/tbellg/kmate

PYTHON=/global/home/users/tbellg/miniforge3/envs/hapfm/bin/python
OUT=data/block_haplotype_cn/chr1_full_clean_dose.npz

echo "[$(date)] starting dose-aware clean reconstruction"
$PYTHON -u scripts/build_block_haplotype_cn_clean.py \
    --block-index sims/visor_freqk/chr1_only_panel/hapfire_block_index_chr1.npz \
    --cn-var data/cn_var_231_v2.cn_var.npz \
    --cn-var-meta data/cn_var_231_v2.meta.npz \
    --vcf imputation/work/test_fix/merged_v2_chr.vcf.gz \
    --ref /global/scratch/users/tbellg/pang/pang_1001gplus/20260209_Exposito-Alonso/chr_only/TAIR10.chr.fa \
    --out $OUT \
    --chrom-filter Chr1 --threads 8 \
    --dose-aware

echo "[$(date)] DONE → $OUT"
ls -lh $OUT
