#!/bin/bash
#SBATCH --job-name=build_clean
#SBATCH --partition=bse
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=4:00:00
#SBATCH --output=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/poolfreq/data/block_haplotype_cn/build_clean_%j.out
#SBATCH --error=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/poolfreq/data/block_haplotype_cn/build_clean_%j.err

set -uo pipefail
cd /carnegie/nobackup/scratch/tbellagio/hapfire_sv

PYTHON=/home/tbellagio/miniforge3/envs/hapfm/bin/python
MAX_BLOCKS=${MAX_BLOCKS:-0}
EXTRA=""
if [[ "$MAX_BLOCKS" -gt 0 ]]; then
    EXTRA="--max-blocks $MAX_BLOCKS"
    OUT=poolfreq/data/block_haplotype_cn/chr1_first${MAX_BLOCKS}_clean.npz
else
    OUT=poolfreq/data/block_haplotype_cn/chr1_full_clean.npz
fi

echo "[$(date)] starting clean reconstruction (max_blocks=$MAX_BLOCKS)"
$PYTHON -u poolfreq/scripts/build_block_haplotype_cn_clean.py \
    --block-index sims/visor_freqk/chr1_only_panel/hapfire_block_index_chr1.npz \
    --cn-var poolfreq/data/cn_var_231_v2.cn_var.npz \
    --cn-var-meta poolfreq/data/cn_var_231_v2.meta.npz \
    --vcf imputation/work/test_fix/merged_v2_chr.vcf.gz \
    --ref /home/tbellagio/scratch/pang/pang_1001gplus/20260209_Exposito-Alonso/chr_only/TAIR10.chr.fa \
    --out $OUT \
    --chrom-filter Chr1 --threads 8 \
    $EXTRA

echo "[$(date)] DONE → $OUT"
ls -lh $OUT
