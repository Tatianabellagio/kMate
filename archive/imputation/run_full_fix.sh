#!/bin/bash
#SBATCH --job-name=full_fix
#SBATCH --partition=bse
#SBATCH --cpus-per-task=8
#SBATCH --mem=256G
#SBATCH --time=18:00:00
#SBATCH --output=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/imputation/logs/full_fix_%j.out
#SBATCH --error=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/imputation/logs/full_fix_%j.err

# =============================================================================
# run_full_fix.sh
# Tier 1 root-cause fix, full genome (all 5 chromosomes).
# Validated on Chr4 (job 57779): slope 1.43 → 1.003, R² 0.68 → 0.994.
#
# Steps:
#   1. bcftools merge ref_80 (with cactus SNPs!) + imputed_151
#      → merged_231_v2.vcf.gz (replaces old merged_231 that was built from
#      cactus_svs only)
#   2. Rename chrom 1..5 → Chr1..Chr5 (cactus convention for cn build)
#   3. Build cn_kmer_v2 per chromosome (5 × ~1.5h = ~7.5h)
#   4. Build cn_var_v2 (~30 min)
# =============================================================================
set -euo pipefail
cd /carnegie/nobackup/scratch/tbellagio/hapfire_sv/imputation/work

BCF=/home/tbellagio/miniforge3/envs/sequencing_pipeline/bin/bcftools
TABIX=/home/tbellagio/miniforge3/envs/sequencing_pipeline/bin/tabix
PYTHON=/home/tbellagio/miniforge3/envs/hapfm/bin/python
POOLFREQ=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/poolfreq
mkdir -p test_fix

echo "[$(date)] === Step 1: full-genome merge ref_80 + imputed_151 ==="
if [ ! -s test_fix/merged_v2_full.vcf.gz ]; then
    $BCF merge ref_80.vcf.gz imputed_151.vcf.gz \
        -Oz -o test_fix/merged_v2_full.vcf.gz --threads 8
    $TABIX -p vcf test_fix/merged_v2_full.vcf.gz
fi
echo "  records: $($BCF view -H test_fix/merged_v2_full.vcf.gz | wc -l)"
echo "  samples: $($BCF query -l test_fix/merged_v2_full.vcf.gz | wc -l)"

echo ""
echo "[$(date)] === Step 2: rename chroms 1..5 → Chr1..Chr5 ==="
if [ ! -s test_fix/merged_v2_chr.vcf.gz ]; then
    echo -e "1\tChr1\n2\tChr2\n3\tChr3\n4\tChr4\n5\tChr5" > test_fix/rename_full.txt
    $BCF annotate --rename-chrs test_fix/rename_full.txt test_fix/merged_v2_full.vcf.gz \
        -Oz -o test_fix/merged_v2_chr.vcf.gz --threads 8
    $TABIX -p vcf test_fix/merged_v2_chr.vcf.gz
fi

echo ""
echo "[$(date)] === Step 3: build cn_kmer_v2 per chromosome ==="
mkdir -p /carnegie/nobackup/scratch/tbellagio/hapfire_sv/poolfreq/data/cn_full_231_v2
for CHROM in Chr1 Chr2 Chr3 Chr4 Chr5; do
    OUT=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/poolfreq/data/cn_full_231_v2/cn_${CHROM}
    if [ -s ${OUT}.cn.npz ]; then
        echo "  $CHROM already built, skip"; continue
    fi
    echo "[$(date)]   building $CHROM..."
    $PYTHON ${POOLFREQ}/src/build_kmer_cn.py \
        --kmers /carnegie/nobackup/scratch/tbellagio/hapfire_sv/pangenie_test/pangenie_idx_rawv2_${CHROM}_kmers.tsv.gz \
        --vcf   test_fix/merged_v2_chr.vcf.gz \
        --ref   /home/tbellagio/scratch/pang/pang_1001gplus/20260209_Exposito-Alonso/chr_only/TAIR10.chr.fa \
        --chrom $CHROM \
        --out   $OUT
done

echo ""
echo "[$(date)] === Step 4: build cn_var_v2 ==="
$PYTHON ${POOLFREQ}/src/build_cn_var.py \
    --vcf test_fix/merged_v2_chr.vcf.gz \
    --out /carnegie/nobackup/scratch/tbellagio/hapfire_sv/poolfreq/data/cn_var_231_v2

echo ""
echo "[$(date)] DONE — corrected 231-founder cn matrices ready"
ls -la /carnegie/nobackup/scratch/tbellagio/hapfire_sv/poolfreq/data/cn_full_231_v2/ \
       /carnegie/nobackup/scratch/tbellagio/hapfire_sv/poolfreq/data/cn_var_231_v2*
