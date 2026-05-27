#!/bin/bash
#SBATCH --job-name=cn_231
#SBATCH --partition=bse
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=12:00:00
#SBATCH --output=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/imputation/logs/cn_231_%j.out
#SBATCH --error=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/imputation/logs/cn_231_%j.err

# Post-imputation: build cn_kmer_231 and cn_var_231 by merging the 80-founder
# cactus VCF with the 151-founder imputed VCF, renaming chromosomes back to
# Chr1/2/3/4/5 (cactus convention), then running build_kmer_cn.py and
# build_cn_var.py.

set -euo pipefail
cd /carnegie/nobackup/scratch/tbellagio/hapfire_sv/imputation/work
mkdir -p cn_231

BCF=/home/tbellagio/miniforge3/envs/sequencing_pipeline/bin/bcftools
TABIX=/home/tbellagio/miniforge3/envs/sequencing_pipeline/bin/tabix
PYTHON=/home/tbellagio/miniforge3/envs/hapfm/bin/python

# 1. Merge: 80-founder cactus SVs + 151-founder imputed VCF → 231-founder unified VCF
echo "[$(date)] Step 1: merge cactus_svs_renamed (80) + imputed_151 (151) → 231-founder VCF"
$BCF merge cactus_svs_renamed.vcf.gz imputed_151.vcf.gz \
    -Oz -o cn_231/merged_231.vcf.gz --threads 8
$TABIX -p vcf cn_231/merged_231.vcf.gz

# 2. Rename chroms 1→Chr1 etc. (cactus convention used by build_kmer_cn.py)
echo "[$(date)] Step 2: rename chroms 1→Chr1"
echo -e "1\tChr1\n2\tChr2\n3\tChr3\n4\tChr4\n5\tChr5" > cn_231/chrom_back.txt
$BCF annotate --rename-chrs cn_231/chrom_back.txt cn_231/merged_231.vcf.gz \
    -Oz -o cn_231/merged_231_chr.vcf.gz --threads 8
$TABIX -p vcf cn_231/merged_231_chr.vcf.gz

# 3. Build cn_kmer_231 per chromosome (using existing build_kmer_cn.py)
echo "[$(date)] Step 3: build cn_kmer_231 per chromosome"
mkdir -p /carnegie/nobackup/scratch/tbellagio/hapfire_sv/poolfreq/data/cn_full_231
for chrom in Chr1 Chr2 Chr3 Chr4 Chr5; do
    echo "[$(date)]   $chrom"
    $PYTHON /carnegie/nobackup/scratch/tbellagio/hapfire_sv/poolfreq/src/build_kmer_cn.py \
        --kmers /carnegie/nobackup/scratch/tbellagio/hapfire_sv/pangenie_test/pangenie_idx_rawv2_${chrom}_kmers.tsv.gz \
        --vcf   cn_231/merged_231_chr.vcf.gz \
        --ref   /home/tbellagio/scratch/pang/pang_1001gplus/20260209_Exposito-Alonso/chr_only/TAIR10.chr.fa \
        --chrom $chrom \
        --out   /carnegie/nobackup/scratch/tbellagio/hapfire_sv/poolfreq/data/cn_full_231/cn_${chrom}
done

# 4. Build cn_var_231 (founder × biallelic-VCF-record)
echo "[$(date)] Step 4: build cn_var_231"
$PYTHON /carnegie/nobackup/scratch/tbellagio/hapfire_sv/poolfreq/src/build_cn_var.py \
    --vcf cn_231/merged_231_chr.vcf.gz \
    --out /carnegie/nobackup/scratch/tbellagio/hapfire_sv/poolfreq/data/cn_var_231

echo "[$(date)] DONE — 231-founder cn matrices built"
ls -la /carnegie/nobackup/scratch/tbellagio/hapfire_sv/poolfreq/data/cn_full_231/
ls -la /carnegie/nobackup/scratch/tbellagio/hapfire_sv/poolfreq/data/cn_var_231*
