#!/bin/bash
#SBATCH --job-name=cn_pg
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=12:00:00
#SBATCH --output=logs/cn_pg_%j.out
#SBATCH --error=logs/cn_pg_%j.err

# =============================================================================
# build_cn_from_pangenie.sh
# Stage 6: build cn_kmer + cn_var matrices from the PanGenie-derived
# 231-founder catalog (output of merge_vcfs.sh).
#
# Inputs:
#   1. founders_231_chr.vcf.gz  ← from merge_vcfs.sh (Stage 4)
#   2. pang69_pangenie_index_<chrom>_kmers.tsv.gz   ← from build_pangenie_index.sh
#   3. TAIR10.chr.fa            ← cactus reference
#
# Outputs (under pangenie_genotyping/data/cn/):
#   - cn_kmer_<chrom>.cn.npz, .meta.npz   (one per chrom; founder × k-mer)
#   - cn_var_231.cn.npz, .meta.npz        (founder × biallelic-VCF-record)
#
# Mirrors imputation/06_build_231_cn_matrices.sh but plugged into the
# PanGenie-genotyping panel (no Beagle imputation in the chain).
# =============================================================================
mkdir -p logs
set -euo pipefail

BASE=/global/scratch/users/tbellg/kmate/pangenie_genotyping
MERGED_VCF=$BASE/data/merged/founders_231_chr.vcf.gz
INDEX_PREFIX=$BASE/data/pang69_pangenie_index   # produced by build_pangenie_index.sh
REF=/global/scratch/users/tbellg/pang/pang_1001gplus/20260209_Exposito-Alonso/chr_only/TAIR10.chr.fa

OUT_DIR=$BASE/data/cn
mkdir -p $OUT_DIR

BCF=/global/home/users/tbellg/miniforge3/envs/sequencing_pipeline/bin/bcftools
TABIX=/global/home/users/tbellg/miniforge3/envs/sequencing_pipeline/bin/tabix
PYTHON=/global/home/users/tbellg/miniforge3/envs/hapfm/bin/python

POOLFREQ_SRC=/global/scratch/users/tbellg/kmate/src

# ---- Sanity checks on inputs --------------------------------------------------
for f in $MERGED_VCF $REF; do
    [ -s "$f" ] || { echo "ERROR: missing input $f" >&2; exit 1; }
done
for chrom in Chr1 Chr2 Chr3 Chr4 Chr5; do
    KFILE=${INDEX_PREFIX}_${chrom}_kmers.tsv.gz
    [ -s "$KFILE" ] || { echo "ERROR: missing PanGenie kmers $KFILE — run build_pangenie_index.sh" >&2; exit 1; }
done

# ---- Step 1: biallelic-normalize the merged 231-founder VCF ------------------
# build_cn_var.py expects one ALT per record; build_kmer_cn.py also walks
# bubble-by-bubble and works correctly when multi-ALT bubbles are split.
NORM_VCF=$OUT_DIR/founders_231_chr.norm.vcf.gz
if [ ! -s $NORM_VCF ]; then
    echo "[$(date)] Step 1: bcftools norm -m- (split multi-allelic into biallelic)"
    $BCF norm -m- -f $REF $MERGED_VCF -Oz -o $NORM_VCF --threads 8
    $TABIX -p vcf $NORM_VCF
fi
echo "  norm records: $($BCF view -H $NORM_VCF | wc -l)"
echo "  norm samples: $($BCF query -l $NORM_VCF | wc -l)"

# ---- Step 2: build cn_kmer per chromosome ------------------------------------
echo "[$(date)] Step 2: build cn_kmer_<chrom> using PanGenie index k-mers"
for chrom in Chr1 Chr2 Chr3 Chr4 Chr5; do
    KFILE=${INDEX_PREFIX}_${chrom}_kmers.tsv.gz
    OUT_PREFIX=$OUT_DIR/cn_kmer_${chrom}
    if [ -s ${OUT_PREFIX}.cn.npz ]; then
        echo "  [skip] $chrom (cn_kmer already exists)"
        continue
    fi
    echo "  [run] $chrom"
    $PYTHON $POOLFREQ_SRC/build_kmer_cn.py \
        --kmers $KFILE \
        --vcf   $NORM_VCF \
        --ref   $REF \
        --chrom $chrom \
        --out   $OUT_PREFIX
done

# ---- Step 3: build cn_var (founder × biallelic-VCF-record) -------------------
echo "[$(date)] Step 3: build cn_var_231"
CN_VAR_PREFIX=$OUT_DIR/cn_var_231
if [ ! -s ${CN_VAR_PREFIX}.cn.npz ]; then
    $PYTHON $POOLFREQ_SRC/build_cn_var.py \
        --vcf $NORM_VCF \
        --out $CN_VAR_PREFIX
fi

echo "[$(date)] DONE — cn matrices in $OUT_DIR"
ls -lh $OUT_DIR/cn_kmer_Chr*.cn.npz $OUT_DIR/cn_var_231*.npz 2>/dev/null
