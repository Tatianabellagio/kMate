#!/bin/bash
#SBATCH --job-name=p231_a3_cn
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=8:00:00
#SBATCH --output=logs/03_cn_full_%j.out
#SBATCH --error=logs/03_cn_full_%j.err

# =============================================================================
# control_p231 Phase A3 -- cn_full k-mer index from the ARCH3 canonical VCF
# (merged_231_chr1_final.vcf.gz), for single-source provenance.
#
# Identical recipe to poolfreq/scripts/build_cn_full_v3qc_v3_chr1.sh EXCEPT
# --vcf points at arch3/chr1/merged_231_chr1_final.vcf.gz instead of
# pangenie_genotyping/data/v3qc_v3/founders_231_v3qc_v3.haploid.vcf.gz.
# Same pang_135 k-mer dictionary (135-asm graph, matches merged_231's annotation
# topology), same ref, same --treat-missing-as-n. cn_full is consensus-derived so
# this should reproduce the production cn_full nearly exactly (validated by
# 03c_compare). bubble_id comes from the pang_135 dictionary (for ω=1/m_b).
# =============================================================================
mkdir -p logs
set -euo pipefail
BASE=/global/scratch/users/tbellg/hapfire_sv
PY=/global/home/users/tbellg/miniforge3/envs/hapfm/bin/python
CHR=Chr1
KMERS=$BASE/pangenie_genotyping/data/pang_135_pangenie_index_${CHR}_kmers.tsv.gz
VCF=$BASE/arch3/chr1/merged_231_chr1_final.vcf.gz
REF=/global/scratch/users/tbellg/pang/pang_1001gplus/20260209_Exposito-Alonso/chr_only/TAIR10.chr.iupacN.fa
OUT_DIR=$BASE/control_p231/data/cn_full_p231
mkdir -p $OUT_DIR
OUT_PREFIX=$OUT_DIR/cn_${CHR}

[ -s "$VCF" ]   || { echo "ERROR: missing $VCF"; exit 1; }
[ -s "$KMERS" ] || { echo "ERROR: missing $KMERS"; exit 1; }
[ ! -s "${OUT_PREFIX}.cn.npz" ] || { echo "exists, skipping"; exit 0; }

echo "[$(date)] $CHR cn_full p231 from arch3 merged_231 (--treat-missing-as-n)"
$PY -u $BASE/poolfreq/src/build_kmer_cn.py \
    --kmers "$KMERS" --vcf "$VCF" --ref "$REF" \
    --chrom "$CHR" --out "$OUT_PREFIX" \
    --treat-missing-as-n
echo "[$(date)] $CHR DONE"
ls -lh ${OUT_PREFIX}.cn.npz ${OUT_PREFIX}.meta.npz
