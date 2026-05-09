#!/bin/bash
#SBATCH --job-name=cn_full_v3
#SBATCH --partition=bse
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=8:00:00
#SBATCH --output=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/poolfreq/logs/cn_full_v3_%A_%a.out
#SBATCH --error=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/poolfreq/logs/cn_full_v3_%A_%a.err

# =============================================================================
# build_cn_full_v3_one.sh
#
# Build cn_full_231_v3 for one chrom (selected by SLURM_ARRAY_TASK_ID, 1..5).
#
# Uses path B: v2's on-the-fly per-bubble reconstruction
# (poolfreq/src/build_kmer_cn.py) with the new haploid VCF as input.
# Path A (FASTA-keyed lookup) was attempted first and produced a 250×
# density drop because bcftools consensus output uses founder-local
# coordinates; lookups at TAIR10 ref coords land on the wrong sequence
# once cumulative indel offsets accumulate (Chr1 example: 100001's
# consensus is 424 kb shorter than TAIR10).
#
# v3 vs v2: same builder, same TAIR10 ref, same pang_135 panel index;
# different VCF (haploid unimputed vs Beagle-imputed).
# =============================================================================
set -euo pipefail

BASE=/carnegie/nobackup/scratch/tbellagio/hapfire_sv
PY=/home/tbellagio/miniforge3/envs/hapfm/bin/python

T=${SLURM_ARRAY_TASK_ID:-1}
CHR="Chr${T}"

KMERS=$BASE/pangenie_genotyping/data/pang_135_pangenie_index_${CHR}_kmers.tsv.gz
VCF=$BASE/pangenie_genotyping/data/merged/founders_231_chr.haploid.vcf.gz
REF=/home/tbellagio/scratch/pang/pang_1001gplus/20260209_Exposito-Alonso/chr_only/TAIR10.chr.iupacN.fa
OUT_DIR=$BASE/poolfreq/data/cn_full_231_v3
OUT_PREFIX=$OUT_DIR/cn_${CHR}
mkdir -p $OUT_DIR $BASE/poolfreq/logs

[ -s "$KMERS" ] || { echo "ERROR: missing $KMERS" >&2; exit 1; }
[ -s "$VCF"   ] || { echo "ERROR: missing $VCF"   >&2; exit 1; }
[ -s "$REF"   ] || { echo "ERROR: missing $REF"   >&2; exit 1; }
[ ! -s "${OUT_PREFIX}.cn.npz" ] || {
    echo "[$(date)] $CHR: cn_full_v3 already exists at ${OUT_PREFIX}.cn.npz — skipping"
    exit 0
}

echo "[$(date)] $CHR: build cn_full_v3 (path B, on-the-fly reconstruction)"
echo "  kmers: $KMERS"
echo "  vcf:   $VCF"
echo "  ref:   $REF"
echo "  out:   ${OUT_PREFIX}.{cn,meta}.npz"

$PY -u $BASE/poolfreq/src/build_kmer_cn.py \
    --kmers "$KMERS" \
    --vcf   "$VCF" \
    --ref   "$REF" \
    --chrom "$CHR" \
    --out   "$OUT_PREFIX"

echo ""
echo "[$(date)] $CHR: DONE"
ls -lh ${OUT_PREFIX}.cn.npz ${OUT_PREFIX}.meta.npz

# Sanity smoke: density should be in the 1-15% range for cactus_em panel
$PY <<EOF
import scipy.sparse as sp
M = sp.load_npz("${OUT_PREFIX}.cn.npz")
density = M.nnz / (M.shape[0] * M.shape[1])
print(f"  sanity: shape={M.shape}  nnz={M.nnz:,}  density={density*100:.3f}%")
if density < 0.005:
    print("  ⚠️  density below 0.5% — possible builder regression; investigate")
EOF
