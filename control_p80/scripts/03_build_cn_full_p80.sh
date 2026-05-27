#!/bin/bash
#SBATCH --job-name=p80_a3_cnfull
#SBATCH --partition=bse
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=8:00:00
#SBATCH --output=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/control_p80/logs/03_cn_full_%j.out
#SBATCH --error=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/control_p80/logs/03_cn_full_%j.err

# =============================================================================
# Phase A3 -- Build cn_full_p80/cn_Chr1.{cn,meta}.npz
# Path B (on-the-fly per-bubble reconstruction; matches v3qc cn_full build).
#
# K-mer index: pang_135 PG-index (production, built on the 135-assembly graph).
# Matches the graph used by the arch3 A1 catalog that produced our canonical
# biallelic VCF. Using a different graph's PG-index (e.g. control_p82's 82-acc
# index) would mismatch bubble topology and silently produce wrong cn_full
# reconstructions.
# =============================================================================
set -euo pipefail

CTRL=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/control_p80
BASE=/carnegie/nobackup/scratch/tbellagio/hapfire_sv
PY=/home/tbellagio/miniforge3/envs/hapfm/bin/python

KMERS=$BASE/pangenie_genotyping/data/pang_135_pangenie_index_Chr1_kmers.tsv.gz
VCF=$CTRL/data/pangenome_p80_chr1.vcf.gz
REF=/home/tbellagio/scratch/pang/pang_1001gplus/20260209_Exposito-Alonso/chr_only/TAIR10.chr.iupacN.fa
OUT_DIR=$CTRL/data/cn_full_p80
OUT_PREFIX=$OUT_DIR/cn_Chr1
mkdir -p $OUT_DIR

for f in "$KMERS" "$VCF" "$REF"; do
    [ -s "$f" ] || { echo "ERROR: missing $f" >&2; exit 1; }
done

if [ -s "${OUT_PREFIX}.cn.npz" ] && [ -s "${OUT_PREFIX}.meta.npz" ]; then
    echo "[$(date)] cn_full_p80 already present -- skip"
    ls -lh ${OUT_PREFIX}.cn.npz ${OUT_PREFIX}.meta.npz
    exit 0
fi

echo "[$(date)] build cn_full_p80 (path B, on-the-fly reconstruction)"
echo "  kmers: $KMERS"
echo "  vcf:   $VCF"
echo "  ref:   $REF"
echo "  out:   ${OUT_PREFIX}.{cn,meta}.npz"

$PY -u $BASE/poolfreq/src/build_kmer_cn.py \
    --kmers "$KMERS" \
    --vcf   "$VCF" \
    --ref   "$REF" \
    --chrom "Chr1" \
    --out   "$OUT_PREFIX"

echo ""
echo "[$(date)] DONE"
ls -lh ${OUT_PREFIX}.cn.npz ${OUT_PREFIX}.meta.npz

$PY <<EOF
import scipy.sparse as sp, numpy as np
M = sp.load_npz("${OUT_PREFIX}.cn.npz")
m = np.load("${OUT_PREFIX}.meta.npz", allow_pickle=True)
density = M.nnz / (M.shape[0] * M.shape[1])
print(f"sanity: shape={M.shape}  nnz={M.nnz:,}  density={density*100:.3f}%")
print(f"founders: {len(m['founders'])} -- first 5: {list(m['founders'][:5])}")
if density < 0.005:
    print("  WARNING: density below 0.5% -- investigate")
EOF
