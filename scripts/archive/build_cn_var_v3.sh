#!/bin/bash
#SBATCH --job-name=cn_var_v3
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=2:00:00
#SBATCH --output=logs/cn_var_v3_%j.out
#SBATCH --error=logs/cn_var_v3_%j.err

# =============================================================================
# build_cn_var_v3.sh
#
# Build cn_var_231_v3 from the haploid (unimputed) golden-standard VCF.
# Wraps src/build_cn_var.py — the existing builder handles haploid
# GTs correctly (any-non-zero-allele = carrier).
#
# Output:
#   data/cn_var_231_v3.cn_var.npz   (CSR, 231 × N_records int8 binary)
#   data/cn_var_231_v3.meta.npz     (founders, chrom, pos, ref_len, alt_len)
#
# Source VCF is the new haploid sibling of the production catalog:
#   panel/pangenie_genotyping/data/merged/founders_231_chr.haploid.vcf.gz
# (multi-allelics already decomposed; cells uniform {0, 1, .}; no Beagle).
#
# v3 vs v2: same builder, new VCF input.
# =============================================================================
set -euo pipefail

BASE=/global/scratch/users/tbellg/kmate
PY=/global/home/users/tbellg/miniforge3/envs/hapfm/bin/python
mkdir -p $BASE/logs

VCF=$BASE/panel/pangenie_genotyping/data/merged/founders_231_chr.haploid.vcf.gz
OUT_PREFIX=$BASE/data/cn_var_231_v3

[ -s "$VCF" ] || { echo "ERROR: missing $VCF" >&2; exit 1; }
[ ! -s "${OUT_PREFIX}.cn_var.npz" ] || {
    echo "[$(date)] cn_var_v3 already exists at ${OUT_PREFIX}.cn_var.npz — exiting"
    exit 0
}

echo "[$(date)] build cn_var_v3"
echo "  vcf: $VCF"
echo "  out: ${OUT_PREFIX}.{cn_var,meta}.npz"

$PY -u $BASE/src/build_cn_var.py --vcf "$VCF" --out "$OUT_PREFIX"

echo ""
echo "=== sanity check ==="
$PY <<EOF
import numpy as np, scipy.sparse as sp
M = sp.load_npz("${OUT_PREFIX}.cn_var.npz")
m = np.load("${OUT_PREFIX}.meta.npz", allow_pickle=True)
print(f"cn_var: shape={M.shape}  dtype={M.dtype}  nnz={M.nnz:,}  density={M.nnz/(M.shape[0]*M.shape[1])*100:.2f}%")
print(f"founders: {len(m['founders'])}")
print(f"per-chrom records: {dict(zip(*np.unique(m['chrom'], return_counts=True)))}")
EOF

echo ""
echo "[$(date)] DONE -> ${OUT_PREFIX}.{cn_var,meta}.npz"
ls -lh ${OUT_PREFIX}.cn_var.npz ${OUT_PREFIX}.meta.npz
