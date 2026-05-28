#!/bin/bash
#SBATCH --job-name=p80_a4_cnvar
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=2:00:00
#SBATCH --output=logs/04_cn_var_%j.out
#SBATCH --error=logs/04_cn_var_%j.err

# =============================================================================
# Phase A4 -- Build cn_var_p80.{cn_var,meta}.npz
# Wraps src/build_cn_var.py with the canonical p80 VCF.
# =============================================================================
mkdir -p logs
set -euo pipefail

CTRL=/global/scratch/users/tbellg/kmate/benchmarks/p80
BASE=/global/scratch/users/tbellg/kmate
PY=/global/home/users/tbellg/miniforge3/envs/hapfm/bin/python

VCF=$CTRL/data/pangenome_p80_chr1.vcf.gz
OUT_PREFIX=$CTRL/data/cn_var_p80
mkdir -p $CTRL/data

[ -s "$VCF" ] || { echo "ERROR: missing $VCF -- run 01 first" >&2; exit 1; }

if [ -s "${OUT_PREFIX}.cn_var.npz" ] && [ -s "${OUT_PREFIX}.meta.npz" ]; then
    echo "[$(date)] cn_var_p80 already present -- skip"
    ls -lh ${OUT_PREFIX}.cn_var.npz ${OUT_PREFIX}.meta.npz
    exit 0
fi

echo "[$(date)] build cn_var_p80"
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
print(f"founders: {len(m['founders'])} -- first 5: {list(m['founders'][:5])}")
print(f"per-chrom records: {dict(zip(*np.unique(m['chrom'], return_counts=True)))}")
EOF

echo ""
echo "[$(date)] DONE -> ${OUT_PREFIX}.{cn_var,meta}.npz"
ls -lh ${OUT_PREFIX}.cn_var.npz ${OUT_PREFIX}.meta.npz
