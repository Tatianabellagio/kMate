#!/bin/bash
#SBATCH --job-name=p80_a4_cnvar
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=2:00:00
#SBATCH --output=logs/04_var_pa_%j.out
#SBATCH --error=logs/04_var_pa_%j.err

# =============================================================================
# Phase A4 -- Build var_pa_p80.{var_pa,meta}.npz
# Wraps src/build_var_pa.py with the canonical p80 VCF.
# =============================================================================
mkdir -p logs
set -euo pipefail

CTRL=/global/scratch/users/tbellg/kmate/benchmarks/p80
BASE=/global/scratch/users/tbellg/kmate
PY=/global/home/users/tbellg/miniforge3/envs/hapfm/bin/python

VCF=$CTRL/data/pangenome_p80_chr1.vcf.gz
OUT_PREFIX=$CTRL/data/var_pa_p80
mkdir -p $CTRL/data

[ -s "$VCF" ] || { echo "ERROR: missing $VCF -- run 01 first" >&2; exit 1; }

if [ -s "${OUT_PREFIX}.var_pa.npz" ] && [ -s "${OUT_PREFIX}.meta.npz" ]; then
    echo "[$(date)] var_pa_p80 already present -- skip"
    ls -lh ${OUT_PREFIX}.var_pa.npz ${OUT_PREFIX}.meta.npz
    exit 0
fi

echo "[$(date)] build var_pa_p80"
echo "  vcf: $VCF"
echo "  out: ${OUT_PREFIX}.{var_pa,meta}.npz"

$PY -u $BASE/src/build_var_pa.py --vcf "$VCF" --out "$OUT_PREFIX"

echo ""
echo "=== sanity check ==="
$PY <<EOF
import numpy as np, scipy.sparse as sp
M = sp.load_npz("${OUT_PREFIX}.var_pa.npz")
m = np.load("${OUT_PREFIX}.meta.npz", allow_pickle=True)
print(f"var_pa: shape={M.shape}  dtype={M.dtype}  nnz={M.nnz:,}  density={M.nnz/(M.shape[0]*M.shape[1])*100:.2f}%")
print(f"founders: {len(m['founders'])} -- first 5: {list(m['founders'][:5])}")
print(f"per-chrom records: {dict(zip(*np.unique(m['chrom'], return_counts=True)))}")
EOF

echo ""
echo "[$(date)] DONE -> ${OUT_PREFIX}.{var_pa,meta}.npz"
ls -lh ${OUT_PREFIX}.var_pa.npz ${OUT_PREFIX}.meta.npz
