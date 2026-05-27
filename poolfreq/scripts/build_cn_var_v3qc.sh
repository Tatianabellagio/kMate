#!/bin/bash
#SBATCH --job-name=cn_var_v3qc
#SBATCH --partition=bse
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=2:00:00
#SBATCH --output=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/poolfreq/logs/cn_var_v3qc_%j.out
#SBATCH --error=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/poolfreq/logs/cn_var_v3qc_%j.err

# =============================================================================
# build_cn_var_v3qc.sh
# Build cn_var_231_v3qc from founders_231_v3qc.haploid.vcf.gz.
# Output: poolfreq/data/cn_var_231_v3qc.{cn_var,meta}.npz
# =============================================================================
set -euo pipefail

BASE=/carnegie/nobackup/scratch/tbellagio/hapfire_sv
PY=/home/tbellagio/miniforge3/envs/hapfm/bin/python
mkdir -p $BASE/poolfreq/logs

VCF=$BASE/pangenie_genotyping/data/v3qc/founders_231_v3qc.haploid.vcf.gz
OUT_PREFIX=$BASE/poolfreq/data/cn_var_231_v3qc

[ -s "$VCF" ] || { echo "ERROR: missing $VCF" >&2; exit 1; }
[ ! -s "${OUT_PREFIX}.cn_var.npz" ] || {
    echo "[$(date)] cn_var_v3qc exists — exiting"
    exit 0
}

echo "[$(date)] build cn_var_v3qc"
$PY -u $BASE/poolfreq/src/build_cn_var.py --vcf "$VCF" --out "$OUT_PREFIX"

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

echo "[$(date)] DONE"
ls -lh ${OUT_PREFIX}.cn_var.npz ${OUT_PREFIX}.meta.npz
