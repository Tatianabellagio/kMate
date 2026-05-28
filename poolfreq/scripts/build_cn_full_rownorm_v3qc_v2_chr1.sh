#!/bin/bash
#SBATCH --job-name=rownorm_v3qc_v2
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --time=1:00:00
#SBATCH --output=logs/rownorm_v3qc_v2_%j.out
#SBATCH --error=logs/rownorm_v3qc_v2_%j.err

# Build a row-normalized cn_full from cn_full_231_v3qc_v2 — each founder's row
# sums to 1. Saved as float32. Then EM with this cn_full uses equal per-founder
# "evidence budget" so PG founders aren't penalized for lower K_f.
#
# Test on Chr1 only first to see if it rebalances h.
mkdir -p logs
set -euo pipefail
PY=/global/home/users/tbellg/miniforge3/envs/hapfm/bin/python
$PY << 'EOF'
import numpy as np
from scipy.sparse import load_npz, save_npz, csr_matrix
from pathlib import Path

SRC_DIR = Path('/global/scratch/users/tbellg/kmate/poolfreq/data/cn_full_231_v3qc_v2')
OUT_DIR = Path('/global/scratch/users/tbellg/kmate/poolfreq/data/cn_full_231_v3qc_v2_rownorm')
OUT_DIR.mkdir(exist_ok=True)

chrom = 'Chr1'
print(f'[{chrom}] loading cn_full_v3qc_v2 binary matrix')
cn = load_npz(SRC_DIR / f'cn_{chrom}.cn.npz').tocsr()
F, K = cn.shape
print(f'  shape: ({F}, {K:,})  nnz: {cn.nnz:,}')

# Per-founder K_f
Kf = np.asarray(cn.sum(axis=1)).flatten().astype(np.float32)
print(f'  K_f range: [{int(Kf.min()):,}, {int(Kf.max()):,}]')

# Row-normalize: cn_norm[f, k] = cn[f, k] / K_f
# scipy.sparse: use multiply by reciprocal vector
inv_Kf = np.divide(1.0, Kf, where=(Kf != 0), out=np.zeros_like(Kf, dtype=np.float32))
# Multiply each row by 1/K_f (using diagonal matrix multiplication trick)
from scipy.sparse import diags
D = diags(inv_Kf)
cn_norm = (D @ cn).astype(np.float32).tocsr()

# Sanity check
row_sums = np.asarray(cn_norm.sum(axis=1)).flatten()
print(f'  After row-norm: row sums range [{row_sums.min():.6f}, {row_sums.max():.6f}] (should all be ~1.0)')
print(f'  cn_norm dtype: {cn_norm.dtype}, nnz: {cn_norm.nnz:,}')

save_npz(OUT_DIR / f'cn_{chrom}.cn.npz', cn_norm)

# Copy meta unchanged
import shutil
shutil.copy(SRC_DIR / f'cn_{chrom}.meta.npz', OUT_DIR / f'cn_{chrom}.meta.npz')
print(f'  saved to {OUT_DIR}/cn_{chrom}.{{cn,meta}}.npz')
EOF
echo "[$(date)] DONE"
ls -lh /global/scratch/users/tbellg/kmate/poolfreq/data/cn_full_231_v3qc_v2_rownorm/
