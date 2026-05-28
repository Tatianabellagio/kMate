#!/bin/bash
#SBATCH --job-name=rn_v3qc_v2
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --time=1:00:00
#SBATCH --output=logs/rn_v3qc_v2_%A_%a.out
#SBATCH --error=logs/rn_v3qc_v2_%A_%a.err

# Row-normalize cn_full_231_v3qc_v2 → cn_full_231_v3qc_v2_rownorm (per-chrom array)
# Each founder's row sums to 1.0 → equal "evidence budget" per founder.
# Chr1 already built by build_cn_full_rownorm_v3qc_v2_chr1.sh — this generalizes.
mkdir -p logs
set -euo pipefail
PY=/global/home/users/tbellg/miniforge3/envs/hapfm/bin/python
CHR="Chr${SLURM_ARRAY_TASK_ID:-2}"

SRC=/global/scratch/users/tbellg/kmate/poolfreq/data/cn_full_231_v3qc_v2
OUT=/global/scratch/users/tbellg/kmate/poolfreq/data/cn_full_231_v3qc_v2_rownorm
mkdir -p $OUT

[ -s "$SRC/cn_${CHR}.cn.npz" ] || { echo "ERROR: missing $SRC/cn_${CHR}.cn.npz"; exit 1; }
[ ! -s "$OUT/cn_${CHR}.cn.npz" ] || { echo "$CHR exists, skipping"; exit 0; }

echo "[$(date)] $CHR: row-norm cn_full_v3qc_v2"
$PY << EOF
import numpy as np
from scipy.sparse import load_npz, save_npz, diags
from pathlib import Path
import shutil

SRC = Path("${SRC}")
OUT = Path("${OUT}")
chrom = "${CHR}"

cn = load_npz(SRC / f"cn_{chrom}.cn.npz").tocsr()
F, K = cn.shape
print(f"  shape: ({F}, {K:,})  nnz: {cn.nnz:,}")
Kf = np.asarray(cn.sum(axis=1)).flatten().astype(np.float32)
print(f"  K_f range: [{int(Kf.min()):,}, {int(Kf.max()):,}]")
inv_Kf = np.divide(1.0, Kf, where=(Kf != 0), out=np.zeros_like(Kf, dtype=np.float32))
D = diags(inv_Kf)
cn_norm = (D @ cn).astype(np.float32).tocsr()
row_sums = np.asarray(cn_norm.sum(axis=1)).flatten()
print(f"  row sums after norm: [{row_sums.min():.6f}, {row_sums.max():.6f}]")
save_npz(OUT / f"cn_{chrom}.cn.npz", cn_norm)
shutil.copy(SRC / f"cn_{chrom}.meta.npz", OUT / f"cn_{chrom}.meta.npz")
print(f"  saved {OUT}/cn_{chrom}.{{cn,meta}}.npz")
EOF
echo "[$(date)] $CHR: DONE"
ls -lh $OUT/cn_${CHR}.cn.npz $OUT/cn_${CHR}.meta.npz
