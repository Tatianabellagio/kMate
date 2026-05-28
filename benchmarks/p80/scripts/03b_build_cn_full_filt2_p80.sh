#!/bin/bash
#SBATCH --job-name=p80_a3b_filt2
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=1:00:00
#SBATCH --output=logs/03b_filt2_%j.out
#SBATCH --error=logs/03b_filt2_%j.err

# =============================================================================
# Phase A3b -- Build cn_full_p80_filt2 by post-filtering cn_full_p80.
# Drops k-mer columns where ac_k = cn.sum(axis=0) < 2 (singletons).
# meta arrays whose first dim equals K are subset by the same mask.
# Idempotent: skips if outputs already exist.
# =============================================================================
mkdir -p logs
set -uo pipefail

PY=/global/home/users/tbellg/miniforge3/envs/hapfm/bin/python
CTRL=/global/scratch/users/tbellg/kmate/benchmarks/p80
SRC=$CTRL/data/cn_full_p80
OUT=$CTRL/data/cn_full_p80_filt2
mkdir -p $OUT

if [ -s "$OUT/cn_Chr1.cn.npz" ] && [ -s "$OUT/cn_Chr1.meta.npz" ]; then
    echo "[$(date)] cn_full_p80_filt2 already present -- skip"
    ls -lh $OUT/cn_Chr1.cn.npz $OUT/cn_Chr1.meta.npz
    exit 0
fi

echo "[$(date)] filt2 (ac>=2) build from $SRC -> $OUT"

$PY << EOF
import numpy as np
from scipy.sparse import load_npz, save_npz
from pathlib import Path

chrom = "Chr1"
SRC = Path("${SRC}")
OUT = Path("${OUT}")
OUT.mkdir(exist_ok=True)

cn = load_npz(SRC / f"cn_{chrom}.cn.npz")
meta = np.load(SRC / f"cn_{chrom}.meta.npz", allow_pickle=True)
F, K = cn.shape
print(f"[{chrom}] cn shape ({F}, {K:,}) nnz={cn.nnz:,}", flush=True)

ac = np.asarray(cn.sum(axis=0)).flatten().astype(np.int32)
keep = ac >= 2
print(f"[{chrom}] keep ac>=2: {keep.sum():,}/{K:,} (ac=0: {(ac==0).sum():,}, ac=1: {(ac==1).sum():,})", flush=True)

cn_f = cn.tocsc()[:, keep].tocsr()
new_meta = {}
for k in meta.keys():
    a = meta[k]
    if a.shape and a.shape[0] == K:
        new_meta[k] = a[keep]
    else:
        new_meta[k] = a

save_npz(OUT / f"cn_{chrom}.cn.npz", cn_f)
np.savez(OUT / f"cn_{chrom}.meta.npz", **new_meta)
print(f"[{chrom}] DONE filtered nnz={cn_f.nnz:,} shape={cn_f.shape}", flush=True)
EOF

echo ""
echo "[$(date)] DONE"
ls -lh $OUT/cn_Chr1.cn.npz $OUT/cn_Chr1.meta.npz
