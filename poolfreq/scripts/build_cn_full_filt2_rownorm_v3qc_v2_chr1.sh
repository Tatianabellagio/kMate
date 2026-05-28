#!/bin/bash
#SBATCH --job-name=f2rn_chr1
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --time=1:00:00
#SBATCH --output=logs/f2rn_chr1_%j.out
#SBATCH --error=logs/f2rn_chr1_%j.err

# filt2 + row-norm on cn_full_v3qc_v2 Chr1.
# Order: filt2 (drop ac_k<2 columns) FIRST, then row-norm (so post-filter rows sum to 1).
# This kills per-founder private-kmer uniqueness asymmetry AND total-budget asymmetry.
mkdir -p logs
set -euo pipefail
PY=/global/home/users/tbellg/miniforge3/envs/hapfm/bin/python

$PY << 'EOF'
import numpy as np
from scipy.sparse import load_npz, save_npz, diags
from pathlib import Path
import shutil

SRC = Path('/global/scratch/users/tbellg/kmate/poolfreq/data/cn_full_231_v3qc_v2')
OUT = Path('/global/scratch/users/tbellg/kmate/poolfreq/data/cn_full_231_v3qc_v2_filt2_rownorm')
OUT.mkdir(exist_ok=True)

chrom = 'Chr1'
cn = load_npz(SRC / f'cn_{chrom}.cn.npz')
meta = np.load(SRC / f'cn_{chrom}.meta.npz', allow_pickle=True)
F, K = cn.shape
print(f'[{chrom}] input shape ({F}, {K:,}) nnz={cn.nnz:,}', flush=True)

# --- Step 1: filt2 ---
ac = np.asarray(cn.sum(axis=0)).flatten().astype(np.int32)
keep = ac >= 2
print(f'[{chrom}] filt2 keep ac>=2: {keep.sum():,}/{K:,}', flush=True)
print(f'           ac=0: {(ac==0).sum():,}, ac=1: {(ac==1).sum():,}', flush=True)
print(f'           ac=2: {(ac==2).sum():,}, ac=3-10: {((ac>=3)&(ac<=10)).sum():,}', flush=True)

cn_f = cn.tocsc()[:, keep].tocsr()
new_meta = {}
for k in meta.keys():
    a = meta[k]
    if a.shape and a.shape[0] == K:
        new_meta[k] = a[keep]
    else:
        new_meta[k] = a
print(f'[{chrom}] post-filt2 nnz={cn_f.nnz:,}', flush=True)

# --- Step 2: row-norm AFTER filt2 ---
Kf = np.asarray(cn_f.sum(axis=1)).flatten().astype(np.float32)
print(f'[{chrom}] K_f after filt2: min={int(Kf.min()):,}, max={int(Kf.max()):,}, median={int(np.median(Kf)):,}', flush=True)
# Split by side (assume first 78 are cactus, last 153 are PG — verify via meta if needed)
cactus_mask = np.zeros(F, dtype=bool); cactus_mask[:78] = True
print(f'  K_f cactus median={int(np.median(Kf[cactus_mask])):,}, PG median={int(np.median(Kf[~cactus_mask])):,}, ratio={np.median(Kf[cactus_mask])/np.median(Kf[~cactus_mask]):.3f}', flush=True)

inv_Kf = np.divide(1.0, Kf, where=(Kf != 0), out=np.zeros_like(Kf, dtype=np.float32))
D = diags(inv_Kf)
cn_norm = (D @ cn_f).astype(np.float32).tocsr()
row_sums = np.asarray(cn_norm.sum(axis=1)).flatten()
print(f'[{chrom}] row sums after norm: [{row_sums.min():.6f}, {row_sums.max():.6f}]', flush=True)

save_npz(OUT / f'cn_{chrom}.cn.npz', cn_norm)
np.savez(OUT / f'cn_{chrom}.meta.npz', **new_meta)
print(f'[{chrom}] DONE, saved to {OUT}', flush=True)
EOF

echo "[$(date)] DONE"
ls -lh /global/scratch/users/tbellg/kmate/poolfreq/data/cn_full_231_v3qc_v2_filt2_rownorm/
