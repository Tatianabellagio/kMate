#!/bin/bash
#SBATCH --job-name=filt2_v3qc
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --time=1:00:00
#SBATCH --output=logs/filt2_v3qc_%A_%a.out
#SBATCH --error=logs/filt2_v3qc_%A_%a.err

# Drop k-mer columns where ac_k < 2 from cn_full_231_v3qc. One chrom per task.
mkdir -p logs
set -uo pipefail
PY=/global/home/users/tbellg/miniforge3/envs/hapfm/bin/python
CHR=${SLURM_ARRAY_TASK_ID:-1}

$PY << EOF
import numpy as np
from scipy.sparse import load_npz, save_npz
from pathlib import Path

chrom = f"Chr${CHR}"
SRC = Path('/global/scratch/users/tbellg/kmate/poolfreq/data/cn_full_231_v3qc')
OUT = Path('/global/scratch/users/tbellg/kmate/poolfreq/data/cn_full_231_v3qc_filt2')
OUT.mkdir(exist_ok=True)

cn = load_npz(SRC / f'cn_{chrom}.cn.npz')
meta = np.load(SRC / f'cn_{chrom}.meta.npz', allow_pickle=True)
F, K = cn.shape
print(f'[{chrom}] cn shape ({F}, {K:,}) nnz={cn.nnz:,}', flush=True)

ac = np.asarray(cn.sum(axis=0)).flatten().astype(np.int32)
keep = ac >= 2
print(f'[{chrom}] keep ac>=2: {keep.sum():,}/{K:,} (ac=0: {(ac==0).sum():,}, ac=1: {(ac==1).sum():,})', flush=True)

cn_f = cn.tocsc()[:, keep].tocsr()
new_meta = {}
for k in meta.keys():
    a = meta[k]
    if a.shape and a.shape[0] == K:
        new_meta[k] = a[keep]
    else:
        new_meta[k] = a

save_npz(OUT / f'cn_{chrom}.cn.npz', cn_f)
np.savez(OUT / f'cn_{chrom}.meta.npz', **new_meta)
print(f'[{chrom}] DONE filtered nnz={cn_f.nnz:,} shape={cn_f.shape}', flush=True)
EOF
