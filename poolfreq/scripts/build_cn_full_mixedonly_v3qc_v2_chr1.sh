#!/bin/bash
#SBATCH --job-name=mixonly
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=2
#SBATCH --mem=48G
#SBATCH --time=1:00:00
#SBATCH --output=logs/mixonly_%A_%a.out
#SBATCH --error=logs/mixonly_%A_%a.err

# Build cn_full subsets filtered to "both-sides-represented" k-mers.
# Array task -> threshold variant:
#   0: mixed-loose   ac_cactus>=1 AND ac_pg>=1 AND ac>=2 (~4.6M Chr1 k-mers)
#   1: mixed-strict  ac_cactus>=5 AND ac_pg>=5            (~1.9M)
#   2: mixed-conserv ac_cactus>=10 AND ac_pg>=10          (~1.6M)
mkdir -p logs
set -euo pipefail
PY=/global/home/users/tbellg/miniforge3/envs/hapfm/bin/python

case $SLURM_ARRAY_TASK_ID in
  0) TAG=mixedloose;   THRESH_C=1;  THRESH_P=1  ;;
  1) TAG=mixedstrict;  THRESH_C=5;  THRESH_P=5  ;;
  2) TAG=mixedconserv; THRESH_C=10; THRESH_P=10 ;;
  *) echo unknown; exit 1 ;;
esac

OUT_DIR=/global/scratch/users/tbellg/hapfire_sv/poolfreq/data/cn_full_231_v3qc_v2_${TAG}
mkdir -p $OUT_DIR

$PY << EOF
import numpy as np, json
from scipy.sparse import load_npz, save_npz
from pathlib import Path

CHR = 'Chr1'
SRC = Path('/global/scratch/users/tbellg/hapfire_sv/poolfreq/data/cn_full_231_v3qc_v2')
OUT = Path('${OUT_DIR}')

cn = load_npz(SRC / f'cn_{CHR}.cn.npz').tocsr()
meta = np.load(SRC / f'cn_{CHR}.meta.npz', allow_pickle=True)
founders = np.asarray(meta['founders']).astype(str)
F, K = cn.shape
print(f'cn_full input: ({F}, {K:,}) nnz={cn.nnz:,}', flush=True)

with open('/global/scratch/users/tbellg/hapfire_sv/data/founder_split_cactus_pg.json') as fp:
    split = json.load(fp)
cactus_set = set(map(str, split['cactus']))
pg_set = set(map(str, split['PG']))
is_c = np.array([f in cactus_set for f in founders])
is_p = np.array([f in pg_set for f in founders])

ac_c = np.asarray(cn[is_c, :].sum(axis=0)).flatten()
ac_p = np.asarray(cn[is_p, :].sum(axis=0)).flatten()

THC = ${THRESH_C}
THP = ${THRESH_P}
keep = (ac_c >= THC) & (ac_p >= THP)
print(f'  threshold: cactus>={THC}, PG>={THP}', flush=True)
print(f'  keep: {keep.sum():,}/{K:,} ({keep.mean()*100:.2f}%)', flush=True)

cn_filt = cn.tocsc()[:, keep].tocsr()
print(f'  filtered cn nnz: {cn_filt.nnz:,}', flush=True)

# K_f check
Kf = np.asarray(cn_filt.sum(axis=1)).flatten()
print(f'  K_f median cactus={int(np.median(Kf[is_c])):,}, PG={int(np.median(Kf[is_p])):,}, ratio={np.median(Kf[is_c])/max(np.median(Kf[is_p]),1):.3f}', flush=True)

save_npz(OUT / f'cn_{CHR}.cn.npz', cn_filt)

# Copy & subset meta arrays that are aligned with k-mer dimension
new_meta = {}
for k in meta.files:
    a = meta[k]
    if a.ndim == 1 and a.shape[0] == K:
        new_meta[k] = a[keep]
    else:
        new_meta[k] = a
np.savez(OUT / f'cn_{CHR}.meta.npz', **new_meta)
print(f'saved {OUT}/cn_{CHR}.{{cn,meta}}.npz', flush=True)
EOF
ls -lh $OUT_DIR/
