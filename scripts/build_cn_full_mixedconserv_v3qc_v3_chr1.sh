#!/bin/bash
#SBATCH --job-name=mixcons_v3
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=2
#SBATCH --mem=48G
#SBATCH --time=1:00:00
#SBATCH --output=logs/mixcons_v3_%j.out
#SBATCH --error=logs/mixcons_v3_%j.err
# Subset cn_full_v3qc_v3 to k-mers carried by ≥10 cactus + ≥10 PG founder.
mkdir -p logs
set -euo pipefail
PY=/global/home/users/tbellg/miniforge3/envs/hapfm/bin/python
OUT_DIR=/global/scratch/users/tbellg/kmate/data/cn_full_231_v3qc_v3_mixedconserv
mkdir -p $OUT_DIR

$PY << EOF
import numpy as np, json
from scipy.sparse import load_npz, save_npz
from pathlib import Path
CHR = 'Chr1'
SRC = Path('/global/scratch/users/tbellg/kmate/data/cn_full_231_v3qc_v3')
OUT = Path('${OUT_DIR}')
cn = load_npz(SRC / f'cn_{CHR}.cn.npz').tocsr()
meta = np.load(SRC / f'cn_{CHR}.meta.npz', allow_pickle=True)
founders = np.asarray(meta['founders']).astype(str)
F, K = cn.shape
print(f'cn_full input: ({F}, {K:,}) nnz={cn.nnz:,}', flush=True)
with open('/global/scratch/users/tbellg/kmate/data/founder_split_cactus_pg.json') as fp:
    split = json.load(fp)
cactus = set(map(str, split['cactus'])); pg = set(map(str, split['PG']))
is_c = np.array([f in cactus for f in founders])
is_p = np.array([f in pg for f in founders])
ac_c = np.asarray(cn[is_c, :].sum(axis=0)).flatten()
ac_p = np.asarray(cn[is_p, :].sum(axis=0)).flatten()
keep = (ac_c >= 10) & (ac_p >= 10)
print(f'mixed-conserv keep: {keep.sum():,}/{K:,} ({keep.mean()*100:.2f}%)', flush=True)
cn_f = cn.tocsc()[:, keep].tocsr()
Kf = np.asarray(cn_f.sum(axis=1)).flatten()
print(f'K_f median cactus={int(np.median(Kf[is_c])):,}, PG={int(np.median(Kf[is_p])):,}, ratio={np.median(Kf[is_c])/max(np.median(Kf[is_p]),1):.3f}', flush=True)
save_npz(OUT / f'cn_{CHR}.cn.npz', cn_f)
new_meta = {}
for k in meta.files:
    a = meta[k]
    if a.ndim == 1 and a.shape[0] == K:
        new_meta[k] = a[keep]
    else:
        new_meta[k] = a
np.savez(OUT / f'cn_{CHR}.meta.npz', **new_meta)
print('saved', flush=True)
EOF
ls -lh $OUT_DIR/
