"""Build cn_full_231_v3_filt2/cn_Chr{N}.{cn,meta}.npz for one chrom by filtering
columns with ac_k < 2 from cn_full_231_v3.

Usage: python build_filt2_one_chrom.py <chrom_number>
"""
import sys, time
import numpy as np
from scipy.sparse import load_npz, save_npz
from pathlib import Path

chrom_num = int(sys.argv[1])
chrom = f'Chr{chrom_num}'
ROOT = Path('/carnegie/nobackup/scratch/tbellagio/hapfire_sv/poolfreq/data')
SRC_DIR = ROOT / 'cn_full_231_v3'
OUT_DIR = ROOT / 'cn_full_231_v3_filt2'
OUT_DIR.mkdir(exist_ok=True)

print(f'[{chrom}] building filt2...', flush=True)
t = time.time()
cn = load_npz(SRC_DIR / f'cn_{chrom}.cn.npz')
meta = np.load(SRC_DIR / f'cn_{chrom}.meta.npz', allow_pickle=True)
F, K = cn.shape
print(f'[{chrom}] loaded cn ({F}, {K:,}), nnz={cn.nnz:,}; meta keys={list(meta.keys())}', flush=True)

# Column carrier count ac_k
ac = np.asarray(cn.sum(axis=0)).flatten().astype(np.int32)
keep = ac >= 2
print(f'[{chrom}] keep ac>=2: {keep.sum():,} / {K:,} '
      f'(ac=0: {(ac==0).sum():,}; ac=1: {(ac==1).sum():,})', flush=True)

# Filter cn
cn_f = cn.tocsc()[:, keep].tocsr()
print(f'[{chrom}] filtered cn shape: {cn_f.shape}, nnz: {cn_f.nnz:,}', flush=True)

# Filter meta — kmer_index and bubble_id are per-kmer; bubble_chrom/start/end are per-bubble (unchanged)
new_meta = {}
for k in meta.keys():
    a = meta[k]
    if a.shape and a.shape[0] == K:
        new_meta[k] = a[keep]
    else:
        new_meta[k] = a

# Save
out_cn = OUT_DIR / f'cn_{chrom}.cn.npz'
out_meta = OUT_DIR / f'cn_{chrom}.meta.npz'
save_npz(out_cn, cn_f)
np.savez(out_meta, **new_meta)
print(f'[{chrom}] saved {out_cn} + {out_meta} in {time.time()-t:.0f}s', flush=True)
