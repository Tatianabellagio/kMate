#!/usr/bin/env python3
"""
Per-founder discriminating-tag supply from a cn_full matrix.
For each founder: total carried k-mers, and counts of LOW-AC (discriminating)
tags by carrier-count stratum. Low AC = carried by few founders = high
identifying power for h.

Writes TSV: founder  n_total  n_ac2  n_ac_le4  n_ac_le10
Usage: extract_per_founder_tags.py <cn_prefix> <out.tsv>
"""
import sys
import numpy as np
from scipy.sparse import load_npz

cn_prefix, out = sys.argv[1], sys.argv[2]
cn = load_npz(cn_prefix + ".cn.npz").tocsr()
meta = np.load(cn_prefix + ".meta.npz", allow_pickle=True)
founders = np.asarray(meta["founders"]).astype(str)
F, K = cn.shape
ac_k = np.asarray(cn.sum(axis=0)).flatten().astype(np.int32)
print(f"F={F} K={K:,}  ac min={ac_k.min()} max={ac_k.max()}", flush=True)

with open(out, "w") as fh:
    fh.write("founder\tn_total\tn_ac2\tn_ac_le4\tn_ac_le10\n")
    for f in range(F):
        s, e = cn.indptr[f], cn.indptr[f + 1]
        cols = cn.indices[s:e]
        a = ac_k[cols]
        fh.write(f"{founders[f]}\t{len(cols)}\t{int((a==2).sum())}\t"
                 f"{int((a<=4).sum())}\t{int((a<=10).sum())}\n")
print(f"wrote {out}", flush=True)
