#!/usr/bin/env python
"""Flower-weighted merge of per-sample window AF into the 355 gen9 pools.

Mirrors build_pool_matrix's merge (p_pool = sum_t flowers_t p_t / sum_t flowers_t,
NaN-aware) but on the WINDOW-mode per-sample gen9-record vectors, in the EXACT pool
order of the global gen9.pools.csv so rows align to the global gen9 matrices.

Output: results/grenenet_gea/gen9_window/window_gen9_af.npy  [355 x 2,684,700] float32
(columns row-aligned to gen9_window/records.csv = same records as the global gen9 class
matrices, concatenated snp|sv|smallindel)."""
import os, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import lib

GW = f"{lib.GEA}/gen9_window"
CM = f"{lib.GEA}/phase1_replication/class_matrices"
pools = pd.read_csv(f"{CM}/gen9.pools.csv").pool.tolist()      # 355, fixed order
pt = lib.pool_table()
N = np.load(f"{GW}/gather_rows.npy").shape[0]
print(f"{len(pools)} pools x {N:,} records")

out = np.lib.format.open_memmap(f"{GW}/window_gen9_af.npy", mode="w+",
                                dtype=np.float32, shape=(len(pools), N))
cache = {}
def vec(s):
    if s not in cache:
        cache[s] = np.load(f"{GW}/persample/{s}.npy")            # float32, NaN
    return cache[s]

for i, pool in enumerate(pools):
    mem = pt[pt.pool == pool]
    acc = np.zeros(N); wsum = np.zeros(N)
    for s, w in zip(mem.sampleid, mem.flowerscollected.astype(float)):
        f = vec(s); ok = np.isfinite(f)
        acc[ok] += w * f[ok]; wsum[ok] += w
    out[i] = np.where(wsum > 0, acc / np.where(wsum > 0, wsum, 1), np.nan)
    if (i + 1) % 50 == 0:
        print(f"  {i+1}/{len(pools)} pools", flush=True)
out.flush()
print(f"wrote {GW}/window_gen9_af.npy  [{len(pools)} x {N:,}]")
