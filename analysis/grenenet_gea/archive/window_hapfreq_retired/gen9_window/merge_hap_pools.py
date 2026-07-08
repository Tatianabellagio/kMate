#!/usr/bin/env python
"""Flower-weighted site_gen_plot POOL matrices of haplotype-cluster frequencies, per
generation (the Pipeline-B analysis unit). Mirrors build_pool_matrix but on the per-sample
hap-freq vectors (hap/persample/{s}.npy, n_traj=38871). hap freqs have no missingness
(window h is always defined) so the merge is a plain flower-weighted mean.

Outputs (hap/pools/):
  pool_gen{g}_hap.npy   float32 [n_pools_g x 38871]
  pool_gen{g}.meta.csv  pool,site,plot,generation,n_timepoints,total_flowers,mean_coverage,bio1..19
"""
import os, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import lib

GW = "results/grenenet_gea/gen9_window"
PS = f"{GW}/hap/persample"
OUT = f"{GW}/hap/pools"
os.makedirs(OUT, exist_ok=True)
N = pd.read_csv(f"{GW}/hap/registry.csv").shape[0]
pt = lib.pool_table(); clim = lib.load_climate()
cache = {}
def vec(s):
    if s not in cache:
        cache[s] = np.load(f"{PS}/{s}.npy")
    return cache[s]

for g in sorted(pt.generation.unique()):
    sub = pt[pt.generation == g]
    pools = sorted(sub.pool.unique(), key=lambda p: tuple(int(x) for x in p.split("_")))
    mat = np.zeros((len(pools), N), np.float32); rows = []
    for i, pool in enumerate(pools):
        mem = sub[sub.pool == pool]
        acc = np.zeros(N); w = 0.0
        for s, fl in zip(mem.sampleid, mem.flowerscollected.astype(float)):
            acc += fl * vec(s); w += fl
        mat[i] = (acc / w) if w > 0 else np.nan
        rows.append(dict(pool=pool, site=int(mem.site.iloc[0]), plot=int(mem["plot"].iloc[0]),
                         generation=int(g), n_timepoints=len(mem),
                         total_flowers=float(mem.flowerscollected.sum()),
                         mean_coverage=float(mem.coverage.mean())))
    np.save(f"{OUT}/pool_gen{g}_hap.npy", mat)
    pd.DataFrame(rows).join(clim, on="site").to_csv(f"{OUT}/pool_gen{g}.meta.csv", index=False)
    print(f"gen{g}: {len(pools)} pools x {N:,} hap trajectories -> pool_gen{g}_hap.npy", flush=True)
