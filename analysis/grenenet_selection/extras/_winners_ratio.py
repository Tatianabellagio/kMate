#!/usr/bin/env python
"""Raw long-read (cactus) vs short-read (PanGenie) counts among the WINNING founders, per site, at a
few top-quantile cutoffs (robustness to the decile choice). Winners = founders with the highest
per-site selection coefficient s. Env: kmate."""
import os, sys, glob, json
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import lib
os.chdir("/global/scratch/users/tbellg/kmate")
HCACHE = f"{lib.GEA}/r3_persite_gwas/results/ecotype_fitness/sample_global_h.npz"
CHROMS = ["Chr1", "Chr2", "Chr3", "Chr4", "Chr5"]; EPS, P0_FLOOR = 1e-3, 1e-5
PURGE = [4, 32, 43, 60, 26]; TOPQ = [0.90, 0.95, 0.80]   # top 10% / 5% / 20%

def logit(p): p = np.clip(p, EPS, 1 - EPS); return np.log(p / (1 - p))
def seedmix_p0():
    ids = sorted({os.path.basename(p).split("_Chr")[0] for p in glob.glob(f"{lib.SEEDMIX}/*_Chr1.h_per_chrom.npz")})
    founders, reps = None, []
    for s in ids:
        chs = []
        for c in CHROMS:
            z = np.load(f"{lib.SEEDMIX}/{s}_{c}.h_per_chrom.npz", allow_pickle=True)
            if founders is None: founders = z["founders"].astype(str)
            chs.append(z[c].astype(float))
        reps.append(np.vstack(chs).mean(0))
    return founders, np.vstack(reps).mean(0)

founders, p0 = seedmix_p0(); p0f = np.maximum(p0, P0_FLOOR)
z = np.load(HCACHE, allow_pickle=True); Hs = z["H"].astype(float)
smap = {s: i for i, s in enumerate(z["samples"].astype(str))}
cactus = set(map(str, json.load(open("data/founder_split_cactus_pg.json"))["cactus"]))
is_cac = np.array([f in cactus for f in founders])
print(f"panel: {is_cac.sum()} cactus (long-read) : {(~is_cac).sum()} PG (short-read)  "
      f"= {is_cac.mean():.3f} cactus  (this is the 0.35 baseline)")

pt = lib.pool_table(); pt = pt[pt.sampleid.astype(str).isin(smap)].copy()
clim = lib.load_climate()
pool_h, pool_meta = {}, {}
for pool, g in pt.groupby("pool"):
    idx = [smap[str(s)] for s in g.sampleid]
    w = g.flowerscollected.to_numpy(float); w = np.where(np.isfinite(w) & (w > 0), w, 1.0); w /= w.sum()
    pool_h[pool] = (Hs[idx] * w[:, None]).sum(0)
    r = g.iloc[0]; pool_meta[pool] = (int(r["site"]), int(r["generation"]), int(r["plot"]))

def site_s(site):
    by_plot = {}
    for pool, (st, gen, plot) in pool_meta.items():
        if st == site and gen in (1, 2, 3): by_plot.setdefault(plot, {})[gen] = pool_h[pool]
    sl = []
    for plot, cells in by_plot.items():
        gens = sorted(cells)
        if 1 not in gens: continue
        t = np.array([0.0] + [float(g) for g in gens]); tc = t - t.mean()
        Y = np.vstack([p0f] + [cells[g] for g in gens])
        sl.append((tc[:, None] * logit(Y)).sum(0) / (tc @ tc))
    return np.mean(sl, 0) if sl else None

print("\ntop-decile (top 10%) WINNERS: cactus:PG counts + fraction cactus, per purge site")
print(f"{'site':>5} {'bio1':>5}  {'nWin':>4}  {'cactus':>6} {'PG':>4}  {'cac:PG':>7}  {'%cactus':>7}")
allc = allp = 0
for site in PURGE:
    s = site_s(site)
    if s is None: continue
    hi = s >= np.quantile(s, 0.90)
    nc = int(is_cac[hi].sum()); npg = int((~is_cac[hi]).sum()); allc += nc; allp += npg
    print(f"{site:>5} {clim.loc[site,'bio1']:>5.1f}  {hi.sum():>4}  {nc:>6} {npg:>4}  "
          f"{'1:%.1f'%(npg/max(nc,1)):>7}  {nc/hi.sum():>7.2f}")
print(f"{'POOL':>5} {'':>5}  {allc+allp:>4}  {allc:>6} {allp:>4}  {'1:%.1f'%(allp/max(allc,1)):>7}  {allc/(allc+allp):>7.2f}")

print("\nrobustness to the cutoff (pooled over the 5 purge sites):")
for q in TOPQ:
    c = p = 0
    for site in PURGE:
        s = site_s(site)
        if s is None: continue
        hi = s >= np.quantile(s, q); c += int(is_cac[hi].sum()); p += int((~is_cac[hi]).sum())
    lab = f"top {int(round((1-q)*100))}%"
    print(f"  {lab:>7}: {c} cactus : {p} PG  = 1:{p/max(c,1):.1f}  ({c/(c+p):.2f} cactus)  "
          f"vs baseline {is_cac.mean():.2f}")
