#!/usr/bin/env python
"""Collapse pool-level Δp to SITE level (mean Δp across pools within a site) and write an
LFMM stem, to test whether the K-invariant inflation floor is pseudoreplication.

If λ(pool, n=352) floors ~2 but λ(site, n=31) ~ 1, the floor is the design effect from
treating ~11 pools/site as independent -> fix is site-level / site-random-effect, not more K.
Usage: build_site_collapsed.py <class> <climate>
"""
import sys, os
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))))
import lib

CLS, CLIM = sys.argv[1], sys.argv[2]
CM = f"{lib.GEA}/phase1_replication/results/class_matrices"
OUT = f"{lib.GEA}/wza_investigation/deinflate/inputs"
_P0 = {"snp": "p0_snp.npy", "sv": "p0_nonsnp.npy", "smallindel": "p0_nonsnp.npy", "nonsnp": "p0_nonsnp.npy"}

recs = pd.read_csv(f"{CM}/{CLS}_gen9.records.csv")
pools = pd.read_csv(f"{CM}/gen9.pools.csv")
af = np.load(f"{CM}/{CLS}_gen9_af.npy")                       # [pools x rec]
p0 = np.load(f"{lib.AF_STORE}/{_P0[CLS]}")[recs["col"].to_numpy()]
dp = af - p0[None, :]                                         # Δp [pools x rec]

# impute NaN per record (column mean), same as build_lfmm_input
with np.errstate(invalid="ignore"):
    colmean = np.nanmean(dp, axis=0)
colmean = np.where(np.isfinite(colmean), colmean, 0.0)
miss = ~np.isfinite(dp); dp[miss] = np.take(colmean, np.where(miss)[1])

site = pools["site"].to_numpy()
usites = pd.unique(site)
# site-mean Δp
S = np.vstack([dp[site == s].mean(0) for s in usites])       # [n_site x rec]
xclim = pools.groupby("site")[CLIM].first().reindex(usites).to_numpy(float)
z = (xclim - xclim.mean()) / xclim.std(ddof=0)

n_site, n_rec = S.shape
print(f"{CLS} {CLIM}: collapsed {len(pools)} pools -> {n_site} sites x {n_rec:,} records "
      f"(mean {len(pools)/n_site:.1f} pools/site)")
stem = f"{OUT}/lfmm_{CLS}_site_gen9"
S.astype(np.float64).tofile(f"{stem}_Y.f64")
open(f"{stem}_dims.txt", "w").write(f"{n_site} {n_rec}\n")
pd.DataFrame({CLIM: z}).to_csv(f"{stem}_env.csv", index=False)
print(f"  -> {stem}_Y.f64 / _dims.txt / _env.csv")
