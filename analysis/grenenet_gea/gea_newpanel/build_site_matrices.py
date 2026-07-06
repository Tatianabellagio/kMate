#!/usr/bin/env python
"""Site-level LFMM inputs: aggregate the 355 plot pools to 31 SITE means
(flower-weighted), to test whether the honest unit (one climate value = one
independent observation) changes power vs the plot-level (pseudoreplicated) run.

For cls in {snp, nonsnp}: Y = Δp = (flower-weighted site-mean AF) - p0, [31 x rec].
Env = per-site bio (all bio1..bio19), standardized. Writes the same _Y.f64 /
_dims.txt / _env_{bio}.csv layout the LFMM R runner consumes.
"""
from __future__ import annotations
import os, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import lib

CMDIR = f"{lib.GEA}/phase1_replication/class_matrices"
OUT = f"{lib.GEA}/gea_newpanel/lfmm_site"
_P0 = {"snp": "p0_snp.npy", "nonsnp": "p0_nonsnp.npy"}
BIOS = [f"bio{i}" for i in range(1, 20)]


def main():
    os.makedirs(OUT, exist_ok=True)
    pools = pd.read_csv(f"{CMDIR}/gen9.pools.csv")
    sites = np.sort(pools.site.unique())
    w = pools.total_flowers.to_numpy(float)
    # site membership + per-site flower weights
    site_of = pools.site.to_numpy()
    for cls in ["snp", "nonsnp"]:
        recs = pd.read_csv(f"{CMDIR}/{cls}_gen9.records.csv")
        af = np.load(f"{CMDIR}/{cls}_gen9_af.npy")                 # [355 x rec]
        p0 = np.load(f"{lib.AF_STORE}/{_P0[cls]}")[recs['col'].to_numpy()]
        # flower-weighted site-mean AF
        site_af = np.empty((len(sites), af.shape[1]), dtype=np.float64)
        for i, s in enumerate(sites):
            m = site_of == s
            ws = w[m]; ws = ws / ws.sum()
            site_af[i] = ws @ af[m].astype(np.float64)
        dp = site_af - p0[None, :]                                 # [31 x rec]
        stem = f"{OUT}/lfmm_{cls}_site"
        dp.tofile(f"{stem}_Y.f64")
        open(f"{stem}_dims.txt", "w").write(f"{len(sites)} {dp.shape[1]}\n")
        print(f"{cls}: site Δp [{dp.shape[0]} x {dp.shape[1]:,}] -> {stem}_Y.f64", flush=True)

    # per-site env (one row per site), standardized, all bios
    sp = pools.drop_duplicates("site").set_index("site").loc[sites]
    for bio in BIOS:
        if bio in sp:
            x = sp[bio].to_numpy(float); z = (x - x.mean()) / x.std(ddof=0)
            pd.DataFrame({bio: z}).to_csv(f"{OUT}/env_site_{bio}.csv", index=False)
    print(f"env: {sum(b in sp for b in BIOS)} bios, {len(sites)} sites -> {OUT}/env_site_*.csv")


if __name__ == "__main__":
    main()
