#!/usr/bin/env python
"""Shared inputs for the LFMM power experiments.

Builds, for cls in {snp, nonsnp}:
  * SITE-level Δp matrices (flower-weighted per-site mean, 31 sites) -> lfmm_site/lfmm_{cls}_site_Y.f64
    (the honest unit: one climate value = one independent observation)
and the ENV files (standardized) at BOTH resolutions:
  * site  (31 rows):  env_site/env_site_{bio}.csv   for bio1..19 and pc1
  * plot  (355 rows): env_plot/env_plot_{bio}.csv    for bio1..19 and pc1
PC1 = first PC of the 19 standardized bioclim vars (sign-oriented to +bio1).
Plot-level Y already exists (phase1_replication/lfmm/lfmm_{cls}_gen9_Y.f64).
"""
from __future__ import annotations
import os, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import lib

CMDIR = f"{lib.GEA}/phase1_replication/class_matrices"
BASE = f"{lib.GEA}/gea_newpanel"
SITE = f"{BASE}/lfmm_site"
ENV_SITE = f"{BASE}/env_site"
ENV_PLOT = f"{BASE}/env_plot"
_P0 = {"snp": "p0_snp.npy", "nonsnp": "p0_nonsnp.npy", "sv": "p0_nonsnp.npy"}
BIOS = [f"bio{i}" for i in range(1, 20)]


def pc1(mat):
    """First PC scores of columns=vars (rows=obs), standardized; +bio1-oriented."""
    Z = (mat - mat.mean(0)) / mat.std(0, ddof=0)
    U, S, Vt = np.linalg.svd(Z - Z.mean(0), full_matrices=False)
    sc = U[:, 0] * S[0]
    if np.corrcoef(sc, mat[:, 0])[0, 1] < 0:      # orient to +bio1
        sc = -sc
    return (sc - sc.mean()) / sc.std(ddof=0)


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--classes", nargs="+", default=["snp", "nonsnp", "sv"],
                    choices=["snp", "nonsnp", "sv"])
    ap.add_argument("--skip-env", action="store_true", help="don't rebuild env files")
    args = ap.parse_args()
    for d in (SITE, ENV_SITE, ENV_PLOT):
        os.makedirs(d, exist_ok=True)
    pools = pd.read_csv(f"{CMDIR}/gen9.pools.csv")
    sites = np.sort(pools.site.unique())
    w = pools.total_flowers.to_numpy(float)
    site_of = pools.site.to_numpy()
    biocols = [b for b in BIOS if b in pools.columns]

    # --- site-level Δp matrices ---
    for cls in args.classes:
        recs = pd.read_csv(f"{CMDIR}/{cls}_gen9.records.csv")
        af = np.load(f"{CMDIR}/{cls}_gen9_af.npy")
        p0 = np.load(f"{lib.AF_STORE}/{_P0[cls]}")[recs['col'].to_numpy()]
        site_af = np.empty((len(sites), af.shape[1]), dtype=np.float64)
        for i, s in enumerate(sites):
            m = site_of == s; ws = w[m] / w[m].sum()
            site_af[i] = ws @ af[m].astype(np.float64)
        dp = site_af - p0[None, :]
        stem = f"{SITE}/lfmm_{cls}_site"
        dp.tofile(f"{stem}_Y.f64")
        open(f"{stem}_dims.txt", "w").write(f"{len(sites)} {dp.shape[1]}\n")
        print(f"{cls}: site Δp [{dp.shape[0]} x {dp.shape[1]:,}]", flush=True)

    # --- envs, both resolutions, all bios + pc1 ---
    if args.skip_env:
        print("skip-env: envs unchanged"); return
    sp = pools.drop_duplicates("site").set_index("site").loc[sites]
    site_bio = sp[biocols].to_numpy(float)
    plot_bio = pools[biocols].to_numpy(float)
    for name, df_bio, outdir, tag in [(sites, site_bio, ENV_SITE, "site"),
                                      (pools, plot_bio, ENV_PLOT, "plot")]:
        for j, b in enumerate(biocols):
            x = df_bio[:, j]; z = (x - x.mean()) / x.std(ddof=0)
            pd.DataFrame({b: z}).to_csv(f"{outdir}/env_{tag}_{b}.csv", index=False)
        pd.DataFrame({"pc1": pc1(df_bio)}).to_csv(f"{outdir}/env_{tag}_pc1.csv", index=False)
    # report PC1 loadings (site level)
    Z = (site_bio - site_bio.mean(0)) / site_bio.std(0, ddof=0)
    _, _, Vt = np.linalg.svd(Z - Z.mean(0), full_matrices=False)
    load = pd.Series(Vt[0], index=biocols).sort_values(key=abs, ascending=False)
    print(f"envs: {len(biocols)} bios + pc1 at site(31) & plot(355)")
    print("PC1 top loadings:", ", ".join(f"{k}={v:+.2f}" for k, v in load.head(6).items()))


if __name__ == "__main__":
    main()
