#!/usr/bin/env python
"""Are within-site plots independent for a CLIMATE association? Measure it.

Two honest numbers for the candidate blocks (block-mean AF, gen1):
  1. ICC (intraclass correlation) of block-mean AF across plots within site, and the
     implied design effect / effective N. Low ICC -> plots carry independent info
     (user's view); high ICC -> plots are near-replicates.
  2. SITE-permutation null: permute bio1 across the 31 sites (the level climate
     actually varies), keep all 326 pools, recompute Kendall tau. This uses ALL
     plots but counts only the independent climate contrasts -> the honest p.
Compares the honest p to the naive pool-level Kendall p.
"""
from __future__ import annotations
import os
import numpy as np, pandas as pd
from scipy.stats import kendalltau

CM = ("/global/scratch/users/tbellg/kmate/analysis/grenenet_gea/"
      "phase1_replication/class_matrices")
KEN = ("/global/scratch/users/tbellg/kmate/analysis/grenenet_gea/"
       "phase1_replication/kendall/kendall_snp_gen1_bio1.csv")
OUT = "/global/scratch/users/tbellg/kmate/analysis/grenenet_gea/r2_gea_nonsnp/wza_investigation/results"
RNG = np.random.default_rng(0)
BLOCKS = {"4_2519": "Chr4 CRK", "2_1265": "CAM5"}


def icc_oneway(values, groups):
    """One-way random-effects ICC(1) + effective N (design effect)."""
    df = pd.DataFrame({"y": values, "g": groups}).dropna()
    grand = df["y"].mean()
    gb = df.groupby("g")["y"]
    ni = gb.size().to_numpy()
    k = len(ni); N = ni.sum()
    msb = (gb.mean().sub(grand).pow(2) * ni).sum() / (k - 1)
    msw = (gb.apply(lambda x: ((x - x.mean())**2).sum())).sum() / (N - k)
    m0 = (N - (ni**2).sum() / N) / (k - 1)
    icc = (msb - msw) / (msb + (m0 - 1) * msw)
    icc = max(icc, 0.0)
    deff = 1 + (N / k - 1) * icc          # design effect at avg cluster size
    return dict(icc=round(icc, 3), n_pools=int(N), n_sites=k,
                n_eff=round(N / deff, 1), deff=round(deff, 2))


def site_perm_p(af_pool, bio1, site, n=20000):
    ok = np.isfinite(af_pool) & np.isfinite(bio1)
    a, b, s = af_pool[ok], bio1[ok], site[ok]
    obs = kendalltau(b, a).statistic
    # one climate value per site
    sites = np.unique(s)
    site_bio = {si: b[s == si][0] for si in sites}
    assert all((b[s == si] == site_bio[si]).all() for si in sites), "climate not site-constant"
    cnt = 0
    for _ in range(n):
        perm = RNG.permutation(sites)
        mapping = {sites[i]: site_bio[perm[i]] for i in range(len(sites))}
        bp = np.array([mapping[si] for si in s])
        if abs(kendalltau(bp, a).statistic) >= abs(obs):
            cnt += 1
    return obs, (cnt + 1) / (n + 1)


def main():
    pools = pd.read_csv(f"{CM}/gen1.pools.csv")
    bio1 = pools["bio1"].to_numpy(); site = pools["site"].to_numpy()
    rec = pd.read_csv(f"{CM}/snp_gen1.records.csv", usecols=["block"])
    ken = pd.read_csv(KEN, usecols=["block"])
    af = np.load(f"{CM}/snp_gen1_af.npy", mmap_mode="r")

    rows = []
    for blk, lab in BLOCKS.items():
        idx = np.where(rec["block"].to_numpy() == blk)[0]
        pool_af = np.nanmean(af[:, idx].astype(np.float64), axis=1)
        ic = icc_oneway(pool_af, site)
        # naive pool-level p
        ok = np.isfinite(pool_af)
        _, pool_p = kendalltau(bio1[ok], pool_af[ok])
        obs, perm_p = site_perm_p(pool_af, bio1, site)
        rows.append(dict(block=blk, label=lab, **ic, tau=round(obs, 3),
                         pool_kendall_p=f"{pool_p:.2e}", site_perm_p=f"{perm_p:.2e}"))
    res = pd.DataFrame(rows)
    res.to_csv(f"{OUT}/site_structure_summary.csv", index=False)
    print(res.to_string(index=False))


if __name__ == "__main__":
    main()
