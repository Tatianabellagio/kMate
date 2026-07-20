#!/usr/bin/env python
"""Genome-wide PC1 block test with a SITE-level permutation null (gen1, bio1).

Per block: collapse [pools x SNPs] AF to PC1 (dominant haplotype axis, correlation-
PCA, NaN->col mean), average PC1 to the 31 sites, and test the 31 site-mean-PC1 vs
31 site-climate with a Spearman statistic. Null = permute climate across the 31
sites (the level it varies). One test per block -> no SNP-count pseudo-replication;
site permutation -> no pool-within-site pseudo-replication. BH across blocks.

Writes pc1_genomewide_gen1_bio1.csv: block, chrom, pos, m, pc1_ve, r_obs, perm_p, bh_q.
"""
from __future__ import annotations
import numpy as np, pandas as pd
from scipy.stats import rankdata

CM = ("/global/scratch/users/tbellg/kmate/analysis/grenenet_gea/"
      "phase1_replication/class_matrices")
OUT = "/global/scratch/users/tbellg/kmate/analysis/grenenet_gea/wza_investigation/results"
NPERM = 20000
RNG = np.random.default_rng(0)


def pc1_site_means(M, site_idx, n_sites):
    """M = pools x SNPs -> (PC1 pool scores aggregated to site means [n_sites], VE)."""
    col_mean = np.nanmean(M, axis=0)
    nan = np.isnan(M)
    if nan.any():
        M = np.where(nan, col_mean, M)
    M = M - M.mean(0)
    sd = M.std(0); sd[sd == 0] = 1.0
    M = M / sd
    G = M @ M.T                                   # 326x326 Gram
    w, V = np.linalg.eigh(G)
    pc1 = V[:, -1] * np.sqrt(max(w[-1], 0))       # top eigenvector -> pool scores
    ve = w[-1] / w.sum() if w.sum() > 0 else np.nan
    sm = np.zeros(n_sites); cnt = np.zeros(n_sites)
    np.add.at(sm, site_idx, pc1); np.add.at(cnt, site_idx, 1)
    return sm / np.maximum(cnt, 1), float(ve)


def main():
    pools = pd.read_csv(f"{CM}/gen1.pools.csv")
    bio1 = pools["bio1"].to_numpy()
    sites, site_idx = np.unique(pools["site"].to_numpy(), return_inverse=True)
    n_sites = len(sites)
    # one climate value per site
    site_clim = np.array([bio1[site_idx == i][0] for i in range(n_sites)])
    print(f"{len(pools)} pools, {n_sites} sites")

    rec = pd.read_csv(f"{CM}/snp_gen1.records.csv", usecols=["chrom", "pos", "block"])
    af = np.load(f"{CM}/snp_gen1_af.npy", mmap_mode="r")
    blocks = rec["block"].to_numpy()
    # contiguous block ranges
    order = pd.RangeIndex(len(blocks))
    grp = rec.groupby("block", sort=False)

    rows = []
    X = []           # ranked site-mean PC1 per block (for vectorized null)
    keep = []
    for blk, g in grp:
        idx = g.index.to_numpy()
        M = af[:, idx].astype(np.float64)
        sm, ve = pc1_site_means(M, site_idx, n_sites)
        # sign-align PC1 so it positively tracks the block-mean AF (interpretability)
        bm = np.nanmean(M, axis=1)
        bms = np.zeros(n_sites); c = np.zeros(n_sites)
        np.add.at(bms, site_idx, np.nan_to_num(bm)); np.add.at(c, site_idx, 1)
        bms /= np.maximum(c, 1)
        if np.corrcoef(sm, bms)[0, 1] < 0:
            sm = -sm
        rows.append((blk, g["chrom"].iloc[0], int(np.median(g["pos"])), len(idx), ve))
        X.append(rankdata(sm))
        keep.append(blk)
    X = np.asarray(X)                                       # n_blocks x n_sites (ranks)
    meta = pd.DataFrame(rows, columns=["block", "chrom", "pos", "m", "pc1_ve"])
    print(f"{len(meta)} blocks; PC1 computed")

    # vectorized site-permutation Spearman null
    cr = rankdata(site_clim)
    Xs = (X - X.mean(1, keepdims=True)); Xs /= np.linalg.norm(Xs, axis=1, keepdims=True)
    cs = (cr - cr.mean()); cs /= np.linalg.norm(cs)
    r_obs = Xs @ cs                                         # n_blocks
    perms = np.array([RNG.permutation(cr) for _ in range(NPERM)])   # NPERM x n_sites
    Ps = perms - perms.mean(1, keepdims=True); Ps /= np.linalg.norm(Ps, axis=1, keepdims=True)
    R_null = Xs @ Ps.T                                      # n_blocks x NPERM
    perm_p = (1 + (np.abs(R_null) >= np.abs(r_obs)[:, None]).sum(1)) / (NPERM + 1)
    meta["r_obs"] = r_obs
    meta["perm_p"] = perm_p
    # BH
    o = np.argsort(perm_p); n = len(perm_p)
    q = np.empty(n); q[o] = np.minimum.accumulate((perm_p[o] * n / np.arange(1, n + 1))[::-1])[::-1]
    meta["bh_q"] = q
    meta.to_csv(f"{OUT}/pc1_genomewide_gen1_bio1.csv", index=False)

    print(f"\nBH q<0.05: {(q < 0.05).sum()} blocks | min perm_p {perm_p.min():.2e}")
    print("\ntop 12 blocks:")
    print(meta.nsmallest(12, "perm_p")[["block", "chrom", "pos", "m", "pc1_ve", "r_obs", "perm_p", "bh_q"]].to_string(index=False))
    print("\ncandidate blocks:")
    print(meta[meta.block.isin(["4_2519", "2_1265"])][["block", "chrom", "pos", "m", "pc1_ve", "r_obs", "perm_p", "bh_q"]].to_string(index=False))


if __name__ == "__main__":
    main()
