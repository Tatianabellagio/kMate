#!/usr/bin/env python
"""PC1 block test: collapse a block's [pools x SNPs] AF matrix to its dominant
haplotype axis (PC1), correlate PC1 with climate, site-permutation p. One test per
block regardless of m -> structurally removes the SNP-count pseudo-replication.

Choices (per user spec):
  - NaN imputed with the SNP's across-pool mean
  - correlation-PCA (each SNP scaled to unit variance) so rare/common contribute alike
  - report PC1 variance-explained (diagnostic: is the block one haplotype?)
Validated here on the two candidate blocks (gen1, bio1); compares PC1 vs block-mean.
"""
from __future__ import annotations
import numpy as np, pandas as pd
from scipy.stats import kendalltau

CM = ("/global/scratch/users/tbellg/kmate/analysis/grenenet_gea/"
      "phase1_replication/class_matrices")
OUT = "/global/scratch/users/tbellg/kmate/analysis/grenenet_gea/r2_gea_nonsnp/wza_investigation/results"
RNG = np.random.default_rng(0)
BLOCKS = {"4_2519": "Chr4 CRK", "2_1265": "CAM5"}


def pc1_scores(M):
    """M = pools x SNPs. NaN->col mean, scale to unit var, return PC1 scores + VE."""
    M = M.copy()
    col_mean = np.nanmean(M, axis=0)
    inds = np.where(np.isnan(M))
    M[inds] = np.take(col_mean, inds[1])
    M = M - M.mean(0)
    sd = M.std(0); sd[sd == 0] = 1.0
    M = M / sd                                   # correlation-PCA
    U, S, Vt = np.linalg.svd(M, full_matrices=False)
    ve = S[0]**2 / (S**2).sum()
    return U[:, 0] * S[0], ve, Vt[0]


def site_perm_p(score, bio1, site, n=20000):
    ok = np.isfinite(score) & np.isfinite(bio1)
    a, b, s = score[ok], bio1[ok], site[ok]
    obs = kendalltau(b, a).statistic
    sites = np.unique(s)
    site_bio = {si: b[s == si][0] for si in sites}
    cnt = 0
    for _ in range(n):
        perm = RNG.permutation(sites)
        mp = {sites[i]: site_bio[perm[i]] for i in range(len(sites))}
        bp = np.array([mp[si] for si in s])
        if abs(kendalltau(bp, a).statistic) >= abs(obs):
            cnt += 1
    return obs, (cnt + 1) / (n + 1)


def main():
    pools = pd.read_csv(f"{CM}/gen1.pools.csv")
    bio1 = pools["bio1"].to_numpy(); site = pools["site"].to_numpy()
    rec = pd.read_csv(f"{CM}/snp_gen1.records.csv", usecols=["block"])
    af = np.load(f"{CM}/snp_gen1_af.npy", mmap_mode="r")

    rows = []
    for blk, lab in BLOCKS.items():
        idx = np.where(rec["block"].to_numpy() == blk)[0]
        M = af[:, idx].astype(np.float64)
        pc1, ve, load = pc1_scores(M)
        bmean = np.nanmean(M, axis=1)
        # align PC1 sign so its correlation with block-mean is positive (interpretability)
        if np.corrcoef(pc1, np.nan_to_num(bmean - np.nanmean(bmean)))[0, 1] < 0:
            pc1 = -pc1
        tau_pc1, p_pc1 = site_perm_p(pc1, bio1, site)
        tau_bm, p_bm = site_perm_p(bmean, bio1, site)
        r_pc1_bm = np.corrcoef(pc1, bmean)[0, 1]
        rows.append(dict(block=blk, label=lab, m_SNPs=len(idx),
                         PC1_var_explained=round(ve, 3),
                         tau_PC1=round(tau_pc1, 3), site_perm_p_PC1=f"{p_pc1:.2e}",
                         tau_blockmean=round(tau_bm, 3), site_perm_p_blockmean=f"{p_bm:.2e}",
                         corr_PC1_blockmean=round(r_pc1_bm, 3)))
    res = pd.DataFrame(rows)
    res.to_csv(f"{OUT}/pc1_block_summary.csv", index=False)
    print(res.to_string(index=False))


if __name__ == "__main__":
    main()
