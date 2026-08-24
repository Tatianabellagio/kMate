#!/usr/bin/env python
"""Step back from WZA: look at the RAW allele-frequency vs climate signal for the
two candidate blocks, before any p-value pooling.

Two questions:
  1. Do the block's allele frequencies actually track climate (bio1)? -> block-mean
     AF vs bio1 across pools/sites, with the sign/strength.
  2. Is the block ONE coherent signal (so the SNPs are NOT independent)? -> the
     distribution of per-SNP Kendall-tau within the block. A tight one-sided tau
     distribution = one haplotype's signal replicated across correlated SNPs
     (exactly why WZA's independence denominator over-counts).
"""
from __future__ import annotations
import os, sys
import numpy as np, pandas as pd
from scipy.stats import kendalltau
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

CM = ("/global/scratch/users/tbellg/kmate/analysis/grenenet_selection/"
      "phase1_replication/class_matrices")
KEN = ("/global/scratch/users/tbellg/kmate/analysis/grenenet_selection/"
       "phase1_replication/kendall/kendall_snp_gen1_bio1.csv")
OUT = "/global/scratch/users/tbellg/kmate/analysis/grenenet_selection/r2_gea_nonsnp/wza_investigation/results"
BLOCKS = {"4_2519": "Chr4 CRK cluster (top block)", "2_1265": "CAM5"}


def main():
    pools = pd.read_csv(f"{CM}/gen1.pools.csv")
    bio1 = pools["bio1"].to_numpy()
    site = pools["site"].to_numpy()
    rec = pd.read_csv(f"{CM}/snp_gen1.records.csv", usecols=["block"])
    ken = pd.read_csv(KEN, usecols=["pos", "tau", "pval", "block"])
    af = np.load(f"{CM}/snp_gen1_af.npy", mmap_mode="r")
    assert af.shape[0] == len(pools), "af rows must align with pools.csv"

    fig, axes = plt.subplots(2, 2, figsize=(14, 9))
    summ = []
    for col, (blk, title) in enumerate(BLOCKS.items()):
        idx = np.where(rec["block"].to_numpy() == blk)[0]
        sub = af[:, idx].astype(np.float64)            # pools x SNPs
        ktau = ken["tau"].to_numpy()[idx]              # per-SNP tau vs bio1 (precomputed)

        # (1) block-mean AF per pool, vs bio1; aggregate to site means too
        mean_af = np.nanmean(sub, axis=0)              # per-SNP mean (unused)
        pool_af = np.nanmean(sub, axis=1)              # per-pool block-mean AF
        ok = np.isfinite(pool_af) & np.isfinite(bio1)
        tb, pb = kendalltau(bio1[ok], pool_af[ok])
        sdf = pd.DataFrame({"site": site[ok], "bio1": bio1[ok], "af": pool_af[ok]})
        site_means = sdf.groupby("site").mean()
        ts, ps = kendalltau(site_means["bio1"], site_means["af"])

        ax = axes[0, col]
        ax.scatter(bio1[ok], pool_af[ok], s=12, alpha=.3, color="grey", label="pools")
        ax.scatter(site_means["bio1"], site_means["af"], s=70, color="firebrick",
                   edgecolor="k", zorder=5, label="site means")
        ax.set(xlabel="bio1 (annual mean temp)", ylabel="block-mean ALT freq",
               title=f"{blk}  {title}\nblock-mean AF vs climate: pool τ={tb:.2f} (p={pb:.1e}); site τ={ts:.2f}")
        ax.legend(fontsize=8)

        # (2) per-SNP tau distribution -> coherence
        ax2 = axes[1, col]
        ax2.hist(ktau, bins=40, color="steelblue")
        ax2.axvline(0, color="red", ls="--")
        ax2.axvline(np.nanmedian(ktau), color="k", lw=2, label=f"median τ={np.nanmedian(ktau):.2f}")
        frac_neg = np.mean(ktau < 0)
        ax2.set(xlabel="per-SNP Kendall τ (AF vs bio1)", ylabel="# SNPs",
                title=f"within-block per-SNP τ — coherence\n{frac_neg*100:.0f}% negative, "
                      f"|τ| median={np.nanmedian(np.abs(ktau)):.2f}")
        ax2.legend(fontsize=8)
        summ.append(dict(block=blk, N=len(idx), pool_tau=round(tb, 3), pool_p=pb,
                         site_tau=round(ts, 3), site_p=ps,
                         frac_snps_same_sign=round(max(frac_neg, 1-frac_neg), 3),
                         median_abs_snp_tau=round(np.nanmedian(np.abs(ktau)), 3),
                         n_sites=site_means.shape[0]))
    fig.tight_layout(); fig.savefig(f"{OUT}/fig11_raw_signal.png", dpi=130)
    s = pd.DataFrame(summ); s.to_csv(f"{OUT}/raw_signal_summary.csv", index=False)
    print(s.to_string(index=False))
    print(f"\nwrote fig11_raw_signal.png")


if __name__ == "__main__":
    main()
