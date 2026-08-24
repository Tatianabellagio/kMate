#!/usr/bin/env python
"""Test the WZA pseudo-replication / panel-density hypothesis on block 4_2519.

WZA's denominator is sqrt(Σ pq²) — the INDEPENDENCE assumption (drops cross-SNP
correlation). For N SNPs in tight LD the correct Stouffer-with-correlation
denominator is sqrt(pq' R pq); dropping R inflates Z by ~sqrt(N/M_eff). So a denser
panel that packs more *correlated* SNPs into a high-LD block inflates the block-Z.

We measure, for the Chr4 peak 4_2519 (1169 SNPs) and CAM5 2_1265 (13 SNPs):
  - within-block LD (mean |r| from the 326-pool AF matrix)
  - effective # independent SNPs M_eff (Li & Ji 2005, eigenvalues of R)
  - WZA-Z (independence) vs correlated-Z (uses R) -> the inflation factor
  - how raw WZA-Z scales as we subsample N SNPs (pseudo-replication curve)
"""
from __future__ import annotations
import os, sys
import numpy as np, pandas as pd
from scipy.stats import norm
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

CM = ("/global/scratch/users/tbellg/kmate/analysis/grenenet_selection/"
      "phase1_replication/class_matrices")
KEN = ("/global/scratch/users/tbellg/kmate/analysis/grenenet_selection/"
       "phase1_replication/kendall/kendall_snp_gen1_bio1.csv")
OUT = "/global/scratch/users/tbellg/kmate/analysis/grenenet_selection/wza/investigation/results"
RNG = np.random.default_rng(0)


def corr_matrix(A):
    """Pairwise Pearson corr across rows (pools) for columns (SNPs), NaN-aware via
    pairwise-complete means. A is pools x SNPs. Returns SNP x SNP R."""
    M = np.ma.masked_invalid(A)
    R = np.ma.corrcoef(M, rowvar=False)
    return np.asarray(R.filled(0.0))


def m_eff(R):
    """Li & Ji (2005) effective number of independent tests from eigenvalues."""
    ev = np.linalg.eigvalsh(R)
    ev = np.clip(ev, 0, None)
    return float(np.sum((ev >= 1) + (ev - np.floor(ev)) * (ev > 0)))


def analyze(block, rec, ken, af):
    # af columns align POSITIONALLY with records rows (both N=1,989,384 in order);
    # the `col` field is a stale pre-filter global index, not the af column.
    idx = np.where(rec["block"].to_numpy() == block)[0]
    sub = af[:, idx].astype(np.float64)                   # pools x SNPs
    chk = np.nanmean(sub, axis=0)
    assert np.allclose(chk, rec["p_bar"].to_numpy()[idx], atol=1e-6, equal_nan=True), "row misalign"
    k = ken.iloc[idx]                                     # same row order
    assert (k["block"] == block).all(), "kendall/records order mismatch"
    p = k["pval"].clip(1e-15, 1 - 1e-3).to_numpy()
    z = norm.ppf(1 - p)
    maf = k["MAF"].to_numpy(); pq = maf * (1 - maf)
    N = len(z)
    wza_num = (pq * z).sum()
    wza_Z = wza_num / np.sqrt((pq**2).sum())              # independence denom
    R = corr_matrix(sub)
    meanabs_r = (np.abs(R)[np.triu_indices(N, 1)]).mean() if N > 1 else np.nan
    corr_den = np.sqrt(pq @ R @ pq)                       # Stouffer-with-correlation
    corr_Z = wza_num / corr_den
    Meff = m_eff(R)
    return dict(block=block, N=N, mean_abs_r=round(meanabs_r, 3), M_eff=round(Meff, 1),
                WZA_Z_indep=round(wza_Z, 2), corr_Z=round(corr_Z, 2),
                inflation=round(wza_Z / corr_Z, 2),
                sqrt_N_over_Meff=round(np.sqrt(N / max(Meff, 1)), 2)), (z, pq, sub)


def subsample_curve(z, pq, label, grid):
    rows = []
    N = len(z)
    for n in grid:
        if n > N: continue
        zz = []
        for _ in range(60):
            idx = RNG.choice(N, n, replace=False)
            zz.append((pq[idx] * z[idx]).sum() / np.sqrt((pq[idx]**2).sum()))
        rows.append((label, n, np.mean(zz)))
    return rows


def main():
    rec = pd.read_csv(f"{CM}/snp_gen1.records.csv", usecols=["p_bar", "block"])
    ken = pd.read_csv(KEN, usecols=["pval", "MAF", "block"])
    assert len(rec) == len(ken), "records/kendall length mismatch"
    af = np.load(f"{CM}/snp_gen1_af.npy", mmap_mode="r")
    print("AF matrix:", af.shape)

    res = []
    curves = []
    for blk, grid in [("4_2519", [10, 20, 50, 100, 200, 400, 700, 1169]),
                      ("2_1265", [2, 4, 8, 13])]:
        info, (z, pq, sub) = analyze(blk, rec, ken, af)
        res.append(info)
        curves += subsample_curve(z, pq, blk, grid)
    summ = pd.DataFrame(res)
    summ.to_csv(f"{OUT}/pseudoreplication_summary.csv", index=False)
    print("\n=== WITHIN-BLOCK LD & PSEUDO-REPLICATION ===")
    print(summ.to_string(index=False))
    cur = pd.DataFrame(curves, columns=["block", "n", "raw_WZA_Z"])
    cur.to_csv(f"{OUT}/pseudoreplication_subsample.csv", index=False)

    # figure
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(13, 5))
    c = cur[cur.block == "4_2519"]
    ax[0].plot(c.n, c.raw_WZA_Z, "o-", label="observed raw WZA-Z (subsampled)")
    n0, z0 = c.n.iloc[0], c.raw_WZA_Z.iloc[0]
    ax[0].plot(c.n, z0*np.sqrt(c.n/n0), "k--", label="∝ √N (perfect-LD expectation)")
    ax[0].plot(c.n, z0*(c.n/n0), "r:", label="∝ N (full pseudo-replication)")
    ax[0].set(xlabel="# SNPs in block (subsampled)", ylabel="raw WZA-Z",
              title="Block 4_2519 (Chr4 CRK): Z grows with SNP count")
    ax[0].legend(fontsize=8)
    # inflation bar
    s = summ.set_index("block")
    ax[1].bar(["4_2519\n(1169 SNPs)", "2_1265\nCAM5 (13)"],
              [s.loc["4_2519", "inflation"], s.loc["2_1265", "inflation"]],
              color=["firebrick", "steelblue"])
    ax[1].axhline(1, color="k", lw=1)
    ax[1].set(ylabel="WZA-Z / correlated-Z", title="LD inflation of WZA-Z (independence vs R-aware)")
    for i, b in enumerate(["4_2519", "2_1265"]):
        ax[1].text(i, s.loc[b, "inflation"]+.05, f"M_eff={s.loc[b,'M_eff']}\nof N={s.loc[b,'N']}",
                   ha="center", fontsize=9)
    fig.tight_layout(); fig.savefig(f"{OUT}/plots/fig10_pseudoreplication.png", dpi=130)
    print(f"\nwrote fig10_pseudoreplication.png")


if __name__ == "__main__":
    main()
