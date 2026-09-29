#!/usr/bin/env python
"""Does a candidate's climate gradient survive removing the GENOME-WIDE garden profile?

38 bio15 candidates on different chromosomes share one garden profile (median |r| 0.69
between unlinked pairs vs 0.19 for random variants; one component explains 71%). That is
a genome-wide effect -- founders of some ancestry doing well in some gardens and carrying
their whole genome -- and it makes every variant they carry look climate-graded. A
per-locus gradient says nothing about THAT locus unless it survives removing it.

So: take the top K principal components of 3,000 random testable variants' garden
profiles (the genome-wide axes, the garden-level analogue of LFMM's latent factors),
regress each candidate's profile on them, and correlate the RESIDUAL with its own axis.
Background residuals are scored the same way, each on its best axis (selection-aware).

  r_resid, pct_resid   locus-specific gradient and its percentile
  gw_share             fraction of the candidate's garden-profile variance explained by
                       the genome-wide components (high = it mostly IS the genome-wide
                       pattern)
Writes the columns into results/candidate_evidence.csv.
"""
import os, sys
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(HERE, "..", ".."))); sys.path.insert(0, HERE)
import lib, axis_clusters as ac                                  # noqa: E402
import candidate_evidence as CE                                  # noqa: E402

OUT = f"{HERE}/results"


def main(K_PC=3):
    C = pd.read_csv(f"{OUT}/candidate_evidence.csv")
    PM = pd.read_csv(ac.POOLMETA); site = PM.site.values
    A = ac.site_climate(); A.index = sorted(PM.site.astype(int).unique())
    p0 = np.load(f"{lib.AF_STORE}/p0_nonsnp.npy")
    M = np.load(f"{lib.GEA}/common/results/pool_matrices/pool_gen9_nonsnp_af.npy", mmap_mode="r")
    rng = np.random.default_rng(6)
    B = pd.read_csv(f"{OUT}/repeat_variants_background_v2.csv.gz",
                    usecols=["chrom", "pos", "ref_len", "alt_len"])
    idx = np.load(f"{lib.AF_STORE}/index_nonsnp.npz", allow_pickle=True)
    K = pd.DataFrame({"chrom": idx["chrom"].astype(str), "pos": idx["pos"],
                      "ref_len": idx["ref_len"], "alt_len": idx["alt_len"]})
    K["row"] = np.arange(len(K))
    brow = np.sort(rng.choice(K.merge(B, on=["chrom", "pos", "ref_len", "alt_len"])["row"].values,
                              3000, replace=False))
    z = lambda G: (G - G.mean()) / G.std()
    Gb = z(CE.site_means(M, p0, brow, site)).fillna(0)
    U, S, _ = np.linalg.svd(Gb.values, full_matrices=False)       # sites x sites
    ev = S ** 2 / np.sum(S ** 2)
    P = U[:, :K_PC]                                               # genome-wide axes
    print("genome-wide garden components, variance explained:",
          ", ".join(f"{100*e:.0f}%" for e in ev[:5]))
    H = P @ np.linalg.pinv(P)                                     # projection
    def resid(G):
        X = z(G).fillna(0).values
        fit = H @ X
        return pd.DataFrame(X - fit, index=G.index, columns=G.columns), \
               (fit ** 2).sum(0) / np.maximum((X ** 2).sum(0), 1e-12)
    Rb, _ = resid(CE.site_means(M, p0, brow, site))
    bg = np.sort(CE.corr_axes(Rb, A).abs().max(axis=1).values)
    Gc = CE.site_means(M, p0, C.store_row.values, site); Gc.columns = C.index
    Rc, share = resid(Gc)
    C["r_resid"] = [abs(CE.corr_axes(Rc[[i]], A[[c.axis if c.axis in A.columns
                                                 else c.best_axis_all]]).iloc[0, 0])
                    for i, c in C.iterrows()]
    C["gw_share"] = share
    C["pct_resid"] = 100 * np.searchsorted(bg, C.r_resid.values) / len(bg)
    C.to_csv(f"{OUT}/candidate_evidence.csv", index=False)
    print("GPX6:", C.loc[C.sym == "GPX6", ["r_own", "gw_share", "r_resid", "pct_sel",
                                          "pct_resid"]].round(3).to_string(index=False))


if __name__ == "__main__":
    main()
