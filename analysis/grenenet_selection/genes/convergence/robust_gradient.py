#!/usr/bin/env python
"""Climate gradients that survive losing the axis's high-leverage gardens.

bio15 is not spread evenly over the 31 gardens: 27 sit at 12-65 and four jump to 87-101
(gardens 4, 43, 26, 60). Two of those are the garden-4 founder-sweep set the review brief
grades WEAK regardless of axis, and two only reached generation 1. So a bio15 "gradient"
is largely a contrast between those four and the rest, and leave-ONE-out cannot see it --
the other three still carry it. Every axis gets the same treatment here:

  rho_all     Spearman rho over all 31 gardens (ranks, so a tail cannot dominate)
  r_trim      Pearson r after dropping the axis's high-leverage gardens (|z| > 2), jointly
  n_trim      how many gardens that dropped
  pct_rho / pct_trim   percentiles against 3,000 matched background records, each on its
              own best axis and computed the same way (selection-aware, like-for-like)

Writes the columns into results/candidate_evidence.csv.
"""
import os, sys
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(HERE, "..", ".."))); sys.path.insert(0, HERE)
import lib, axis_clusters as ac                                  # noqa: E402
import candidate_evidence as CE                                  # noqa: E402

OUT = f"{HERE}/results"


def lever(a):
    z = (a - a.mean()) / a.std()
    return list(a.index[np.abs(z) > 2])


def stats(G, a):
    """G sites x vars, a site vector -> (|rho|, |r_trim|) per var."""
    a = a.reindex(G.index)
    rk = G.rank(); rho = CE.corr_axes(rk, a.rank().to_frame("x"))["x"].abs().values
    keep = [s for s in G.index if s not in lever(a)]
    rt = CE.corr_axes(G.loc[keep], a.loc[keep].to_frame("x"))["x"].abs().values
    return rho, rt


def main():
    C = pd.read_csv(f"{OUT}/candidate_evidence.csv")
    p0 = np.load(f"{lib.AF_STORE}/p0_nonsnp.npy")
    M = np.load(f"{lib.GEA}/common/results/pool_matrices/pool_gen9_nonsnp_af.npy", mmap_mode="r")
    PM = pd.read_csv(ac.POOLMETA); site = PM.site.values
    A = ac.site_climate(); A.index = sorted(PM.site.astype(int).unique())
    print("high-leverage gardens per axis (|z|>2):")
    for ax in ["bio15", "pc1", "pc3", "bio6", "bio11", "bio10", "bio3", "bio8"]:
        print(f"  {ax:<6} {lever(A[ax])}")

    Gc = CE.site_means(M, p0, C.store_row.values, site); Gc.columns = C.index
    rho, rt, ntr = [], [], []
    for i, c in C.iterrows():
        ax = c.axis if c.axis in A.columns else c.best_axis_all
        a, b = stats(Gc[[i]], A[ax]); rho.append(a[0]); rt.append(b[0]); ntr.append(len(lever(A[ax])))
    C["rho_all"], C["r_trim"], C["n_trim"] = rho, rt, ntr

    rng = np.random.default_rng(2)
    B = pd.read_csv(f"{OUT}/repeat_variants_background_v2.csv.gz",
                    usecols=["chrom", "pos", "ref_len", "alt_len"])
    idx = np.load(f"{lib.AF_STORE}/index_nonsnp.npz", allow_pickle=True)
    K = pd.DataFrame({"chrom": idx["chrom"].astype(str), "pos": idx["pos"],
                      "ref_len": idx["ref_len"], "alt_len": idx["alt_len"]})
    K["row"] = np.arange(len(K))
    brow = np.sort(rng.choice(K.merge(B, on=["chrom", "pos", "ref_len", "alt_len"])["row"].values,
                              3000, replace=False))
    Gb = CE.site_means(M, p0, brow, site)
    # selection-aware: each background record on ITS best axis under each statistic
    R_rho = pd.DataFrame({ax: stats(Gb, A[ax])[0] for ax in A.columns})
    R_tr = pd.DataFrame({ax: stats(Gb, A[ax])[1] for ax in A.columns})
    bg_rho, bg_tr = np.sort(R_rho.max(1).values), np.sort(R_tr.max(1).values)
    C["pct_rho"] = 100 * np.searchsorted(bg_rho, C.rho_all.values) / len(bg_rho)
    C["pct_trim"] = 100 * np.searchsorted(bg_tr, C.r_trim.values) / len(bg_tr)
    C.to_csv(f"{OUT}/candidate_evidence.csv", index=False)
    print("GPX6:", C.loc[C.sym == "GPX6", ["r_own", "rho_all", "r_trim", "n_trim",
                                          "pct_sel", "pct_rho", "pct_trim"]].round(3).to_string(index=False))


if __name__ == "__main__":
    main()
