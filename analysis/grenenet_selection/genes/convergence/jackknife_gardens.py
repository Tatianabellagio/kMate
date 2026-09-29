#!/usr/bin/env python
"""Is a candidate's climate gradient carried by the gardens, or by one garden?

candidate_evidence.py aggregates to the 31 gardens, which removes pool pseudo-replication
but NOT a single-garden effect: a founder sweep in one garden (the garden-4 pattern the
review brief grades WEAK regardless of axis) produces a large garden-level r on whatever
axis that garden happens to be extreme on. So recompute the own-axis r 31 times, leaving
each garden out, and keep:

  r_jk_min     the smallest |r| over the 31 leave-one-out fits (worst case)
  jk_drop_g    the garden whose removal hurts most
  r_no_g4      |r| with garden 4 removed
  pct_jk       r_jk_min against the same selection-aware background, jackknifed the same
               way, so the comparison stays like-for-like

Writes the columns back into results/candidate_evidence.csv.
"""
import os, sys
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(HERE, "..", ".."))); sys.path.insert(0, HERE)
import lib, axis_clusters as ac                                  # noqa: E402
import candidate_evidence as CE                                  # noqa: E402

OUT = f"{HERE}/results"


def jk(G, a):
    """G: sites x vars site-mean table; a: site vector. Returns (min |r|, argmin site, per-site)."""
    sites = list(G.index); res = []
    for s in sites:
        keep = [x for x in sites if x != s]
        R = CE.corr_axes(G.loc[keep], a.loc[keep].to_frame("ax"))["ax"].abs().values
        res.append(R)
    res = np.array(res)                                  # sites x vars
    return res.min(0), np.array(sites)[res.argmin(0)], res


def main():
    C = pd.read_csv(f"{OUT}/candidate_evidence.csv")
    p0 = np.load(f"{lib.AF_STORE}/p0_nonsnp.npy")
    M = np.load(f"{lib.GEA}/common/results/pool_matrices/pool_gen9_nonsnp_af.npy", mmap_mode="r")
    PM = pd.read_csv(ac.POOLMETA); site = PM.site.values
    A = ac.site_climate(); A.index = sorted(PM.site.astype(int).unique())
    Gc = CE.site_means(M, p0, C.store_row.values, site)
    Gc.columns = C.index
    mins, drop, no4 = [], [], []
    for i, c in C.iterrows():
        ax = c.axis if c.axis in A.columns else c.best_axis_all
        m, d, res = jk(Gc[[i]], A[ax])
        mins.append(m[0]); drop.append(int(d[0]))
        no4.append(res[list(Gc.index).index(4), 0] if 4 in Gc.index else np.nan)
    C["r_jk_min"], C["jk_drop_g"], C["r_no_g4"] = mins, drop, no4

    # background, jackknifed the same way on each record's best axis
    rng = np.random.default_rng(1)
    B = pd.read_csv(f"{OUT}/repeat_variants_background_v2.csv.gz",
                    usecols=["chrom", "pos", "ref_len", "alt_len"])
    idx = np.load(f"{lib.AF_STORE}/index_nonsnp.npz", allow_pickle=True)
    K = pd.DataFrame({"chrom": idx["chrom"].astype(str), "pos": idx["pos"],
                      "ref_len": idx["ref_len"], "alt_len": idx["alt_len"]})
    K["row"] = np.arange(len(K))
    brow = np.sort(rng.choice(K.merge(B, on=["chrom", "pos", "ref_len", "alt_len"])["row"].values,
                              1500, replace=False))
    Gb = CE.site_means(M, p0, brow, site)
    best = CE.corr_axes(Gb, A).abs().idxmax(1)
    bj = []
    for ax, cols in best.groupby(best).groups.items():
        m, _, _ = jk(Gb[list(cols)], A[ax]); bj.extend(m)
    bj = np.sort(np.array(bj))
    C["pct_jk"] = 100 * np.searchsorted(bj, C.r_jk_min.values) / len(bj)
    C.to_csv(f"{OUT}/candidate_evidence.csv", index=False)
    print(f"background jackknifed best-axis |r|: median {np.median(bj):.3f}, p97 "
          f"{np.percentile(bj, 97):.3f}, p99 {np.percentile(bj, 99):.3f}")
    print("GPX6:", C.loc[C.sym == "GPX6", ["r_own", "r_jk_min", "jk_drop_g", "r_no_g4",
                                          "pct_sel", "pct_jk"]].round(3).to_string(index=False))


if __name__ == "__main__":
    main()
