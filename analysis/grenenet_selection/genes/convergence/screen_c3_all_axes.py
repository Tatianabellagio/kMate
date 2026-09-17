#!/usr/bin/env python
"""Re-test the climate criterion (C3) of `screen_3criteria.py` on EVERY climate axis, and
report it on the axis the GEA actually found each candidate on.

Why: the first screen tested per-garden allele-frequency change against bio1 and bio12
only. Of the 89 candidates, 3 were found by the GEA on bio1 and none on bio12; 39 were
found on bio15, 22 on pc3, 7 on bio3, 6 on bio10, 4 on bio6. The FLC-block lead
(AT5G10130) has p 1.6e-8 on bio6 and p 0.02 on bio12 -- judging it on bio12 could only
fail. So: same per-garden delta (gen9 mean - founding p0, 31 gardens), Spearman against
all 22 axes (bio1-19 + pc1-3 recomputed exactly as the GEA did, via axis_clusters).

Output -> results/screen_c3_all_axes.csv, one row per screen candidate:
  C3own_axis / C3own_rho / C3own_p      on `gea_best_axis` from the master table
  C3all_best_axis / _rho / _p           best over all 22
  C3all_same_cluster                    is the empirical best axis in the GEA axis's cluster
  rho_<axis>, p_<axis>                  all 22, so the table can be re-cut

Env: kmate. Compute node. ~2-8 min (Lustre read of the gen9 pool matrix dominates).
"""
from __future__ import annotations
import os
import sys
import numpy as np
import pandas as pd
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(HERE, "..", "..")))
import lib                                                       # noqa: E402
import axis_clusters as ac                                       # noqa: E402

PM = f"{lib.GEA}/common/results/pool_matrices"
OUT = f"{HERE}/results"
GEN = "gen9"


def garden_means(cols: np.ndarray) -> pd.DataFrame:
    M = np.load(f"{PM}/pool_{GEN}_nonsnp_af.npy", mmap_mode="r")
    meta = pd.read_csv(f"{PM}/pool_{GEN}_nonsnp.meta.csv")
    A = np.empty((M.shape[0], len(cols)), dtype=np.float32)
    for i in range(M.shape[0]):                     # row-at-a-time: see screen_3criteria
        A[i] = np.asarray(M[i])[cols]
    df = pd.DataFrame(A)
    df["site"] = meta.site.astype(int).values
    return df.groupby("site").mean()


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--all", action="store_true",
                    help="score EVERY record in screen_3criteria.csv (2,529), not just the 89 in "
                         "screen_top_loci.csv; output -> screen_c3_all_axes_full.csv")
    a = ap.parse_args()
    S = pd.read_csv(f"{OUT}/screen_3criteria.csv" if a.all else f"{OUT}/screen_top_loci.csv")
    suffix = "_full" if a.all else ""
    M = pd.read_csv(f"{OUT}/master_candidate_genes.csv").set_index("gene")
    S["gea_best_axis"] = S.target_gene.map(M.gea_best_axis)
    S["gea_sig_axes"] = S.target_gene.map(M.gea_sig_axes)

    # 31 sites x 22 axes, in sorted-site order (axis_clusters drops the index)
    A = ac.site_climate()
    sites_sorted = sorted(pd.read_csv(ac.POOLMETA).site.astype(int).unique())
    A.index = sites_sorted
    clus = ac.clusters()

    gm = garden_means(S.store_row.to_numpy())
    sites = gm.index.to_numpy()
    delta = gm.to_numpy(np.float32) - S.p0.to_numpy(np.float32)[None, :]
    Ax = A.loc[sites]

    rows = []
    for j, r in enumerate(S.itertuples()):
        d = delta[:, j]
        m = np.isfinite(d)
        rec = {"target_gene": r.target_gene, "chrom": r.chrom, "pos": r.pos,
               "store_row": r.store_row, "gea_best_axis": r.gea_best_axis,
               "gea_sig_axes": r.gea_sig_axes}
        best = (None, 0.0, 1.0)
        for ax in ac.AXES:
            res = stats.spearmanr(d[m], Ax[ax].to_numpy(float)[m])
            rec[f"rho_{ax}"], rec[f"p_{ax}"] = res.statistic, res.pvalue
            if abs(res.statistic) > abs(best[1]):
                best = (ax, res.statistic, res.pvalue)
        own = r.gea_best_axis if isinstance(r.gea_best_axis, str) else "bio1"
        rec["C3own_axis"] = own
        rec["C3own_rho"], rec["C3own_p"] = rec[f"rho_{own}"], rec[f"p_{own}"]
        rec["C3all_best_axis"], rec["C3all_best_rho"], rec["C3all_best_p"] = best
        rec["C3all_same_cluster"] = clus.get(best[0]) == clus.get(own)
        rec["C3_bio1_rho"], rec["C3_bio12_rho"] = rec["rho_bio1"], rec["rho_bio12"]
        rows.append(rec)
    R = pd.DataFrame(rows)
    R["PASS_C3own"] = R.C3own_p < 0.05
    # carry the original bio1/bio12 pass flags so the two C3 definitions can be compared
    for c in ("ftier", "PASS_C1", "PASS_C2", "PASS_C3", "n_pass", "p0", "C2_n_big",
              "C2_consistency", "key_ambiguous"):
        if c in S.columns:
            R[c] = S[c].to_numpy()
    R.to_csv(f"{OUT}/screen_c3_all_axes{suffix}.csv", index=False)
    print(f"{len(R)} candidates; own-axis C3 passes: {int(R.PASS_C3own.sum())}; "
          f"best axis in same cluster as GEA axis: {int(R.C3all_same_cluster.sum())}")
    if not a.all:
        print(R[["target_gene", "gea_best_axis", "C3own_rho", "C3own_p", "C3all_best_axis",
                 "C3all_best_rho", "C3_bio1_rho", "C3_bio12_rho"]].round(3).to_string())
    else:
        ok = R.PASS_C1 & R.PASS_C2
        print(f"C1+C2 pass: {int(ok.sum())} records, {R[ok].target_gene.nunique()} genes")
        print(f"  + C3 on bio1/bio12 (old): {int((ok & R.PASS_C3).sum())} records, "
              f"{R[ok & R.PASS_C3].target_gene.nunique()} genes")
        print(f"  + C3 on own GEA axis   : {int((ok & R.PASS_C3own).sum())} records, "
              f"{R[ok & R.PASS_C3own].target_gene.nunique()} genes")
        new = ok & R.PASS_C3own & ~R.PASS_C3
        print(f"  gained by own-axis test : {int(new.sum())} records, {R[new].target_gene.nunique()} genes")


if __name__ == "__main__":
    main()
