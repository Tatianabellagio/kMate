#!/usr/bin/env python
"""Empirical correlation clusters of the 22 climate axes.

Why this exists
---------------
"Recurrent across multiple bioclim axes" is the headline triage signal for the GEA
candidate pool, but the bioclim axes are strongly correlated: bio1/5/9/10/11 are all
warm-season temperature, bio12/13/16/19 are all precipitation amount. A variant
significant on five temperature axes has ONE line of evidence counted five times, not
five. `n_axes` (the column the older `dissection/screen_sig_blocks.py` reports) has
exactly this problem.

So recurrence is counted over correlated CLUSTERS of axes, derived from the data
rather than assumed: average-linkage hierarchical clustering on 1 - |pearson r| of the
19 bioclim variables + pc1-3, at the SITE level (31 sites, one climate vector each --
pools within a site share a climate, so pooling would only reweight by plot count).

The PCs are recomputed here exactly as `build_class_matrices.py` does (SVD of the
z-scored 19 bioclim variables; pc1 oriented to +bio1, pc2/pc3 oriented so the largest
absolute loading is positive), so the cluster assignment matches the axes the LFMM
scans were actually run on.

env: kmate.
"""
from __future__ import annotations
import os
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform

HERE = os.path.dirname(os.path.abspath(__file__))
POOLMETA = os.path.abspath(os.path.join(
    HERE, "..", "..", "common", "results", "pool_matrices",
    "pool_gen9_nonsnp.meta.csv"))

BIOS = [f"bio{i}" for i in range(1, 20)]
AXES = BIOS + ["pc1", "pc2", "pc3"]

# cluster if mean |r| >= this. 0.6 gives 7 clusters and is the setting used for the
# candidate pool; 0.7 splits bio15 out of the dry-precipitation cluster (9 clusters).
R_CUT = 0.6


def site_climate() -> pd.DataFrame:
    """31 sites x 22 axes, PCs recomputed as in build_class_matrices.py."""
    d = pd.read_csv(POOLMETA)
    S = d.groupby("site")[BIOS].first()
    m = S.to_numpy(float)
    Z = (m - m.mean(0)) / m.std(0, ddof=0)
    U, Sv, Vt = np.linalg.svd(Z - Z.mean(0), full_matrices=False)
    pcs = {}
    for k in range(3):
        sc = U[:, k] * Sv[k]
        flip = (np.corrcoef(sc, m[:, 0])[0, 1] < 0 if k == 0
                else Vt[k][np.argmax(np.abs(Vt[k]))] < 0)
        if flip:
            sc = -sc
        pcs[f"pc{k + 1}"] = (sc - sc.mean()) / sc.std(ddof=0)
    return pd.concat([S.reset_index(drop=True), pd.DataFrame(pcs)], axis=1)[AXES]


def clusters(r_cut: float = R_CUT) -> dict[str, str]:
    """axis -> cluster label. Labels are the cluster's most-central axis."""
    A = site_climate()
    C = A.corr()
    D = (1.0 - C.abs()).to_numpy().copy()
    np.fill_diagonal(D, 0.0)
    D = (D + D.T) / 2.0
    lab = fcluster(linkage(squareform(D, checks=False), "average"),
                   1.0 - r_cut, "distance")
    out = {}
    for cid in np.unique(lab):
        members = [a for a, l in zip(C.columns, lab) if l == cid]
        # name the cluster after its most-central member (highest mean |r| within)
        sub = C.loc[members, members].abs()
        head = sub.mean(axis=1).idxmax() if len(members) > 1 else members[0]
        for a in members:
            out[a] = head
    return out


def main():
    cl = clusters()
    rows = [(a, cl[a]) for a in AXES]
    df = pd.DataFrame(rows, columns=["axis", "cluster"])
    os.makedirs(f"{HERE}/results", exist_ok=True)
    df.to_csv(f"{HERE}/results/axis_clusters.csv", index=False)
    print(f"{df.cluster.nunique()} clusters at |r| >= {R_CUT}")
    for c, g in df.groupby("cluster", sort=False):
        print(f"  {c:6s} <- {', '.join(g.axis)}")
    print(f"wrote {HERE}/results/axis_clusters.csv")


if __name__ == "__main__":
    main()
