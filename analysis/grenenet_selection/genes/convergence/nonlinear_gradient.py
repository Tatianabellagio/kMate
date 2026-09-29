#!/usr/bin/env python
"""Climate response that is not a slope: humps and U-shapes, selection-aware.

A linear garden-level correlation assumes the allele is favoured at one climatic extreme.
An allele favoured at intermediate (temperate) sites rises in the middle gardens and falls
at both ends -- a hump with ~zero linear slope, invisible to the linear test (the CRY2-side
1,634 bp insertion shows the reverse shape in founder origins). So per variant and per axis:

  fit   delta_garden ~ z + z^2   (z = standardised axis value across the 31 gardens)
  R2_quad   variance explained by the quadratic model
  curvature sign of the z^2 coefficient (negative = hump / intermediate optimum,
            positive = U / favoured at both extremes)
  shape     'linear' if the linear term dominates, else 'hump' or 'U'

Best axis by R2_quad, compared with background records' best-of-22 R2_quad
(selection-aware), and also on the genome-wide-residual profile (locus-specific).

Adds columns to results/evidence_matrix_A.csv: quad_axis, r2_quad, curvature, shape,
pct_quad, pct_quad_resid.
"""
import os, sys
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(HERE, "..", ".."))); sys.path.insert(0, HERE)
import lib, axis_clusters as ac                                  # noqa: E402
import candidate_evidence as CE                                  # noqa: E402

OUT = f"{HERE}/results"
K4 = ["chrom", "pos", "ref_len", "alt_len"]


def quad(G, A):
    """G sites x vars; A sites x axes -> per var: best axis, R2, curvature, linear share."""
    A = A.reindex(G.index)
    Y = G.values - np.nanmean(G.values, 0)
    Y = np.nan_to_num(Y)
    best = dict(ax=None, r2=np.full(Y.shape[1], -1.0), curv=np.zeros(Y.shape[1]),
                lin=np.zeros(Y.shape[1]))
    axes = np.empty(Y.shape[1], dtype=object)
    tss = (Y ** 2).sum(0) + 1e-12
    for ax in A.columns:
        z = (A[ax] - A[ax].mean()) / A[ax].std()
        X = np.column_stack([np.ones(len(z)), z, z ** 2 - (z ** 2).mean()])
        beta, *_ = np.linalg.lstsq(X, Y, rcond=None)
        r2 = 1 - ((Y - X @ beta) ** 2).sum(0) / tss
        lin_ss = (np.outer(z, beta[1]) ** 2).sum(0)
        quad_ss = (np.outer(X[:, 2], beta[2]) ** 2).sum(0)
        m = r2 > best["r2"]
        best["r2"][m] = r2[m]; best["curv"][m] = beta[2][m]
        best["lin"][m] = lin_ss[m] / np.maximum(lin_ss[m] + quad_ss[m], 1e-12)
        axes[m] = ax
    return axes, best["r2"], best["curv"], best["lin"]


def main():
    V = pd.read_csv(f"{OUT}/evidence_matrix_A.csv")
    PM = pd.read_csv(ac.POOLMETA); site = PM.site.values
    A = ac.site_climate(); A.index = sorted(PM.site.astype(int).unique())
    p0 = np.load(f"{lib.AF_STORE}/p0_nonsnp.npy")
    M = np.load(f"{lib.GEA}/common/results/pool_matrices/pool_gen9_nonsnp_af.npy", mmap_mode="r")
    idx = np.load(f"{lib.AF_STORE}/index_nonsnp.npz", allow_pickle=True)
    K = pd.DataFrame({"chrom": idx["chrom"].astype(str), "pos": idx["pos"],
                      "ref_len": idx["ref_len"], "alt_len": idx["alt_len"]})
    K["row"] = np.arange(len(K))
    rng = np.random.default_rng(12)
    Bg = pd.read_csv(f"{OUT}/repeat_variants_background_v2.csv.gz", usecols=K4)
    brow = np.sort(rng.choice(K.merge(Bg, on=K4)["row"].values, 4000, replace=False))
    z = lambda G: ((G - G.mean()) / G.std()).fillna(0)
    Gb = CE.site_means(M, p0, brow, site)
    U, S, _ = np.linalg.svd(z(Gb).values, full_matrices=False)
    P = U[:, :3]; H = P @ np.linalg.pinv(P)
    resid = lambda G: pd.DataFrame(z(G).values - H @ z(G).values, index=G.index, columns=G.columns)
    _, bg_r2, _, _ = quad(Gb, A); bg_r2 = np.sort(bg_r2)
    _, bg_r2r, _, _ = quad(resid(Gb), A); bg_r2r = np.sort(bg_r2r)

    Gc = CE.site_means(M, p0, V.store_row.values, site); Gc.columns = V.index
    ax, r2, curv, lin = quad(Gc, A)
    _, r2r, _, _ = quad(resid(Gc), A)
    V["quad_axis"], V["r2_quad"], V["curvature"] = ax, r2, curv
    V["shape"] = np.where(lin > 0.5, "linear", np.where(curv < 0, "hump", "U"))
    V["pct_quad"] = 100 * np.searchsorted(bg_r2, r2) / len(bg_r2)
    V["pct_quad_resid"] = 100 * np.searchsorted(bg_r2r, r2r) / len(bg_r2r)
    V.to_csv(f"{OUT}/evidence_matrix_A.csv", index=False)
    s = V[V.pct_quad >= 97]
    print(f"strong quadratic (>=97th, selection-aware): {len(s)}; shapes {s['shape'].value_counts().to_dict()}")
    print(f"  of which NOT caught by the linear test (pct_sel < 97): "
          f"{int((s.pct_sel < 97).sum())}  -> shapes {s[s.pct_sel < 97]['shape'].value_counts().to_dict()}")


if __name__ == "__main__":
    main()
