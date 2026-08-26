#!/usr/bin/env python
"""How much of the founder insertion-load x climate effect is provenance, really?

`_founder_climate_confound.py` controlled founder provenance with home bio1 + bio18 only,
and the insertion-fraction effect survived (rho -0.537 -> partial -0.363). But bio1+bio18
is a coarse two-number summary of where a founder is from, and founder home climate is the
single strongest predictor of the founder climate response (rho +0.632 with gamma_bio1).
If richer provenance control kills the residual, the SV reading dissolves into ordinary
local adaptation.

Three escalating provenance controls:
  (P1) home bio1 + bio18                 (what the previous script used)
  (P2) K home-climate PCs over all 19 bioclim axes (K = 2,3,5,8)
  (P3) all 19 home bioclim axes, raw

Plus a clade jackknife: drop each clade in turn and recompute, to check the effect is not
carried by one group of founders.

Env: kmate. Reads founder_climate_confound.npz + class_grms.npz + the 1001g bioclim table.
Writes results/sv_adaptive/founder_provenance_depth.csv + printed report.
"""
import os, sys
import numpy as np, pandas as pd
from scipy import stats
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import lib

GEA = lib.GEA
OUT = f"{GEA}/r1_sv_negative_selection/results/sv_adaptive"
ECO = ("/global/scratch/users/tbellg/gea_grene-net/key_files/"
       "1001g_regmap_grenet_ecotype_info_corrected_bioclim_2024May16.csv")
BIOS = [f"bio{i}" for i in range(1, 20)]


def partial_spearman(x, y, Z):
    rx = stats.rankdata(x); ry = stats.rankdata(y)
    if Z is None or Z.shape[1] == 0:
        r = float(stats.spearmanr(x, y).correlation)
        dof = len(rx) - 2
    else:
        RZ = np.column_stack([stats.rankdata(Z[:, j]) for j in range(Z.shape[1])])
        A = np.column_stack([np.ones(len(rx)), RZ])
        bx, *_ = np.linalg.lstsq(A, rx, rcond=None)
        by, *_ = np.linalg.lstsq(A, ry, rcond=None)
        ex = rx - A @ bx; ey = ry - A @ by
        r = float(np.corrcoef(ex, ey)[0, 1])
        dof = len(rx) - 2 - Z.shape[1]
    t = r * np.sqrt(dof / max(1e-12, 1 - r ** 2))
    return r, float(2 * stats.t.sf(abs(t), dof)), dof


def main():
    C = np.load(f"{OUT}/founder_climate_confound.npz", allow_pickle=True)
    G = np.load(f"{GEA}/r3_persite_gwas/results/varexp/class_grms.npz", allow_pickle=True)
    founders = C["founders"].astype("U6")
    keep = C["keep"].astype(bool)
    ins = C["ins_frac"][keep]
    g1 = C["gamma_bio1"][keep]; g18 = C["gamma_bio18"][keep]

    eco = pd.read_csv(ECO)
    eco["ecotypeid"] = eco["ecotypeid"].astype(str)
    em = eco.set_index("ecotypeid")
    H = np.column_stack([[float(em.loc[f, b]) if f in em.index else np.nan
                          for f in founders[keep]] for b in BIOS])
    assert np.isfinite(H).all(), "missing home bioclim"
    Hz = (H - H.mean(0)) / H.std(0)
    U, sv, _ = np.linalg.svd(Hz, full_matrices=False)
    varexp = sv ** 2 / (sv ** 2).sum()
    print(f"[home PCA] variance explained by first 8 PCs: "
          f"{np.round(varexp[:8], 3)}  (cum {varexp[:8].sum():.3f})")

    rows = []
    print("\n" + "=" * 78)
    print("P1-P3  insertion-fraction effect under escalating provenance control")
    for nm, gg in (("gamma_bio1", g1), ("gamma_bio18", g18)):
        print(f"\n  --- {nm} ---")
        specs = [("none", None),
                 ("home bio1+bio18", Hz[:, [0, 17]])]
        for k in (2, 3, 5, 8):
            specs.append((f"home PC1-{k}", U[:, :k] * sv[:k]))
        specs.append(("all 19 home bioclim", Hz))
        for lbl, Zc in specs:
            r, p, dof = partial_spearman(ins, gg, Zc)
            flag = "***" if p < 1e-3 else "**" if p < 0.01 else "*" if p < 0.05 else "ns"
            print(f"    partial(gamma, ins_frac | {lbl:22s}) = {r:+.3f}  p={p:.2e} {flag}")
            rows.append(dict(test="provenance", axis=nm, control=lbl, rho=r, p=p, dof=dof))

    # ---------------------------------------------------------------- clade jackknife
    print("\n" + "=" * 78)
    print("Clade jackknife (k=20): drop each clade, recompute rho(ins_frac, gamma_bio1)")
    K = G["K_snp"].astype(np.float64); K = K / np.mean(np.diag(K))
    Ks = K[np.ix_(keep, keep)]
    dist = np.add.outer(np.diag(Ks), np.diag(Ks)) - 2 * Ks
    np.fill_diagonal(dist, 0.0); dist = np.maximum(dist, 0.0)
    lab = fcluster(linkage(squareform(dist, checks=False), method="average"),
                   20, criterion="maxclust")
    full = stats.spearmanr(ins, g1).correlation
    jk = []
    for c in np.unique(lab):
        m = lab != c
        if m.sum() < 20:
            continue
        r = stats.spearmanr(ins[m], g1[m]).correlation
        jk.append((c, int((lab == c).sum()), r))
    jk.sort(key=lambda t: t[2])
    print(f"  full rho={full:+.3f}")
    print(f"  jackknife range: {jk[0][2]:+.3f} (drop clade {jk[0][0]}, n={jk[0][1]}) .. "
          f"{jk[-1][2]:+.3f} (drop clade {jk[-1][0]}, n={jk[-1][1]})")
    print(f"  all {len(jk)} jackknife estimates negative: {all(t[2] < 0 for t in jk)}")
    for c, sz, r in jk[:3]:
        print(f"    most influential: drop clade {c} (n={sz}) -> rho={r:+.3f}")
    rows.append(dict(test="jackknife", axis="gamma_bio1", control="drop-1-clade",
                     rho=full, p=np.nan, dof=len(jk)))

    pd.DataFrame(rows).to_csv(f"{OUT}/founder_provenance_depth.csv", index=False)
    print(f"\n[wrote] {OUT}/founder_provenance_depth.csv")


if __name__ == "__main__":
    main()
