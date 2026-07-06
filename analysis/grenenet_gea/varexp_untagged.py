#!/usr/bin/env python
"""Does the SNP-UNTAGGED non-SNP layer explain selection-trait variance beyond K_snp?

The direct test of the VAREXP_SELECTION_HANDOFF.md caveat: the genome-wide K_nonsnp is
redundant with K_snp (corr=0.998) because it's dominated by markers SNPs already tag via LD --
that answers "is non-SNP relatedness different" (no, by construction), not "what does kMate's
non-SNP layer add that SNPs can't see". Here K_nonsnp is replaced by K_untagged_{r02,r05}
(build_untagged_grm.py: indel+SV markers with founder-LD r2 < 0.2 / 0.5 against every panel SNP
within +/-50kb) and re-run through the same three estimators as varexp_selection.py (marginal
h2, joint 2-GRM REML + LRT, GBLUP CV gain), reusing its functions directly.

Reads results/grenenet_gea/varexp/{selection_s_matrix.npz, class_grms.npz, untagged_grms.npz}.
Writes varexp_untagged.csv. Env: kmate.
"""
from __future__ import annotations
import os, sys, json
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lib
from varexp_selection import h2_1grm, joint_2grm, gblup_cv, rint

OUT = f"{lib.GEA}/varexp"
THRESHOLDS = ["r02", "r05"]


def main():
    S = np.load(f"{OUT}/selection_s_matrix.npz", allow_pickle=True)
    G = np.load(f"{OUT}/class_grms.npz", allow_pickle=True)
    U = np.load(f"{OUT}/untagged_grms.npz", allow_pickle=True)
    nmark = json.loads(str(U["n_markers"]))
    tf = S["founders"].astype(str)
    gf = G["founders"].astype(str)
    uf = U["founders"].astype(str)
    assert (gf == uf).all(), "class_grms / untagged_grms founder order mismatch"
    ana = S["analyzable"].astype(bool)
    gpos = {f: i for i, f in enumerate(gf)}
    keep = np.array([ana[i] and tf[i] in gpos for i in range(len(tf))])
    order = [gpos[tf[i]] for i in range(len(tf)) if keep[i]]
    n = len(order)
    K_snp = G["K_snp"][np.ix_(order, order)]
    K_u = {k: U[f"K_untagged_{k}"][np.ix_(order, order)] for k in THRESHOLDS}
    Smat = S["S"][:, keep]
    sites = S["sites"]; bio1 = S["bio1"]
    print(f"founders analyzable & genotyped: {n}   sites: {len(sites)}", flush=True)
    print("untagged marker counts:", nmark, flush=True)
    for k in THRESHOLDS:
        c = np.corrcoef(K_snp[np.triu_indices(n, 1)], K_u[k][np.triu_indices(n, 1)])[0, 1]
        print(f"corr(K_snp, K_untagged_{k}) offdiag = {c:.3f}", flush=True)

    zone = pd.qcut(bio1, 3, labels=["cold", "mid", "hot"]).astype(str)
    axes = {"w_global": Smat.mean(0)}
    for z in ("cold", "mid", "hot"):
        axes[f"w_{z}"] = Smat[zone == z].mean(0)

    rows = []
    for name, y_raw in axes.items():
        for transform, y in (("raw", y_raw), ("rint", rint(y_raw))):
            y = np.asarray(y, float)
            for k in THRESHOLDS:
                h2, se, p = h2_1grm(y, K_u[k])
                rows.append(dict(trait=name, transform=transform, thresh=k, method="marginal",
                                 cls="untagged", value=h2, se=se, p=p))
                frac, p_add = joint_2grm(y, K_snp, K_u[k])
                for cl, v in zip(["snp", "untagged", "resid"], frac):
                    rows.append(dict(trait=name, transform=transform, thresh=k, method="joint2",
                                     cls=cl, value=v, se=np.nan,
                                     p=(p_add if cl == "untagged" else np.nan)))
                r_s, r2_s = gblup_cv(y, [K_snp])
                r_u, r2_u = gblup_cv(y, [K_u[k]])
                r_b, r2_b = gblup_cv(y, [K_snp, K_u[k]])
                for cl, rr, r2 in [("snp", r_s, r2_s), ("untagged", r_u, r2_u), ("both", r_b, r2_b)]:
                    rows.append(dict(trait=name, transform=transform, thresh=k, method="predcv",
                                     cls=cl, value=r2, se=np.nan, p=rr))
                rows.append(dict(trait=name, transform=transform, thresh=k, method="predcv",
                                 cls="gain_untagged", value=r2_b - r2_s, se=np.nan, p=r_b - r_s))
            print(f"{name}/{transform} done", flush=True)

    df = pd.DataFrame(rows)
    df.to_csv(f"{OUT}/varexp_untagged.csv", index=False)

    print("\n===== SUMMARY (w_global, raw) =====")
    wg = df[(df["trait"] == "w_global") & (df["transform"] == "raw")]
    for k in THRESHOLDS:
        sub = wg[wg["thresh"] == k]
        m = {r.cls: round(r.value, 4) for r in sub[sub["method"] == "marginal"].itertuples()}
        j = {r.cls: round(r.value, 4) for r in sub[sub["method"] == "joint2"].itertuples()}
        pcv = {r.cls: round(r.value, 4) for r in sub[sub["method"] == "predcv"].itertuples()}
        p_add = sub[(sub["method"] == "joint2") & (sub["cls"] == "untagged")]["p"].iloc[0]
        print(f"[{k}, n_markers={nmark[k]['total']}] marginal_h2(untagged)={m.get('untagged')} "
              f"joint_frac={j} p_add={p_add:.3g}  predR2={pcv}")
    print(f"wrote {OUT}/varexp_untagged.csv")


if __name__ == "__main__":
    main()
