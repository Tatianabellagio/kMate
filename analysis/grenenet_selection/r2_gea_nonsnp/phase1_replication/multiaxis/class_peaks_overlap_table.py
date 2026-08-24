#!/usr/bin/env python
"""Before/after-WZA peak counts AND pairwise class overlap, all 4 classes, per model.

Extends snp_nonsnp_peaks_table.py's before/after count (rec_* = per-record BEFORE
WZA, blk_* = per-block AFTER WZA isotonic) to cover snp/nonsnp/sv/smallindel in one
table, then adds the piece that table doesn't have: pairwise overlap of AFTER-WZA
Bonferroni-significant blocks between classes.

Overlap unit = (axis, block) pairs, pooled over the 20 axes, per model. Two classes
"share" a hit when the same block is Bonferroni-sig in both classes at the same axis.
This matches how class_specific() in _build_snp_vs_nonsnp_viz_nb.py already counts
"new" (class-only) peaks -- this script adds the reciprocal (shared) side and puts
both in one table instead of only printing class-specific unions.
"""
from __future__ import annotations
import os, sys, itertools
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))))
import lib

MA = f"{lib.GEA}/phase1_replication/results/multiaxis"
WD = f"{MA}/wza"
AXES = [f"bio{i}" for i in range(1, 20)] + ["pc1"]
MODELS = ["kendall", "lfmm", "binomial"]
CLASSES = ["snp", "nonsnp", "sv", "smallindel"]


def bh(p):
    p = np.asarray(p, float); n = len(p)
    if n == 0: return p
    o = np.argsort(p); q = np.empty(n)
    q[o] = (p[o] * n) / (np.arange(n) + 1)
    q[o] = np.minimum.accumulate(q[o][::-1])[::-1]
    return np.clip(q, 0, 1)


def counts(p):
    p = np.asarray(p, float); p = p[np.isfinite(p)]
    n = len(p)
    if n == 0: return 0, 0, 0
    nb = int((p < 0.05 / n).sum())
    nf = int((bh(p) < 0.05).sum())
    return n, nb, nf


def load_wza(model, cls, axis):
    f = f"{WD}/wza_{model}_{cls}_gen9_{axis}_isotonic.csv"
    if not os.path.exists(f): return None
    w = pd.read_csv(f).rename(columns={"index": "block"})
    w["block"] = w["block"].astype(str)
    w = w[w["Z_pVal"].notna() & (w["Z_pVal"] > 0)].copy()
    return w


def bonf_sig_blocks(w):
    n = len(w)
    if n == 0: return set()
    return set(w.loc[w["Z_pVal"] < 0.05 / n, "block"])


def main():
    # --- part 1: before/after WZA counts, all 4 classes, per axis x model ---
    rows = []
    for axis in AXES:
        for model in MODELS:
            for cls in CLASSES:
                fin = f"{MA}/wza_in/{model}_{cls}_gen9_{axis}.csv"
                fout = f"{WD}/wza_{model}_{cls}_gen9_{axis}_isotonic.csv"
                if not (os.path.exists(fin) and os.path.exists(fout)):
                    rows.append(dict(axis=axis, model=model, cls=cls, status="MISSING")); continue
                rn, rb, rf = counts(pd.read_csv(fin)["pval"].to_numpy())
                w = pd.read_csv(fout)
                bn, bb, bf = counts(w["Z_pVal"].to_numpy())
                rows.append(dict(axis=axis, model=model, cls=cls, status="ok",
                                  rec_n=rn, rec_bonf=rb, rec_fdr=rf,
                                  blk_n=bn, blk_bonf=bb, blk_fdr=bf))
    df = pd.DataFrame(rows)
    df.to_csv(f"{MA}/class_peaks_beforeafter_4class.csv", index=False)
    ok = df[df.status == "ok"].copy()
    print(f"present: {len(ok)}/{len(df)} (axis x model x class) cells\n")

    agg = ok.groupby(["model", "cls"])[["rec_bonf", "rec_fdr", "blk_bonf", "blk_fdr"]].sum()
    agg = agg.reindex(pd.MultiIndex.from_product([MODELS, CLASSES], names=["model", "cls"]))
    print("=== TOTAL significant over 20 axes (bio1-19 + pc1), per model x class ===")
    print("   before-WZA = per-record ; after-WZA = per-block (isotonic)\n")
    hdr = f"{'model':9s} {'class':11s} | {'rec_Bonf':>9} {'rec_FDR':>8} | {'blk_Bonf':>9} {'blk_FDR':>8}"
    print(hdr); print("-" * len(hdr))
    for (m, c), r in agg.iterrows():
        print(f"{m:9s} {c:11s} | {int(r.rec_bonf):9d} {int(r.rec_fdr):8d} | {int(r.blk_bonf):9d} {int(r.blk_fdr):8d}")

    # --- part 2: pairwise class overlap of AFTER-WZA Bonferroni-sig (axis,block) hits ---
    sig = {}  # (model, cls) -> set of (axis, block)
    for model in MODELS:
        for cls in CLASSES:
            s = set()
            for axis in AXES:
                w = load_wza(model, cls, axis)
                if w is None: continue
                s |= {(axis, b) for b in bonf_sig_blocks(w)}
            sig[(model, cls)] = s

    orows = []
    for model in MODELS:
        for a, b in itertools.combinations(CLASSES, 2):
            sa, sb = sig[(model, a)], sig[(model, b)]
            shared = sa & sb
            union = sa | sb
            jac = len(shared) / len(union) if union else np.nan
            orows.append(dict(model=model, class_a=a, class_b=b,
                               n_a=len(sa), n_b=len(sb), n_shared=len(shared),
                               n_a_only=len(sa - sb), n_b_only=len(sb - sa),
                               jaccard=round(jac, 3) if union else np.nan))
    odf = pd.DataFrame(orows)
    odf.to_csv(f"{MA}/class_peaks_overlap_bonf.csv", index=False)
    print("\n=== Pairwise class overlap, Bonferroni-sig (axis,block) hits pooled over 20 axes ===")
    print("   n_a/n_b = total sig hits per class ; n_shared = same block sig in both classes at same axis\n")
    hdr2 = f"{'model':9s} {'A':11s} {'B':11s} | {'n_A':>5} {'n_B':>5} {'shared':>7} {'A_only':>7} {'B_only':>7} {'jaccard':>8}"
    print(hdr2); print("-" * len(hdr2))
    for _, r in odf.iterrows():
        print(f"{r.model:9s} {r.class_a:11s} {r.class_b:11s} | {r.n_a:5d} {r.n_b:5d} {r.n_shared:7d} "
              f"{r.n_a_only:7d} {r.n_b_only:7d} {r.jaccard:8.3f}")

    print(f"\n-> {MA}/class_peaks_beforeafter_4class.csv")
    print(f"-> {MA}/class_peaks_overlap_bonf.csv")


if __name__ == "__main__":
    main()
