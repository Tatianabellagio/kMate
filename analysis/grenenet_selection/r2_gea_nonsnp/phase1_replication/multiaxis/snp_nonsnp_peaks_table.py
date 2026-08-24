#!/usr/bin/env python
"""SNP vs non-SNP peaks — before/after WZA significance table, per model, all axes.

For every axis (bio1..bio19 + pc1) x model (kendall/lfmm/binomial) x class, count
significant tests BEFORE WZA (per-record, the wza_in `pval` that WZA aggregates) and
AFTER WZA (per-block, the isotonic `Z_pVal`), each at Bonferroni (p < 0.05/n) and
BH-FDR (q < 0.05). Same MAF>0.05 in-block record set feeds both, so the contrast
isolates block-aggregation.

--classes controls the class set: default the clean 3-class (snp/sv/smallindel).
Pass --classes snp nonsnp to build the pooled view (needs the nonsnp re-run present).

Outputs a tidy master CSV + prints per-(model,class) axis pivots for blk_fdr/blk_bonf.
"""
from __future__ import annotations
import argparse, os, sys, glob
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))))
import lib

MA = f"{lib.GEA}/phase1_replication/results/multiaxis"
AXES = [f"bio{i}" for i in range(1, 20)] + ["pc1"]
MODELS = ["kendall", "lfmm", "binomial"]


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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--classes", nargs="+", default=["snp", "sv", "smallindel"])
    ap.add_argument("--regime", default="isotonic")
    ap.add_argument("--out", default=f"{MA}/snp_nonsnp_peaks_beforeafter_3class.csv")
    args = ap.parse_args()

    rows = []
    for axis in AXES:
        for model in MODELS:
            for cls in args.classes:
                fin = f"{MA}/wza_in/{model}_{cls}_gen9_{axis}.csv"          # before (per-record)
                fout = f"{MA}/wza/wza_{model}_{cls}_gen9_{axis}_{args.regime}.csv"  # after (per-block)
                if not (os.path.exists(fin) and os.path.exists(fout)):
                    rows.append(dict(axis=axis, model=model, cls=cls, status="MISSING")); continue
                rn, rb, rf = counts(pd.read_csv(fin)["pval"].to_numpy())
                w = pd.read_csv(fout)
                bn, bb, bf = counts(w["Z_pVal"].to_numpy())
                rows.append(dict(axis=axis, model=model, cls=cls, status="ok",
                                 rec_n=rn, rec_bonf=rb, rec_fdr=rf,
                                 blk_n=bn, blk_bonf=bb, blk_fdr=bf))
    df = pd.DataFrame(rows)
    df.to_csv(args.out, index=False)
    ok = df[df.status == "ok"].copy()
    print(f"present: {len(ok)}/{len(df)} (axis x model x class) cells | classes={args.classes}\n")

    # headline: totals summed over the 20 axes, per (model, class)
    agg = ok.groupby(["model", "cls"])[["rec_bonf", "rec_fdr", "blk_bonf", "blk_fdr"]].sum()
    print("=== TOTAL significant over 20 axes (bio1-19 + pc1), per model x class ===")
    print("   before-WZA = per-record ; after-WZA = per-block\n")
    hdr = f"{'model':9s} {'class':11s} | {'rec_Bonf':>9} {'rec_FDR':>8} | {'blk_Bonf':>9} {'blk_FDR':>8}"
    print(hdr); print("-" * len(hdr))
    for (m, c), r in agg.iterrows():
        print(f"{m:9s} {c:11s} | {int(r.rec_bonf):9d} {int(r.rec_fdr):8d} | {int(r.blk_bonf):9d} {int(r.blk_fdr):8d}")

    # per-axis pivot of the after-WZA FDR block counts (the usual headline metric)
    for metric, lab in [("blk_fdr", "after-WZA blocks BH q<0.05"), ("blk_bonf", "after-WZA blocks Bonferroni")]:
        print(f"\n=== {lab}: rows=model x class, cols=axis ===")
        piv = ok.pivot_table(index=["model", "cls"], columns="axis", values=metric, fill_value=0)
        piv = piv.reindex(columns=[a for a in AXES if a in piv.columns])
        print(piv.astype(int).to_string())
    print(f"\n-> {args.out}")


if __name__ == "__main__":
    main()
