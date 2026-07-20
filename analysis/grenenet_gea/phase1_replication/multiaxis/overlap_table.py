#!/usr/bin/env python
"""Cross-axis SNP vs non-SNP BH-hit-block OVERLAP table (all 20 axes).

For every axis (bio1..bio19 + pc1) and model (kendall/lfmm/binomial), take the
clq0.9 blocks with WZA BH-FDR q<0.05 in SNP and in non-SNP, and report the overlap:
  n_snp, n_nonsnp, n_shared, n_snp_only, n_nonsnp_only, jaccard.
Also a per-axis "any-model" view (block counted if BH-sig in >=1 model for that class).

Output:
  multiaxis/overlap_snp_nonsnp_by_axis_model.csv   (axis x model rows)
  multiaxis/overlap_snp_nonsnp_by_axis.csv         (axis rows, union-over-models)

  PY=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python
  $PY overlap_table.py
"""
from __future__ import annotations
import os, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import lib

WD = f"{lib.GEA}/phase1_replication/results/multiaxis/wza"
OUTDIR = f"{lib.GEA}/phase1_replication/results/multiaxis"
AXES = [f"bio{i}" for i in range(1, 20)] + ["pc1"]
MODELS = ["kendall", "lfmm", "binomial"]


def bh_sig_blocks(model, cls, axis, q=0.05):
    f = f"{WD}/wza_{model}_{cls}_gen9_{axis}_deg2.csv"
    if not os.path.exists(f):
        return None
    w = pd.read_csv(f).rename(columns={"index": "block"})
    w["block"] = w["block"].astype(str)
    w = w[w["Z_pVal"].notna()].copy()
    p = w["Z_pVal"].to_numpy(); n = len(p); o = np.argsort(p); qv = np.empty(n)
    qv[o] = (p[o] * n) / (np.arange(n) + 1); qv[o] = np.minimum.accumulate(qv[o][::-1])[::-1]
    return set(w.loc[np.clip(qv, 0, 1) < q, "block"])


def overlap(snp, non):
    sh = snp & non; un = snp | non
    return dict(n_snp=len(snp), n_nonsnp=len(non), n_shared=len(sh),
                n_snp_only=len(snp - non), n_nonsnp_only=len(non - snp),
                jaccard=round(len(sh) / len(un), 3) if un else np.nan)


def main():
    rows, per_axis = [], {}
    for axis in AXES:
        union = {"snp": set(), "nonsnp": set()}
        any_present = False
        for model in MODELS:
            s = bh_sig_blocks(model, "snp", axis); n = bh_sig_blocks(model, "nonsnp", axis)
            if s is None or n is None:
                continue
            any_present = True
            union["snp"] |= s; union["nonsnp"] |= n
            rows.append(dict(axis=axis, model=model, **overlap(s, n)))
        if any_present:
            per_axis[axis] = dict(axis=axis, scope="any_model", **overlap(union["snp"], union["nonsnp"]))

    by_am = pd.DataFrame(rows)
    by_ax = pd.DataFrame(per_axis.values())
    by_am.to_csv(f"{OUTDIR}/overlap_snp_nonsnp_by_axis_model.csv", index=False)
    by_ax.to_csv(f"{OUTDIR}/overlap_snp_nonsnp_by_axis.csv", index=False)
    print("=== per axis x model ===")
    print(by_am.to_string(index=False))
    print("\n=== per axis (union over models) ===")
    print(by_ax.to_string(index=False))
    print(f"\n-> {OUTDIR}/overlap_snp_nonsnp_by_axis_model.csv")
    print(f"-> {OUTDIR}/overlap_snp_nonsnp_by_axis.csv")


if __name__ == "__main__":
    main()
