#!/usr/bin/env python
"""Aggregate the multi-axis clq0.9 WZA outputs into one master summary table.

Scans multiaxis/wza/wza_{model}_{cls}_gen9_{axis}_deg2.csv for all 20 axes
(bio1..bio19 + pc1) x 3 models x 2 classes, computes per-output BH-FDR / Bonferroni
block counts + min block p, and the CAM5 (Chr2_4332) rank/q. Robust to missing
files (batch still running) — reports what's present.

Output: multiaxis/multiaxis_summary.csv  (one row per axis,model,cls)

  PY=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python
  $PY aggregate_multiaxis.py
"""
from __future__ import annotations
import os, sys, glob, re
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import lib

WD = f"{lib.GEA}/phase1_replication/multiaxis/wza"
AXES = [f"bio{i}" for i in range(1, 20)] + ["pc1"]
MODELS = ["kendall", "lfmm", "binomial"]
CLASSES = ["snp", "nonsnp"]
CAM5 = "Chr2_4332"
_RE = re.compile(r"wza_(?P<model>\w+?)_(?P<cls>snp|nonsnp)_gen9_(?P<axis>bio\d+|pc1)_deg2\.csv$")


def bh(p):
    p = np.asarray(p, float); n = len(p); o = np.argsort(p); q = np.empty(n)
    q[o] = (p[o] * n) / (np.arange(n) + 1)
    q[o] = np.minimum.accumulate(q[o][::-1])[::-1]
    return np.clip(q, 0, 1)


def main():
    rows = []
    for axis in AXES:
        for model in MODELS:
            for cls in CLASSES:
                f = f"{WD}/wza_{model}_{cls}_gen9_{axis}_deg2.csv"
                if not os.path.exists(f):
                    rows.append(dict(axis=axis, model=model, cls=cls, status="MISSING"))
                    continue
                w = pd.read_csv(f).rename(columns={"index": "block"})
                w["block"] = w["block"].astype(str)
                w = w[w["Z_pVal"].notna()].copy()
                q = bh(w["Z_pVal"].to_numpy())
                w["q"] = q
                n = len(w)
                cam = w[w["block"] == CAM5]
                cam_q = float(cam["q"].min()) if len(cam) else np.nan
                rows.append(dict(
                    axis=axis, model=model, cls=cls, status="ok", n_blocks=n,
                    n_bh=int((q < 0.05).sum()),
                    n_bonf=int((w["Z_pVal"] < 0.05 / n).sum()),
                    min_Zp=float(w["Z_pVal"].min()),
                    top_block=str(w.loc[w["Z_pVal"].idxmin(), "block"]),
                    cam5_q=cam_q))
    df = pd.DataFrame(rows)
    out = f"{lib.GEA}/phase1_replication/multiaxis/multiaxis_summary.csv"
    df.to_csv(out, index=False)
    ok = df[df.status == "ok"]
    print(f"present: {len(ok)}/{len(df)} outputs")
    if len(ok):
        piv = ok.pivot_table(index=["model", "cls"], columns="axis", values="n_bh")
        piv = piv.reindex(columns=[a for a in AXES if a in piv.columns])
        print("\nBH q<0.05 block counts (rows=model×class, cols=axis):")
        print(piv.to_string())
    print(f"\n-> {out}")


if __name__ == "__main__":
    main()
