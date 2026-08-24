#!/usr/bin/env python
"""Summarize the power experiments: for each LFMM output (raw + GIF p), report GIF,
min p, and Bonferroni / BH-FDR pass counts, so we can see which setup gains power."""
import sys, glob, os
import numpy as np, pandas as pd

def bh(p):
    p = np.asarray(p, float); n = len(p); o = np.argsort(p); q = np.empty(n)
    q[o] = (p[o] * n) / (np.arange(n) + 1); q[o] = np.minimum.accumulate(q[o][::-1])[::-1]
    return np.clip(q, 0, 1)

def line(tag, p, gif):
    p = np.asarray(p, float); p = p[np.isfinite(p)]; n = len(p)
    bonf = int((p < 0.05 / n).sum())
    q = bh(p); f05 = int((q < 0.05).sum()); f10 = int((q < 0.10).sum())
    print(f"  {tag:26s} GIF={gif:>5} minp={p.min():.1e}  Bonf={bonf:>4}  FDR.05={f05:>4}  FDR.10={f10:>4}")

def main():
    files = sorted(sys.argv[1:]) or sorted(glob.glob(
        "/global/scratch/users/tbellg/kmate/analysis/grenenet_gea/r2_gea_nonsnp/phase1_replication/results/**/*_bothp*.csv",
        recursive=True))
    for f in files:
        d = pd.read_csv(f)
        gtxt = f.replace(".csv", ".gif.txt")
        gif = open(gtxt).read().strip()[:5] if os.path.exists(gtxt) else "?"
        name = os.path.basename(f).replace("_bothp.csv", "")
        print(f"\n{name}")
        if "pval_raw" in d:
            line("raw (uncalibrated)", d.pval_raw, gif)
            line("GIF-calibrated", d.pval_gif, "1.00")
        else:
            line("p", d.iloc[:, 0], gif)

if __name__ == "__main__":
    main()
