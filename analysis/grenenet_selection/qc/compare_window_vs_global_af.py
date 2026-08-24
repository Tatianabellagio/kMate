#!/usr/bin/env python
"""Per-sample comparison of kMate per-variant AF: window-mode vs global-mode.

Both runs are row-aligned (identical panel, same variant order) — verified by pos match.
For one sample we read the alt_freq (and n_called) column from each mode across all 5 chroms,
compute delta = window - global, and emit one row of aggregate stats + a |delta| histogram.

Usage:
  python compare_window_vs_global_af.py SAMPLE GLOBAL_DIR WINDOW_DIR OUT_CSV
"""
import sys, os
import numpy as np
import pandas as pd

SAMPLE, GDIR, WDIR, OUT = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]
CHROMS = [f"Chr{i}" for i in range(1, 6)]

# |delta| histogram bin edges (absolute AF difference)
EDGES = np.array([0, 1e-4, 1e-3, 5e-3, 1e-2, 2e-2, 5e-2, 1e-1, 2e-1, 5e-1, 1.0001])

def load_af(d, c):
    f = os.path.join(d, f"{SAMPLE}_{c}.tsv")
    df = pd.read_csv(f, sep="\t", usecols=["alt_freq", "n_called"],
                     dtype={"alt_freq": "float32", "n_called": "int32"})
    return df["alt_freq"].to_numpy(), df["n_called"].to_numpy()

n = 0
sad = 0.0          # sum |d|
ssd = 0.0          # sum d^2
maxad = 0.0
# pearson accumulators
sx = sy = sxx = syy = sxy = 0.0
# threshold exceedance counts
thr = {0.0: 0, 1e-3: 0, 1e-2: 0, 5e-2: 0, 1e-1: 0}
hist = np.zeros(len(EDGES) - 1, dtype=np.int64)
n_finite_both = 0

for c in CHROMS:
    g, gnc = load_af(GDIR, c)
    w, wnc = load_af(WDIR, c)
    if g.shape != w.shape:
        sys.exit(f"shape mismatch {c}: {g.shape} vs {w.shape}")
    m = np.isfinite(g) & np.isfinite(w)
    g, w = g[m].astype(np.float64), w[m].astype(np.float64)
    d = w - g
    ad = np.abs(d)
    n += d.size
    n_finite_both += int(m.sum())
    sad += ad.sum()
    ssd += np.dot(d, d)
    if ad.size:
        maxad = max(maxad, float(ad.max()))
    sx += g.sum(); sy += w.sum(); sxx += np.dot(g, g); syy += np.dot(w, w); sxy += np.dot(g, w)
    for t in thr:
        thr[t] += int((ad > t).sum())
    hist += np.histogram(ad, bins=EDGES)[0]

mean_ad = sad / n if n else np.nan
rms = np.sqrt(ssd / n) if n else np.nan
# pearson r
cov = sxy - sx * sy / n
vx = sxx - sx * sx / n
vy = syy - sy * sy / n
r = cov / np.sqrt(vx * vy) if vx > 0 and vy > 0 else np.nan

row = {
    "sample": SAMPLE,
    "n_variants": n,
    "mean_abs_delta": mean_ad,
    "rms_delta": rms,
    "max_abs_delta": maxad,
    "pearson_r": r,
    "frac_changed": thr[0.0] / n if n else np.nan,
    "frac_delta_gt_1e-3": thr[1e-3] / n if n else np.nan,
    "frac_delta_gt_1e-2": thr[1e-2] / n if n else np.nan,
    "frac_delta_gt_5e-2": thr[5e-2] / n if n else np.nan,
    "frac_delta_gt_1e-1": thr[1e-1] / n if n else np.nan,
}
for i in range(len(hist)):
    row[f"hist_{EDGES[i]:.4g}_{EDGES[i+1]:.4g}"] = int(hist[i])

pd.DataFrame([row]).to_csv(OUT, index=False)
print(f"{SAMPLE}: n={n:,} mean|d|={mean_ad:.5f} rms={rms:.5f} r={r:.5f} "
      f"frac>0.01={row['frac_delta_gt_1e-2']:.4f} frac>0.05={row['frac_delta_gt_5e-2']:.4f}")
