#!/usr/bin/env python
"""Compare LFMM per-record de-inflation schemes on one class/axis.

Schemes:
  raw            : lfmm raw two-sided p (no calibration)
  gif            : genomic control, divide z^2 by lambda=median(z^2)/0.4549 (lfmm built-in)
  median+MAD     : Efron simple empirical null, z' = (z-median)/(1.4826*MAD)
  efron_central  : Efron central-matching empirical null -- fit N(d0,s0) to the log-histogram
                   near the MODE (quadratic fit), so non-null shoulders don't inflate s0.
                   This is the one that can DIFFER from gif/MAD: if the bulk over-dispersion
                   is real polygenic signal (shoulders), s0 -> ~1 and deflation is mild.

For each: report the post-calibration lambda, #Bonferroni hits, #BH-FDR(q<0.05) hits,
and whether the strongest raw peaks survive. Usage: compare_calibrations.py <scores_csv> <class> <axis>
"""
import sys, os
import numpy as np, pandas as pd
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
GEA = os.path.dirname(os.path.dirname(HERE))
CM = f"{GEA}/phase1_replication/results/class_matrices"
WZAIN = f"{GEA}/phase1_replication/results/multiaxis/wza_in_clq09_tile"

scores_csv, CLS, AXIS = sys.argv[1], sys.argv[2], sys.argv[3]

sc = pd.read_csv(scores_csv)                       # score, praw, pgif  (record order = class matrix)
z = sc["score"].to_numpy(float)
n = len(z)

# attach block/chrom/pos via the wza_in file (same record set; align by raw p, which is deterministic)
wz = pd.read_csv(f"{WZAIN}/lfmm_{CLS}_gen9_{AXIS}.csv")
wz = wz[wz["MAF"] > 0.05].reset_index(drop=True)
assert len(wz) == n, (len(wz), n)
# sanity: my rerun raw p must match wza_in raw p to confirm row alignment
align = np.corrcoef(-np.log10(sc["praw"].clip(1e-300)), -np.log10(wz["pval"].clip(1e-300)))[0, 1]
print(f"row-alignment check (raw nlp corr, must be ~1.0): {align:.4f}")
df = wz[["chrom", "pos", "block", "MAF"]].copy()
df["z"] = z


def lam(p):
    p = np.clip(np.asarray(p, float), 1e-300, 1.0)
    return float(np.median(stats.chi2.isf(p, 1)) / stats.chi2.isf(0.5, 1))


def bh_hits(p, q=0.05):
    p = np.asarray(p, float); m = len(p); o = np.argsort(p)
    thr = q * (np.arange(1, m + 1) / m)
    passed = p[o] <= thr
    k = np.where(passed)[0]
    return int(k.max() + 1) if len(k) else 0


def efron_central(z):
    """Efron central-matching: quadratic fit to log density near the mode."""
    lo, hi = np.percentile(z, [1, 99])
    bins = np.linspace(lo, hi, 120)
    c, edges = np.histogram(z, bins=bins)
    x = 0.5 * (edges[:-1] + edges[1:])
    m = c > 0
    x, y = x[m], np.log(c[m])
    # central window: within +-2 of the empirical median
    med = np.median(z)
    w = np.abs(x - med) < 2.0
    b2, b1, b0 = np.polyfit(x[w], y[w], 2)
    s0 = np.sqrt(-1.0 / (2 * b2)) if b2 < 0 else np.nan
    d0 = -b1 / (2 * b2) if b2 < 0 else med
    return d0, s0


# ---- build calibrated p for each scheme ----
med, mad = np.median(z), stats.median_abs_deviation(z)
sig_mad = 1.4826 * mad
d0, s0 = efron_central(z)
lam_gif = lam(sc["praw"])

schemes = {
    "raw":           sc["praw"].to_numpy(float),
    "gif":           sc["pgif"].to_numpy(float),
    "median+MAD":    2 * stats.norm.sf(np.abs((z - med) / sig_mad)),
    "efron_central": 2 * stats.norm.sf(np.abs((z - d0) / s0)),
}

print(f"\nclass={CLS} axis={AXIS}  n={n:,}  blocks={df.block.nunique():,}")
print(f"z: median {med:+.3f}  sd {z.std():.3f}  |  lambda(raw)={lam_gif:.3f}")
print(f"null-scale estimates:  MAD-> s0={sig_mad:.3f}   Efron-central-> d0={d0:+.3f} s0={s0:.3f}")
bonf = 0.05 / n
print(f"\nper-record Bonferroni p = 0.05/{n} = {bonf:.2e}  (-log10={-np.log10(bonf):.2f})")
print(f"{'scheme':14s} {'lambda':>7s} {'minp':>10s} {'#Bonf':>6s} {'#BH.05':>7s}")
res = {}
for name, p in schemes.items():
    p = np.clip(p, 1e-300, 1.0)
    res[name] = p
    print(f"{name:14s} {lam(p):7.3f} {p.min():10.2e} {int((p<bonf).sum()):6d} {bh_hits(p):7d}")

# ---- survival of the strongest raw peaks (block leads) ----
df2 = df.copy()
for name, p in res.items():
    df2[name] = p
lead = df2.loc[df2.groupby("block")["raw"].idxmin()].copy()   # each block's strongest raw record
lead["raw_nlp"] = -np.log10(lead["raw"].clip(1e-300))
top = lead.sort_values("raw").head(12)
print("\nTop-12 raw block-lead SV records -- do they survive each scheme? (nlp = -log10 p)")
cols = ["chrom", "pos", "block", "MAF"]
show = top[cols].copy()
for name in schemes:
    show[name] = (-np.log10(top[name].clip(1e-300))).round(2)
show["Bonf_line"] = round(-np.log10(bonf), 2)
print(show.to_string(index=False))
print(f"\n(a record 'survives' scheme X at Bonferroni if its {chr(39)}X{chr(39)} nlp column exceeds Bonf_line={-np.log10(bonf):.2f})")
