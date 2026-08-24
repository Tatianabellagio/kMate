#!/usr/bin/env python
"""Build + execute wza_investigation.ipynb (run in the `basic` env).

Shows, for gen1 SNP x bio1 Kendall:
  - the REAL SNP-per-block distribution (hapFIRE BigLD blocks)
  - the SD-of-Z vs SNPs/block poly fit (deg-2 vs deg-7) over the FULL range,
    making visible where deg-2 turns negative (-> NaN) and deg-7 wiggles/explodes
  - what each SNP cap actually downsamples (no block deleted)
  - CAM5 adjudicated polynomial-free against the local empirical null

  BPY=/global/home/users/tbellg/miniforge3/envs/basic/bin/python
  $BPY analysis/grenenet_selection/r2_gea_nonsnp/wza_investigation/_build_investigation_nb.py
"""
import os
import nbformat as nbf
from nbconvert.preprocessors import ExecutePreprocessor

HERE = "/global/scratch/users/tbellg/kmate/analysis/grenenet_selection/wza_investigation"
OUT = f"{HERE}/wza_investigation.ipynb"

nb = nbf.v4.new_notebook()
C = []

C.append(nbf.v4.new_markdown_cell(r"""# WZA investigation — SNP-per-block distribution & polynomial fit
**gen1 SNP × bio1 (Kendall-τ), hapFIRE BigLD blocks.**

The WZA SNP-number correction fits a polynomial to *(SNPs-per-block → mean & SD of
the window-Z)*, then computes each block's p as `1 − Φ(Z; mean̂, sd̂)`. It assumes a
bounded, densely-sampled SNP-per-block distribution. Our blocks span **1 → 9,158
SNPs** (centromeric + high-LD arm blocks), so the deg-2 polynomial turns **negative**
in the sparse tail → predicted SD < 0 → **NaN p-values** (the Booker-email bug).

Blocks = hapFIRE two-level LD partition: independent **r²=0.1** (window 100) + fine
**BigLD CLQcut r²=0.5** (density). The cap *downsamples* SNPs inside oversized
blocks (100 resamples averaged) — it does **not** delete any block."""))

C.append(nbf.v4.new_code_cell(r"""import sys, os
sys.path.insert(0, "/global/scratch/users/tbellg/kmate/analysis/grenenet_selection/wza_investigation")
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import norm
import wza_core as wc

KEN = ("/global/scratch/users/tbellg/kmate/analysis/grenenet_selection/"
       "phase1_replication/kendall/kendall_snp_gen1_bio1.csv")
df = pd.read_csv(KEN)
df = df[df.MAF >= 0.05].copy()                 # phase-1 maf05
n = df.groupby("block").size()
print(f"{len(n):,} blocks, {len(df):,} SNPs (MAF>=0.05)")
print(f"SNPs/block: min {n.min()}  median {int(n.median())}  q95 {int(n.quantile(.95))}  "
      f"q99 {int(n.quantile(.99))}  max {n.max()}")
base  = wc.raw_wza(df, cap=None)               # raw weighted-Z per block, no cap
cap2k = wc.raw_wza(df, cap=2000)               # cap at 2000 (downsample big blocks)"""))

C.append(nbf.v4.new_markdown_cell(r"""## 1. The real SNP-per-block distribution
Heavily right-skewed: median ~13 SNPs but a long tail to ~9,000. A handful of giant
blocks hold most of the SNP mass."""))

C.append(nbf.v4.new_code_cell(r"""fig, ax = plt.subplots(1, 2, figsize=(13, 4.5))
ax[0].hist(np.log10(n), bins=70, color="steelblue")
for c, lab in [(53, "q75=53"), (638, "q95=638"), (2000, "cap2000")]:
    ax[0].axvline(np.log10(c), ls="--", lw=1.4, label=lab)
ax[0].set(xlabel="log10(SNPs per block)", ylabel="# blocks",
          title=f"SNP-per-block distribution (n={len(n):,})"); ax[0].legend(fontsize=8)

xs = np.unique(np.round(np.logspace(0, np.log10(n.max()), 200)).astype(int))
ax[1].plot(xs, [(n <= x).mean() for x in xs], label="% of blocks ≤ x SNPs")
ax[1].plot(xs, [n[n <= x].sum()/n.sum() for x in xs], label="% of all SNPs in blocks ≤ x")
for c in (53, 638, 2000): ax[1].axvline(c, ls="--", lw=1, color="grey")
ax[1].set(xscale="log", xlabel="SNPs per block (log)", ylabel="cumulative fraction",
          title="blocks vs SNP mass"); ax[1].legend(fontsize=8)
plt.tight_layout(); plt.show()"""))

C.append(nbf.v4.new_markdown_cell(r"""## 2. The polynomial fit — SD of window-Z vs SNPs/block (the failure)
Grey = empirical rolling-bin SD of Z (the data the correction fits). Curves = deg-2
(canonical Booker) and deg-7 (phase-1) fit on the rolling support, **evaluated over
the full SNP range**. Watch deg-2 turn **negative** past the support (→ NaN), and
deg-7 wiggle then explode. **Right panel:** the cap bounds everything to ≤2000, so
no extrapolation → no NaN, for any degree."""))

C.append(nbf.v4.new_code_cell(r"""def rolling_support(rw, roller=50, minE=40):
    s = rw[~rw.Z.isnull()].sort_values("SNPs")
    var = s.Z.rolling(roller, min_periods=minE).var(); m = ~var.isnull()
    return (s.SNPs.rolling(roller, min_periods=minE).mean()[m].to_numpy(),
            np.sqrt(var[m]).to_numpy())

fig, ax = plt.subplots(1, 2, figsize=(14, 5))
for a, (rw, lab, xmax) in zip(ax, [(base, "NO CAP", n.max()), (cap2k, "CAP 2000", 2000)]):
    x, sd = rolling_support(rw)
    a.scatter(x, sd, s=7, color="grey", alpha=.45, label="rolling-bin SD of Z (data)")
    g = np.linspace(0, xmax, 800)
    for deg, col in [(2, "C0"), (7, "C1")]:
        p = np.poly1d(np.polyfit(x, sd, deg))
        rmse = np.sqrt(np.mean((p(x) - sd)**2))
        a.plot(g, p(g), col, lw=2, label=f"deg-{deg} fit (RMSE {rmse:.2f})")
    a.axhline(0, color="red", ls=":", lw=1)
    a.set(xlabel="SNPs per block", ylabel="SD of window-Z",
          title=f"{lab}: predicted SD vs SNPs/block (full range)",
          ylim=(-30, 60)); a.legend(fontsize=9)
plt.tight_layout(); plt.show()
print("deg-2 with no cap dips below 0 in the tail -> norm.cdf returns NaN.")
print("Capping removes the extrapolation; SD stays positive for any degree.")"""))

C.append(nbf.v4.new_markdown_cell(r"""## 3. What each cap actually touches
The cap **downsamples** SNPs inside oversized blocks (cap SNPs/draw × 100 resamples,
averaged). **No block is deleted** — every block still gets a WZA score."""))

C.append(nbf.v4.new_code_cell(r"""rows = []
for c in (2000, int(n.quantile(.95)), int(n.quantile(.75))):
    above = n[n > c]
    rows.append(dict(cap=c, blocks_downsampled=len(above),
                     pct_blocks=round(len(above)/len(n)*100, 1),
                     SNPs_in_those_blocks=int(above.sum()),
                     pct_of_all_SNPs=round(above.sum()/n.sum()*100, 1)))
print(pd.DataFrame(rows).to_string(index=False))
print("\nNote: cap2000 touches ~1% of blocks but they hold ~26% of all SNPs — a few "
      "giant blocks dominate the SNP mass.")"""))

C.append(nbf.v4.new_markdown_cell(r"""## 4. CAM5 adjudicated — polynomial-free
CAM5 block `2_1265` (13 SNPs, raw Z = 17.75 ≈ 18σ, a real coherent signal). Instead
of trusting either polynomial, estimate the null `(mean, SD)` of Z **directly from
blocks with a similar SNP count** and read off CAM5's p. This is the ground truth
the polynomial is supposed to approximate."""))

C.append(nbf.v4.new_code_cell(r"""cam = base[base.block == "2_1265"].iloc[0]
print(f"CAM5 (2_1265): {int(cam.SNPs_raw)} SNPs, raw Z = {cam.Z:.3f}\n")
print("Polynomial-free local empirical null near 13 SNPs:")
for lo, hi in [(11, 15), (10, 16), (8, 20), (5, 25)]:
    nb_ = base[(base.SNPs_raw >= lo) & (base.SNPs_raw <= hi)]
    m, s = nb_.Z.mean(), nb_.Z.std()
    print(f"  {lo:2d}-{hi:2d} SNPs (n={len(nb_):4d}): mean={m:.3f} SD={s:.3f}"
          f"  ->  CAM5 p = {1-norm.cdf(cam.Z, m, s):.2e}")
print("\nWhat each polynomial said at 13 SNPs:")
for tag, deg, cap in [("deg2_nocap", 2, None), ("deg7_nocap", 7, None), ("deg2_cap2000", 2, 2000)]:
    rw = base if cap is None else cap2k
    mp, sp, _ = wc.fit_correction(rw, deg=deg)
    print(f"  {tag:13}: mean̂={mp(13):.2f} sd̂={sp(13):.2f}  ->  p = {1-norm.cdf(cam.Z, mp(13), sp(13)):.2e}")
print("\nTRUTH (local empirical) ~ 1.6e-4.  deg-2 too conservative (1.5e-3); "
      "deg-7 too liberal (5.7e-5) but closer.  Phase-1's published p~3e-8 is deg-7 tail over-fit.")"""))

C.append(nbf.v4.new_markdown_cell(r"""## 5. Takeaways
1. **The cap is the load-bearing fix.** The NaN is an *extrapolation* artifact (deg-2
   predicting past the ~4,500-SNP rolling support out to the 9,158-SNP blocks). Any
   cap removes the extrapolation → no NaN, regardless of degree.
2. **deg-7 genuinely fits the curve better** (lower RMSE; closer at CAM5) — the author
   was right that a quadratic is too rigid here. Its danger is *tail extrapolation*,
   which the cap eliminates. So the author's *(cap=2000 + deg-7)* is internally
   consistent; deg-7 alone (phase-1, no cap) is the fragile part.
3. **deg-2 is the canonical default and more stable** (Bonferroni count swings
   25→6 between degrees; deg-2 less wiggly), but **slightly conservative** at low SNP
   counts. With the cap, deg-2 vs deg-7 is a judgement call — best resolved by the
   **local empirical null** (§4), which both polynomials only approximate.
4. WZA does **not** over-weight low-SNP windows — small windows are *under*-powered;
   deg-7 is what lifts CAM5-sized windows. CAM5 is a real top-~0.1% block (p≈1.6e-4),
   just over-stated by the uncapped deg-7."""))

nb["cells"] = C
nb.metadata["kernelspec"] = {"name": "python3", "display_name": "Python 3"}
print("executing notebook…")
ep = ExecutePreprocessor(timeout=1200, kernel_name="python3", startup_timeout=180)
ep.preprocess(nb, {"metadata": {"path": HERE}})
with open(OUT, "w") as f:
    nbf.write(nb, f)
print("wrote", OUT)
