#!/usr/bin/env python
"""Build+execute fit_inspection.ipynb -- SHOW every candidate SD and mean fit so the choice
can be made by eye as well as by RMSE. The two disagree, and the disagreement is the point.

Background. `cap_poly_decision.ipynb` §4 chose the mean fit on an UNWEIGHTED out-of-sample
RMSE. That metric is ~99.6% determined by the small-block bulk (85% of rolling points are
below 50 records/block; 0.37% above 700), so it could not see that the empirical mean(Z)
falls from ~0 to -5.5 (snp clq0.9) / -7.7 (snp clq0.5) at large blocks. T. Bellagio spotted
that from the §4 figure immediately.

Re-scoring per block-size stratum (mean_fit_stratified.py, job 35979974) did not resolve it,
because even the top 1% band spans Xf 380 -> 3171 and is itself dominated by its lower edge.
The ~10-15 blocks that live past the rolling-support edge are ~0.02% of rolling points --
invisible to ANY aggregate metric. So for those blocks the choice has to be made on the
extrapolated VALUE, not on RMSE:

    kendall / clq0.9 / snp, mean(Z) at the largest block (20,365 records):
      empirical (blocks past support): -10.53
      deg2_clamp -8.54 | deg5_clamp -6.22 | deg7_clamp -5.62 | isotonic_auto -5.47
      const      -0.18   <-- off by ~10
      deg2 (UNCLAMPED) -289.25 | deg7 (UNCLAMPED) -9.2e+08

Two further findings worth seeing plotted:
  * "just use a higher degree" FAILS. Even clamped, deg15/deg20 oscillate (Runge) near the
    sparse upper end of the fitted range: large-block-band RMSE 9.8 and 149 vs ~1.0-1.2 for
    deg2-deg10. Booker's "overfitting is fine" holds only up to ~deg 5-7 on this data.
  * clamping is what makes any polynomial usable at all -- it is the entire difference
    between deg2_clamp (-8.54) and deg2 (-289).

Runs in `basic` env. No chart titles (repo convention).
"""
import os, sys
import nbformat as nbf
from nbformat.v4 import new_notebook, new_markdown_cell, new_code_cell
from nbconvert.preprocessors import ExecutePreprocessor

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
# one notebook per per-record model; file prefix in wza_in_* is the dict key
MODELS = {"kendall": "Kendall's tau", "lfmm": "LFMM (K=16 latent factors)",
          "binomial": "quasi-binomial (K=16 LF, phi-scaled)"}
MODEL = sys.argv[1] if len(sys.argv) > 1 else "lfmm"
assert MODEL in MODELS, f"model must be one of {list(MODELS)}"
MODEL_LABEL = MODELS[MODEL]
OUT = f"{ROOT}/analysis/grenenet_gea/wza_investigation/fit_inspection_{MODEL}.ipynb"

md_title = r"""# Inspecting every candidate WZA fit (SD and mean), on **__MODEL_LABEL__**, per block definition × class

**Why this notebook exists.** The mean fit was previously chosen on an unweighted
out-of-sample RMSE. That metric is ~99.6% determined by the small-block bulk — 85% of rolling
points sit below 50 records/block, only 0.37% above 700 — so it could not see that the
empirical mean(Z) falls from ≈0 to **−5.5** (snp clq0.9) and **−7.7** (snp clq0.5) at large
blocks. `const` sits at −0.2 there. That was visible in the figure long before it was visible
in the table.

Re-scoring per block-size stratum did not settle it either: even the top-1% band spans
Xf 380 → 3171 and is dominated by its own lower edge. **The 8–15 blocks past the
rolling-support edge are ~0.02% of rolling points — invisible to any aggregate metric.** For
those blocks the choice must be made on the extrapolated *value*, which is what §3 shows.

Two things to look for:
1. **Does the fit track the decline** at large blocks, or flatten through it?
2. **Does it stay sane past the support edge** (dashed vertical line)? Everything to the
   right of that line is extrapolation.

`_clamp` = fitted inside support, held at the boundary value outside. That clamp is the
entire difference between `deg2_clamp` (−8.5 at the largest snp block) and plain `deg2`
(**−289**).

## ⚠ Everything in this notebook is **__MODEL_LABEL__**

One notebook is built per per-record model, because the models differ *qualitatively* here —
so a fit chosen on one does not automatically transfer. The §4 figure in
`cap_poly_decision.ipynb` that prompted this re-examination, with mean(Z) plunging to −5.5,
was **kendall** (that notebook hardcodes `model="kendall"`). LFMM, which corrects for
population structure, is far better behaved. Empirical mean(Z) of the blocks past support:

| class (clq0.9 tile) | kendall | **LFMM** |
|---|---|---|
| snp | **−10.53** | **+3.82** |
| smallindel | −5.78 | −0.30 |
| nonsnp | −5.83 | −0.17 |
| sv | −2.24 | +0.52 |

So the dramatic decline is mostly a kendall phenomenon; on LFMM the mean is near zero in 6 of
8 blockdef × class cases, and *rises* for snp. **Note the opposite sign** — that is why no
hardcoded monotone direction (and no `isotonic(increasing=...)`, and no `increasing="auto"`
inferred from a swamped Spearman rho) can serve all three models, while a clamped polynomial
infers nothing and can.

Companion notebooks: `fit_inspection_kendall.ipynb`, `fit_inspection_lfmm.ipynb`,
`fit_inspection_binomial.ipynb` (quasi-binomial)."""

code_setup = r'''
import os
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
import scipy.stats
from sklearn.isotonic import IsotonicRegression
os.chdir("__ROOT__")
MA = "analysis/grenenet_gea/phase1_replication/results/multiaxis"
ROLLER, MINE = 50, 40
DEFS = {"clq0.9 tile": "clq09_tile", "clq0.5 tile": "clq05_tile"}
CLASSES = ["snp", "sv", "smallindel", "nonsnp"]

COL = {"const": "#588157", "deg2_clamp": "#B7C0CC", "deg5_clamp": "#3A5A98",
       "deg7_clamp": "#C9A227", "deg10_clamp": "#8a2a1a", "isotonic": "#111111",
       "isotonic_auto": "#111111", "deg15_clamp": "#d62728", "deg20_clamp": "#9467bd",
       "deg2": "#999999", "deg7": "#ff7f0e", "interp": "#e377c2"}

def block_Z(tag, cls, model, axis="bio1"):
    d = pd.read_csv(f"{MA}/wza_in_{tag}/{model}_{cls}_gen9_{axis}.csv",
                    usecols=["block","MAF","pval"])
    d = d[d.pval.notna() & d.block.notna()]
    pV = np.clip(d.pval.rank(method="first").to_numpy()/len(d), 1e-15, 1-1e-3)
    z = scipy.stats.norm.ppf(1-pV); m = d.MAF.to_numpy(float); pq = m*(1-m)
    t = pd.DataFrame({"b": d.block.to_numpy(), "num": pq*z, "den": pq**2})
    g = t.groupby("b").agg(num=("num","sum"), den=("den","sum"), n=("num","size"))
    g = g[g.n>=2]
    return pd.DataFrame({"SNPs": g.n.astype(float), "Z": g.num/np.sqrt(g.den)})

def curve(g):
    s = g.sort_values("SNPs")
    v = s.Z.rolling(ROLLER, min_periods=MINE).var(); m = ~v.isnull()
    return (np.asarray(s.SNPs.rolling(ROLLER,min_periods=MINE).mean()[m], float),
            np.asarray(np.sqrt(v)[m], float),
            np.asarray(s.Z.rolling(ROLLER,min_periods=MINE).mean()[m], float))

def fit(Xs, y, target, kind):
    if kind == "const":         return np.full_like(target, y.mean())
    if kind == "interp":        return np.interp(target, Xs, y)
    if kind == "isotonic":      return IsotonicRegression(increasing=True, out_of_bounds="clip").fit(Xs,y).predict(target)
    if kind == "isotonic_auto": return IsotonicRegression(increasing="auto", out_of_bounds="clip").fit(Xs,y).predict(target)
    deg = int(kind.replace("deg","").replace("_clamp",""))
    mu, sg = Xs.mean(), Xs.std() or 1.0
    p = np.poly1d(np.polyfit((Xs-mu)/sg, y, deg))
    t = np.clip(target, Xs.min(), Xs.max()) if kind.endswith("_clamp") else target
    return p((t-mu)/sg)

MODEL = "__MODEL__"
CACHE = {}
def get(tag, cls, model):
    k = (tag, cls, model)
    if k not in CACHE:
        g = block_Z(tag, cls, model); CACHE[k] = (g,) + curve(g)
    return CACHE[k]
print("ready")
'''

md_1 = r"""## 1. SD of Z — all viable fits, full range

Black = empirical rolling SD (the thing being predicted). Dashed vertical = the rolling
support edge; everything right of it is extrapolation. A fit that goes ≤ 0 anywhere is
disqualified outright (NaN p, or fabricated `p==0`)."""

code_1 = r'''
KS = ["deg2_clamp","deg5_clamp","deg7_clamp","deg10_clamp","isotonic"]
for model in [MODEL]:
    fig, axes = plt.subplots(2, 4, figsize=(19, 7.5))
    for i, (lab, tag) in enumerate(DEFS.items()):
        for j, cls in enumerate(CLASSES):
            ax = axes[i][j]; g, X, sd, mn = get(tag, cls, model)
            ax.plot(X, sd, color="k", lw=2.6, zorder=6, label="empirical")
            grid = np.linspace(X.min(), g.SNPs.max(), 700)
            for k in KS:
                ax.plot(grid, fit(X, sd, grid, k), color=COL[k], lw=1.5, label=k)
            ax.axvline(X.max(), ls="--", c="k", lw=.9)
            ax.axhline(0, ls=":", c="r", lw=.9)
            ax.set_xscale("log")
            ax.set_ylim(min(-0.3, sd.min()-0.5), sd.max()*1.5)
            ax.set_xlabel("records/block (log)"); ax.set_ylabel("SD of Z")
            ax.annotate(f"{lab}\n{cls}\n{model}", (0.03,0.97), xycoords="axes fraction",
                        fontsize=8.5, va="top", weight="bold")
            if i==0 and j==0: ax.legend(frameon=False, fontsize=7)
            for sp in ["top","right"]: ax.spines[sp].set_visible(False)
    fig.tight_layout(); plt.show()
'''

md_2 = r"""## 2. Mean of Z — all viable fits, full range

This is the figure that started the re-examination: the trend at large blocks is obvious, and
`const` (green) flattens straight through it. **`deg5_clamp` (thick blue) is production** —
it fits the trend, and is clamped at the support edge so it cannot blow up.

`const` is now **retired**: over all 24 blockdef × model × class cells its mean absolute error
at the largest block is **2.81** (worst cell 10.35) versus **1.60** (worst 4.41) for
`deg5_clamp`. Upstream fits the mean with a deg-2 polynomial too — fitting the trend is the
faithful choice, and `const` was a kMate invention."""

code_2 = r'''
KM = ["deg5_clamp","const","deg2_clamp","deg7_clamp","deg10_clamp","isotonic_auto"]
for model in [MODEL]:
    fig, axes = plt.subplots(2, 4, figsize=(19, 7.5))
    for i, (lab, tag) in enumerate(DEFS.items()):
        for j, cls in enumerate(CLASSES):
            ax = axes[i][j]; g, X, sd, mn = get(tag, cls, model)
            ax.plot(X, mn, color="grey", lw=1.2, alpha=.75, zorder=4, label="empirical")
            grid = np.linspace(X.min(), g.SNPs.max(), 700)
            for k in KM:
                pr = "  <-- PRODUCTION" if k == "deg5_clamp" else ""
                ax.plot(grid, fit(X, mn, grid, k), color=COL[k],
                        lw=3.0 if k == "deg5_clamp" else 1.5,
                        zorder=7 if k == "deg5_clamp" else 5, label=k + pr)
            above = g[g.SNPs > X.max()]
            if len(above) >= 8:
                ax.scatter([g.SNPs.max()], [above.Z.mean()], marker="*", s=210, color="red",
                           zorder=8, label="empirical mean, blocks past support")
            ax.axvline(X.max(), ls="--", c="k", lw=.9); ax.axhline(0, ls=":", c="k", lw=.8)
            ax.set_xscale("log")
            ax.set_xlabel("records/block (log)"); ax.set_ylabel("mean of Z")
            ax.annotate(f"{lab}\n{cls}\n{model}", (0.03,0.06), xycoords="axes fraction",
                        fontsize=8.5, va="bottom", weight="bold")
            if i==0 and j==0: ax.legend(frameon=False, fontsize=6.5, loc="upper right")
            for sp in ["top","right"]: ax.spines[sp].set_visible(False)
    fig.tight_layout(); plt.show()
'''

md_3 = r"""## 3. Zoom on the large-block region — where the decision actually lives

Upper 10% of the support range and beyond, linear x. The **red star** is the empirical
mean(Z) of the blocks that sit *past* the support edge — the value a fit ought to be heading
toward. Read off which curve gets closest to it, and note that only ~8–15 blocks per class
are out there, which is why no RMSE could adjudicate this."""

code_3 = r'''
for model in [MODEL]:
    fig, axes = plt.subplots(2, 4, figsize=(19, 7.5))
    for i, (lab, tag) in enumerate(DEFS.items()):
        for j, cls in enumerate(CLASSES):
            ax = axes[i][j]; g, X, sd, mn = get(tag, cls, model)
            lo = np.quantile(X, 0.90); mx = g.SNPs.max()
            sel = X >= lo
            ax.plot(X[sel], mn[sel], color="grey", lw=1.3, alpha=.8, label="empirical")
            grid = np.linspace(lo, mx, 500)
            for k in KM:
                pr = "  <-- PRODUCTION" if k == "deg5_clamp" else ""
                ax.plot(grid, fit(X, mn, grid, k), color=COL[k],
                        lw=3.2 if k == "deg5_clamp" else 1.6,
                        zorder=7 if k == "deg5_clamp" else 5, label=k + pr)
            above = g[g.SNPs > X.max()]
            if len(above) >= 8:
                ax.scatter([mx], [above.Z.mean()], marker="*", s=260, color="red", zorder=8,
                           label=f"empirical past support (n={len(above)})")
            ax.axvline(X.max(), ls="--", c="k", lw=1.0)
            ax.set_xlabel("records/block"); ax.set_ylabel("mean of Z")
            ax.annotate(f"{lab}\n{cls}\n{model}", (0.03,0.06), xycoords="axes fraction",
                        fontsize=8.5, va="bottom", weight="bold")
            if i==0 and j==0: ax.legend(frameon=False, fontsize=6.5)
            for sp in ["top","right"]: ax.spines[sp].set_visible(False)
    fig.tight_layout(); plt.show()
'''

md_4 = r"""## 4. Why "just use a higher degree" does not work

Booker's *"this will kind of overfit the data, but I think this is ok as the goal is simply
to interpolate"* is right in spirit, but only up to a point on this data. Even **clamped**,
deg-15 and deg-20 oscillate wildly (Runge phenomenon) near the sparse upper end of the fitted
range, where rolling points thin out. Large-block-band out-of-sample RMSE:

| fit | SD | mean |
|---|---|---|
| deg5_clamp | **1.174** | 1.063 |
| deg7_clamp | 1.257 | 1.096 |
| deg10_clamp | 1.301 | 1.185 |
| deg15_clamp | **11.44** | **9.78** |
| deg20_clamp | **216.5** | **149.1** |

Plotted below on the worst case. Also shown: the same fits *unclamped*, which is the −289
failure."""

code_4 = r'''
fig, axes = plt.subplots(1, 3, figsize=(17, 4.6))
g, X, sd, mn = get("clq09_tile", "snp", "kendall")
grid = np.linspace(X.min(), g.SNPs.max(), 900)
for ax, ks, ttl, y in [
        (axes[0], ["deg5_clamp","deg10_clamp","deg15_clamp","deg20_clamp"], "clamped, high degree", mn),
        (axes[1], ["deg2","deg7"], "UNCLAMPED (the -289 failure)", mn),
        (axes[2], ["deg5_clamp","deg15_clamp","deg20_clamp","isotonic"], "SD, clamped high degree", sd)]:
    ax.plot(X, y, color="grey", lw=1.4, alpha=.8, label="empirical")
    for k in ks:
        ax.plot(grid, fit(X, y, grid, k), color=COL[k], lw=1.6, label=k)
    ax.axvline(X.max(), ls="--", c="k", lw=.9)
    ax.set_xscale("log"); ax.set_xlabel("records/block (log)")
    ax.set_ylabel("mean of Z" if y is mn else "SD of Z")
    lim = np.percentile(np.abs(y), 99) * 6
    ax.set_ylim(-lim, lim)
    ax.annotate(f"clq0.9 tile / snp / kendall\n{ttl}\n(y clipped to +/-{lim:.0f})",
                (0.03,0.04), xycoords="axes fraction", fontsize=8, va="bottom")
    ax.legend(frameon=False, fontsize=7.5)
    for sp in ["top","right"]: ax.spines[sp].set_visible(False)
fig.tight_layout(); plt.show()
'''

md_5 = r"""## 5. The numbers behind the zoom

Predicted mean(Z) at each class's largest block, against what those blocks empirically
average. This is the table the decision rests on — RMSE cannot see this region."""

code_5 = r'''
rows = []
for model in [MODEL]:
    for lab, tag in DEFS.items():
        for cls in CLASSES:
            g, X, sd, mn = get(tag, cls, model)
            mx = np.array([g.SNPs.max()], float)
            above = g[g.SNPs > X.max()]
            r = dict(model=model, blockdef=lab, cls=cls, support_to=round(float(X.max())),
                     maxN=int(mx[0]), n_past=len(above),
                     empirical=above.Z.mean() if len(above) >= 8 else np.nan)
            for k in ["const","deg2_clamp","deg5_clamp","deg7_clamp","isotonic_auto","deg2"]:
                r[k] = fit(X, mn, mx, k)[0]
            rows.append(r)
R = pd.DataFrame(rows)
print(R.round(2).to_string(index=False))
print("\nabsolute error vs the empirical past-support mean (lower better):")
E = R.dropna(subset=["empirical"]).copy()
for k in ["const","deg2_clamp","deg5_clamp","deg7_clamp","isotonic_auto"]:
    E["e_"+k] = (E[k] - E.empirical).abs()
print(E.groupby("model")[["e_const","e_deg2_clamp","e_deg5_clamp","e_deg7_clamp",
                          "e_isotonic_auto"]].mean().round(2).to_string())
print("\noverall mean |error|:")
print(E[[c for c in E.columns if c.startswith("e_")]].mean().sort_values().round(2).to_string())
'''

code_setup = code_setup.replace("__ROOT__", ROOT).replace("__MODEL__", MODEL)
md_title = md_title.replace("__MODEL_LABEL__", MODEL_LABEL)
nb = new_notebook(cells=[
    new_markdown_cell(md_title), new_code_cell(code_setup),
    new_markdown_cell(md_1), new_code_cell(code_1),
    new_markdown_cell(md_2), new_code_cell(code_2),
    new_markdown_cell(md_3), new_code_cell(code_3),
    new_markdown_cell(md_4), new_code_cell(code_4),
    new_markdown_cell(md_5), new_code_cell(code_5),
])
ep = ExecutePreprocessor(timeout=7200, kernel_name="basic", startup_timeout=180)
ep.preprocess(nb, {"metadata": {"path": os.path.dirname(OUT)}})
with open(OUT, "w") as f:
    nbf.write(nb, f)
print(f"[built+executed] {OUT}")
