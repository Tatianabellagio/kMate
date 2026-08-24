#!/usr/bin/env python
"""FINAL per-cell WZA fit decision: best SD fit and best mean fit for every
block partition x model x variant class (2 x 3 x 4 = 24 cells).

WHY A REPRODUCIBILITY TEST AND NOT PLAIN OUT-OF-SAMPLE RMSE
-----------------------------------------------------------
Half-split OOS RMSE kept selecting `const` for the mean even where the full-data rolling
curve visibly trends at large blocks. That is a bias-variance artifact, not evidence:
each half contains only ~5 of the ~10 blocks past the support edge, so the VALIDATION
TARGET is extremely noisy out there. A zero-variance predictor (`const`) beats a
trend-follower on RMSE against a noisy target even when the trend is real. In production
the fit uses ALL the blocks, where the trend is far better estimated.

So the question is not "which fit predicts a noisy half" but "DOES THE TREND REPLICATE".
Test: split the blocks in half 20x (40 halves); in each half fit a slope to mean(Z) vs
log10(block size) over the top 10% of the support range; record the slope's sign. If the
sign agrees across nearly all halves, the trend is real and must be fitted. If it is a
coin flip, it is noise and `const` is correct.

DECISION RULE (pre-specified)
  sign agreement >= 0.90  -> trend is REAL     -> fit it: best of {deg2_clamp, deg5_clamp,
                                                  isotonic_auto} by top-band OOS RMSE
  sign agreement <  0.90  -> not established   -> `const` (do not fit noise)

  SD: always `isotonic`. The SD trend is real everywhere by construction (Var of a
  weighted sum of correlated z grows with block size), isotonic tracks it, is monotone,
  and CANNOT go negative -- the failure mode that started this whole investigation. A
  single isotonic rule costs only 1.7% vs per-cell winners, whose margins are 1-3%.

ALL POLYNOMIALS ARE CLAMPED at the support boundary (fit inside, held flat outside).
Unclamped is what produces mean(Z) = -289 at the largest snp block. Degrees above ~10 are
excluded entirely: they Runge-oscillate near the sparse upper end (top-band RMSE 11-216).

Writes wza_fit_decision.csv and prints a paste-ready dict for wza_fit_config.py.
"""
from __future__ import annotations
import os, sys, warnings
import numpy as np
import pandas as pd
from scipy.stats import norm
from sklearn.isotonic import IsotonicRegression

warnings.filterwarnings("ignore")
HERE = os.path.dirname(os.path.abspath(__file__))
GEA = os.path.dirname(HERE)
MA = f"{GEA}/r2_gea_nonsnp/phase1_replication/results/multiaxis"

DEFS = {"clq0.9": "clq09_tile", "clq0.5": "clq05_tile"}
MODELS = ["kendall", "lfmm", "binomial"]
CLASSES = ["snp", "sv", "smallindel", "nonsnp"]
ROLLER, MINE = 50, 40
NSPLIT = 20
AGREE_THRESH = 0.90
TREND_CANDS = ["deg2_clamp", "deg5_clamp", "isotonic_auto"]


def block_Z(tag, cls, model, axis="bio1"):
    d = pd.read_csv(f"{MA}/wza_in_{tag}/{model}_{cls}_gen9_{axis}.csv",
                    usecols=["block", "MAF", "pval"])
    d = d[d.pval.notna() & d.block.notna()]
    pV = np.clip(d.pval.rank(method="first").to_numpy() / len(d), 1e-15, 1 - 1e-3)
    z = norm.ppf(1 - pV)
    m = d.MAF.to_numpy(float)
    pq = m * (1 - m)
    t = pd.DataFrame({"b": d.block.to_numpy(), "num": pq * z, "den": pq ** 2})
    g = t.groupby("b").agg(num=("num", "sum"), den=("den", "sum"), n=("num", "size"))
    g = g[g.n >= 2]
    return pd.DataFrame({"SNPs": g.n.astype(float), "Z": g.num / np.sqrt(g.den)})


def curve(g):
    s = g.sort_values("SNPs")
    v = s.Z.rolling(ROLLER, min_periods=MINE).var()
    m = ~v.isnull()
    return (np.asarray(s.SNPs.rolling(ROLLER, min_periods=MINE).mean()[m], float),
            np.asarray(np.sqrt(v)[m], float),
            np.asarray(s.Z.rolling(ROLLER, min_periods=MINE).mean()[m], float))


def fit(Xs, y, target, kind):
    if kind == "const":
        return np.full_like(target, y.mean())
    if kind == "isotonic":
        return IsotonicRegression(increasing=True, out_of_bounds="clip").fit(Xs, y).predict(target)
    if kind == "isotonic_auto":
        return IsotonicRegression(increasing="auto", out_of_bounds="clip").fit(Xs, y).predict(target)
    deg = int(kind.replace("deg", "").replace("_clamp", ""))
    mu, sg = Xs.mean(), (Xs.std() or 1.0)
    p = np.poly1d(np.polyfit((Xs - mu) / sg, y, deg))
    return p((np.clip(target, Xs.min(), Xs.max()) - mu) / sg)


def main():
    rng = np.random.default_rng(0)
    rows = []
    for bd, tag in DEFS.items():
        for model in MODELS:
            for cls in CLASSES:
                g = block_Z(tag, cls, model)
                X, sd, mn = curve(g)
                lo = np.quantile(X, 0.90)
                top = X >= lo
                slope_full = np.polyfit(np.log10(X[top]), mn[top], 1)[0]

                # --- reproducibility of the large-block trend --------------------------
                slopes, rm = [], {k: [] for k in TREND_CANDS + ["const"]}
                for rep in range(NSPLIT):
                    ix = rng.permutation(len(g))
                    h = len(g) // 2
                    Xa, sda, mna = curve(g.iloc[ix[:h]])
                    Xb, sdb, mnb = curve(g.iloc[ix[h:]])
                    for Xh, mh in ((Xa, mna), (Xb, mnb)):
                        s2 = Xh >= lo
                        if s2.sum() >= 8:
                            slopes.append(np.polyfit(np.log10(Xh[s2]), mh[s2], 1)[0])
                    if len(Xa) < 10 or len(Xb) < 10:
                        continue
                    selb = Xb >= lo
                    if selb.sum() >= 3:
                        for k in rm:
                            rm[k].append(np.sqrt(np.nanmean(
                                (fit(Xa, mna, Xb, k)[selb] - mnb[selb]) ** 2)))
                slopes = np.array(slopes)
                pos = float((slopes > 0).mean())
                agree = max(pos, 1 - pos)
                real = agree >= AGREE_THRESH

                scores = {k: float(np.mean(v)) for k, v in rm.items() if v}
                if real:
                    best = min(TREND_CANDS, key=lambda k: scores.get(k, np.inf))
                else:
                    best = "const"

                above = g[g.SNPs > X.max()]
                emp = above.Z.mean() if len(above) >= 8 else np.nan
                mx = np.array([g.SNPs.max()], float)
                rows.append(dict(
                    blockdef=bd, model=model, cls=cls,
                    slope_full=round(slope_full, 2), sign_agree=round(agree, 2),
                    trend_real=real, mean_fit=best, sd_fit="isotonic",
                    n_past_support=len(above),
                    empirical_tail=None if pd.isna(emp) else round(float(emp), 2),
                    pred_chosen=round(float(fit(X, mn, mx, best)[0]), 2),
                    pred_const=round(float(fit(X, mn, mx, "const")[0]), 2),
                    **{f"rmse_{k}": round(scores.get(k, np.nan), 3) for k in TREND_CANDS + ["const"]}))
                print(f"  {bd:7s} {model:9s} {cls:11s} slope {slope_full:+6.2f} "
                      f"agree {agree:.0%} -> {best}", flush=True)

    D = pd.DataFrame(rows)
    D.to_csv(f"{HERE}/wza_fit_decision.csv", index=False)

    print("\n" + "=" * 118)
    print("FINAL PER-CELL DECISION (SD is isotonic everywhere; mean varies)\n")
    show = ["blockdef", "model", "cls", "slope_full", "sign_agree", "trend_real",
            "mean_fit", "n_past_support", "empirical_tail", "pred_chosen", "pred_const"]
    print(D[show].to_string(index=False))
    print("\nmean-fit family counts:", D.mean_fit.value_counts().to_dict())
    print("cells where the trend replicates (>=90%):", int(D.trend_real.sum()), "/", len(D))
    print("\nby blockdef x model:")
    print(pd.crosstab([D.blockdef, D.model], D.mean_fit).to_string())

    print("\n" + "=" * 118)
    print("PASTE-READY (wza_fit_config.py):\n")
    print("MEAN_FIT = {")
    for _, r in D.iterrows():
        print(f'    ("{r.blockdef}", "{r.model}", "{r.cls}"): "{r.mean_fit}",')
    print("}")
    print('SD_FIT = "isotonic"   # all cells')


if __name__ == "__main__":
    main()
