#!/usr/bin/env python
"""AUDIT: is `isotonic_auto` a SAFE mean fit across all 60 real scans, or does its
direction inference flip on noise?

Why this matters. `IsotonicRegression(increasing="auto")` picks the monotone direction from
a Spearman test on (rolling mean SNP number, rolling mean Z). Both the decision notebook and
the permutation-null jobs emitted:

    UserWarning: Confidence interval of the Spearman correlation coefficient spans zero.
    Determination of ``increasing`` may be suspect.

If that happens on REAL scans, the mean fit's direction is set by noise, and two adjacent
climate axes could get opposite-trending mean corrections -- an arbitrary, invisible
inconsistency across the 60 (model x axis) scans that the candidate-gene unions are built on.

The honest alternative is `const` (global mean of the rolling means): it was the best
out-of-sample mean fit in 9/12 measured cases anyway, has NO direction to infer, and is
bounded by construction. `isotonic_auto` is only preferable if the trend it tracks is REAL
and CONSISTENT -- which is exactly what this script tests.

Verdict rule:
  isotonic_auto is safe iff, across every model x axis x class, the Spearman test on the
  real (unpermuted) rolling curve is (a) unambiguous (CI does not span zero) and (b) the
  SIGN is consistent within a model. Otherwise prefer `const`.

Read-only over results/multiaxis/wza_in_clq09_tile/ AND wza_in_clq05_tile/. Writes mean_fit_direction_audit.csv.
"""
from __future__ import annotations
import os, sys, warnings
import numpy as np
import pandas as pd
from scipy.stats import norm, spearmanr

warnings.filterwarnings("ignore")
HERE = os.path.dirname(os.path.abspath(__file__))
GEA = os.path.dirname(HERE)
MA = f"{GEA}/phase1_replication/results/multiaxis"

AXES = ["bio%d" % i for i in range(1, 20)] + ["pc1"]
MODELS = ["kendall", "lfmm", "binomial"]
CLASSES = ["snp", "sv", "smallindel", "nonsnp"]
# BOTH block definitions: clq0.5 blocks are much larger (p99 941 vs 380 records for snp), so
# the mean-vs-block-size trend could plausibly be stronger/more consistent there. Auditing
# only clq0.9 would decide the mean fit on half the evidence.
TAGS = {"clq0.9 tile": "clq09_tile", "clq0.5 tile": "clq05_tile"}
ROLLER, MINE = 50, 40


def rolling_curve(f):
    """Production block-Z (rank-transformed p, MAF-weighted Stouffer) -> rolling mean/SD."""
    d = pd.read_csv(f, usecols=["block", "MAF", "pval"])
    d = d[d.pval.notna() & d.block.notna()]
    pV = np.clip(d.pval.rank(method="first").to_numpy() / len(d), 1e-15, 1 - 1e-3)
    z = norm.ppf(1 - pV)
    m = d.MAF.to_numpy(float)
    pq = m * (1 - m)
    t = pd.DataFrame({"b": d.block.to_numpy(), "num": pq * z, "den": pq ** 2})
    g = t.groupby("b").agg(num=("num", "sum"), den=("den", "sum"), n=("num", "size"))
    g = g[g.n >= 2]
    res = pd.DataFrame({"SNPs": g.n.astype(float), "Z": g.num / np.sqrt(g.den)})
    s = res.sort_values("SNPs")
    v = s.Z.rolling(ROLLER, min_periods=MINE).var()
    msk = ~v.isnull()
    X = np.asarray(s.SNPs.rolling(ROLLER, min_periods=MINE).mean()[msk], float)
    mn = np.asarray(s.Z.rolling(ROLLER, min_periods=MINE).mean()[msk], float)
    return X, mn


def main():
    rows = []
    for blockdef, TAG in TAGS.items():
      for model in MODELS:
        for cls in CLASSES:
            for axis in AXES:
                f = f"{MA}/wza_in_{TAG}/{model}_{cls}_gen9_{axis}.csv"
                if not os.path.exists(f):
                    rows.append(dict(blockdef=blockdef, model=model, cls=cls, axis=axis,
                                     missing=True))
                    continue
                X, mn = rolling_curve(f)
                rho, p = spearmanr(X, mn)
                # sklearn's own ambiguity criterion (isotonic.check_increasing): Fisher-z CI
                # on the Spearman rho; "suspect" when it spans zero.
                n = len(X)
                F = 0.5 * np.log((1 + rho) / (1 - rho))
                se = 1.0 / np.sqrt(n - 3)
                lo, hi = np.tanh(F - 1.97 * se), np.tanh(F + 1.97 * se)
                rows.append(dict(blockdef=blockdef, model=model, cls=cls, axis=axis,
                                 missing=False,
                                 npts=n, rho=rho, p=p, ci_lo=lo, ci_hi=hi,
                                 ambiguous=bool(lo <= 0 <= hi),
                                 direction="increasing" if rho > 0 else "decreasing"))
                print(f"  {blockdef:12s} {model:9s} {cls:11s} {axis:5s} rho={rho:+.3f} "
                      f"CI[{lo:+.3f},{hi:+.3f}] {'AMBIGUOUS' if lo <= 0 <= hi else ''}", flush=True)
    D = pd.DataFrame(rows)
    D.to_csv(f"{HERE}/mean_fit_direction_audit.csv", index=False)
    ok = D[~D.missing]
    print("\n" + "=" * 78)
    print(f"scans audited: {len(ok)} (missing {int(D.missing.sum())})")
    print(f"AMBIGUOUS (Fisher-z CI on rho spans 0 -> direction set by noise): "
          f"{int(ok.ambiguous.sum())}/{len(ok)}")
    print("\nsign of the trend, per blockdef x model x class (want ONE direction per cell):")
    print(pd.crosstab([ok.blockdef, ok.model, ok.cls], ok.direction).to_string())
    print("\nper blockdef x model x class: min |rho| and #ambiguous")
    g = ok.assign(absrho=ok.rho.abs()).groupby(["blockdef", "model", "cls"]).agg(
        min_absrho=("absrho", "min"), n_ambiguous=("ambiguous", "sum"),
        n_increasing=("direction", lambda s: (s == "increasing").sum()),
        n_decreasing=("direction", lambda s: (s == "decreasing").sum()))
    print(g.to_string())
    flips = g[(g.n_increasing > 0) & (g.n_decreasing > 0)]
    print("\n" + "=" * 78)
    if len(flips) == 0 and ok.ambiguous.sum() == 0:
        print("VERDICT: isotonic_auto is SAFE -- direction unambiguous and consistent everywhere.")
    else:
        print("VERDICT: isotonic_auto is NOT safe. Prefer mean_fit=const.")
        if int(ok.ambiguous.sum()):
            print(f"  - {int(ok.ambiguous.sum())} scans have an ambiguous direction.")
        if len(flips):
            print(f"  - {len(flips)} blockdef x model x class cells contain BOTH directions "
                  f"across the 20 axes (out of {len(g)}):")
            print(flips.to_string())
    print(f"\nwrote {HERE}/mean_fit_direction_audit.csv")


if __name__ == "__main__":
    main()
