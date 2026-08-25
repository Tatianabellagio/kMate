#!/usr/bin/env python
"""RE-DO the WZA SD *and* mean fit choice with (a) higher-order polynomials and (b) a metric
that can actually SEE the large-block region.

WHY THIS EXISTS -- correction to cap_poly_decision.ipynb §4 and the `const` mean decision.

1. THE METRIC WAS BLIND. The earlier selection used an UNWEIGHTED out-of-sample RMSE over
   rolling points. For snp/clq0.9 tile, 85% of rolling points sit below 50 records/block and
   only 0.37% above 700 -- so that RMSE was ~99.6% determined by the small-block bulk, where
   mean(Z) really is ~0 and pure noise. It could not see the large-block region.

2. AND THE MEAN IS CLEARLY NOT CONSTANT THERE. The empirical rolling mean(Z) falls from ~0
   to -5.5 (snp clq0.9) and -7.7 (snp clq0.5). `const` sits at -0.2, i.e. badly wrong exactly
   where blocks are big enough to matter. Same reason the direction audit measured
   |rho|~0.01-0.06: the swamping noisy bulk makes the OVERALL rank correlation tiny even
   though the large-N trend is unmistakable in the plots.

3. BOOKER SAYS OVERFIT ON PURPOSE. "This will kind of overfit the data, but I think this is
   ok as the goal is simply to interpolate across SNP number bins. I found a degree of 7
   worked ok for that." So high-degree polynomials are the author-endorsed approach; we just
   have to pick the degree on OUR block-size distribution instead of inheriting 7.

4. THE ONE HARD CONSTRAINT HE DID NOT FACE. Our target values run far past the fitted range
   (largest snp block 20,365 vs rolling support ending at 3,171) because his windows had
   tight SNP counts. An unbounded polynomial evaluated out there explodes -- that IS the
   deg-2 mean predicting -289 and the deg-7 SD RMSE of 17.7. So a flexible fit must be
   CLAMPED at the support boundary. Clamped, a high-degree poly is exactly the
   overfit-to-interpolate Booker recommends, with the extrapolation risk removed.

NUMERICAL NOTE: np.polyfit on x in [2, 4000] is badly conditioned past ~deg 10, so x is
standardized before fitting and the same transform applied to the target. Without this,
high degrees return garbage rather than a real overfit.

Evaluates BOTH targets (SD of Z, mean of Z) per block-size stratum, plus extrapolation
sanity, plus (for SD) how often a fit goes <= 0, which is disqualifying.
Models: kendall AND lfmm -- their large-N mean trends have OPPOSITE sign, so no hardcoded
monotone direction serves both, while an unconstrained clamped poly can.

Writes mean_fit_stratified.csv / mean_fit_extrapolation.csv.
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

DEFS = {"clq0.9 tile": "clq09_tile", "clq0.5 tile": "clq05_tile"}
CLASSES = ["snp", "sv", "smallindel", "nonsnp"]
MODELS = ["kendall", "lfmm", "binomial"]
ROLLER, MINE, NREP = 50, 40, 10
QBANDS = [(0.0, 0.5), (0.5, 0.9), (0.9, 0.99), (0.99, 1.0)]

DEGREES = [2, 5, 7, 10, 15, 20]
# clamped polys (safe) at every degree, plus the unclamped deg2/deg7 as the historical
# reference, plus the non-parametric options
SD_CANDS = [f"deg{d}_clamp" for d in DEGREES] + ["deg2", "deg7", "isotonic"]
MEAN_CANDS = [f"deg{d}_clamp" for d in DEGREES] + ["deg2", "deg7", "isotonic_auto",
                                                   "const", "interp"]


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
    """WZA's own rolling estimate (wza_script.py:64-68) -> (X, sd, mean)."""
    s = g.sort_values("SNPs")
    v = s.Z.rolling(ROLLER, min_periods=MINE).var()
    m = ~v.isnull()
    return (np.asarray(s.SNPs.rolling(ROLLER, min_periods=MINE).mean()[m], float),
            np.asarray(np.sqrt(v)[m], float),
            np.asarray(s.Z.rolling(ROLLER, min_periods=MINE).mean()[m], float))


def fit(Xs, y, target, kind):
    if kind == "const":
        return np.full_like(target, y.mean())
    if kind == "interp":
        return np.interp(target, Xs, y)
    if kind == "isotonic":
        return IsotonicRegression(increasing=True, out_of_bounds="clip").fit(Xs, y).predict(target)
    if kind == "isotonic_auto":
        return IsotonicRegression(increasing="auto", out_of_bounds="clip").fit(Xs, y).predict(target)
    deg = int(kind.replace("deg", "").replace("_clamp", ""))
    # standardize x so high degrees are numerically well conditioned
    mu, sg = Xs.mean(), Xs.std()
    sg = sg if sg > 0 else 1.0
    p = np.poly1d(np.polyfit((Xs - mu) / sg, y, deg))
    t = np.clip(target, Xs.min(), Xs.max()) if kind.endswith("_clamp") else target
    return p((t - mu) / sg)


def main():
    rows, xrows = [], []
    rng = np.random.default_rng(0)
    for lab, tag in DEFS.items():
        for model in MODELS:
            for cls in CLASSES:
                g = block_Z(tag, cls, model)
                Xf, sdf, mnf = curve(g)
                edges = [np.quantile(Xf, q) for q in [0.0, 0.5, 0.9, 0.99, 1.0]]
                targets = {"sd": (sdf, SD_CANDS), "mean": (mnf, MEAN_CANDS)}
                acc = {(w, k, bi): [] for w in targets for k in targets[w][1]
                       for bi in range(len(QBANDS))}
                for rep in range(NREP):
                    ix = rng.permutation(len(g))
                    h = len(g) // 2
                    Xa, sda, mna = curve(g.iloc[ix[:h]])
                    Xb, sdb, mnb = curve(g.iloc[ix[h:]])
                    if len(Xa) < 10 or len(Xb) < 10:
                        continue
                    for w, (_, cands) in targets.items():
                        ya = sda if w == "sd" else mna
                        yb = sdb if w == "sd" else mnb
                        for k in cands:
                            pred = fit(Xa, ya, Xb, k)
                            for bi in range(len(QBANDS)):
                                sel = (Xb >= edges[bi]) & (Xb <= edges[bi + 1])
                                if sel.sum() >= 3:
                                    acc[(w, k, bi)].append(
                                        np.sqrt(np.nanmean((pred[sel] - yb[sel]) ** 2)))
                for (w, k, bi), v in acc.items():
                    if v:
                        rows.append(dict(blockdef=lab, model=model, cls=cls, target=w, fit=k,
                                         band=f"{QBANDS[bi][0]:.2f}-{QBANDS[bi][1]:.2f}",
                                         band_lo=round(edges[bi]), band_hi=round(edges[bi + 1]),
                                         rmse=np.mean(v)))
                # full-fit diagnostics: extrapolation at the largest block + negative-SD count
                allN = g.SNPs.to_numpy(float)
                mx = np.array([g.SNPs.max()], float)
                above = g[g.SNPs > Xf.max()]
                base = dict(blockdef=lab, model=model, cls=cls, maxN=int(mx[0]),
                            support_to=round(float(Xf.max())), n_past_support=len(above))
                r = dict(base, target="mean",
                         empirical_tail=above.Z.mean() if len(above) >= 8 else np.nan)
                for k in MEAN_CANDS:
                    r[k] = fit(Xf, mnf, mx, k)[0]
                xrows.append(r)
                r = dict(base, target="sd",
                         empirical_tail=above.Z.std() if len(above) >= 15 else np.nan)
                for k in SD_CANDS:
                    r[k] = fit(Xf, sdf, mx, k)[0]
                    r["neg_" + k] = int((fit(Xf, sdf, allN, k) <= 0).sum())
                xrows.append(r)
                print(f"  done {lab} {model} {cls}", flush=True)

    D = pd.DataFrame(rows)
    X = pd.DataFrame(xrows)
    D.to_csv(f"{HERE}/mean_fit_stratified.csv", index=False)
    X.to_csv(f"{HERE}/mean_fit_extrapolation.csv", index=False)

    for w, cands in [("sd", SD_CANDS), ("mean", MEAN_CANDS)]:
        print("\n" + "=" * 118)
        print(f"OUT-OF-SAMPLE RMSE, target = {w.upper()}, BY BLOCK-SIZE STRATUM")
        print("the 0.99-1.00 band is the large-block region the unweighted RMSE was blind to\n")
        for lab in DEFS:
            for model in MODELS:
                sub = D[(D.blockdef == lab) & (D.model == model) & (D.target == w)]
                if not len(sub):
                    continue
                piv = sub.pivot_table(index=["cls", "band"], columns="fit", values="rmse")
                piv = piv[[c for c in cands if c in piv.columns]]
                print(f"--- {lab} / {model} / {w} ---")
                print(piv.round(3).to_string())
                print()
        top = D[(D.target == w) & (D.band == f"{QBANDS[-1][0]:.2f}-{QBANDS[-1][1]:.2f}")]
        v = top.groupby("fit").rmse.mean().sort_values()
        print(f"  >>> {w.upper()} large-block band (0.99-1.00), mean RMSE over all 16 cases:")
        print("      " + v.round(3).to_string().replace("\n", "\n      "))
        allb = D[D.target == w].groupby("fit").rmse.mean().sort_values()
        print(f"  >>> {w.upper()} averaged over ALL bands equally:")
        print("      " + allb.round(3).to_string().replace("\n", "\n      "))

    print("\n" + "=" * 118)
    print("SD FITS THAT EVER GO <= 0 (disqualifying -- yields NaN or fabricated p==0):\n")
    sdx = X[X.target == "sd"]
    negcols = [c for c in sdx.columns if c.startswith("neg_")]
    print(sdx.set_index(["blockdef", "model", "cls"])[negcols].to_string())

    print("\n" + "=" * 118)
    print("EXTRAPOLATION AT THE LARGEST BLOCK vs empirical tail\n")
    for w, cands in [("mean", MEAN_CANDS), ("sd", SD_CANDS)]:
        sub = X[X.target == w]
        cols = ["blockdef", "model", "cls", "support_to", "maxN", "n_past_support",
                "empirical_tail"] + cands
        print(f"--- {w} ---")
        print(sub[cols].round(2).to_string(index=False))
        print()


if __name__ == "__main__":
    main()
