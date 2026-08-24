#!/usr/bin/env python
"""WHY is block-level significance ~170x anti-conservative, and what threshold would be honest?

`perm_null_regime.py` established (jobs 35975820 pool / 35976724 site) that under a
signal-free site-permutation null the settled WZA regime still yields ~8.5 Bonferroni
blocks per replicate when the nominal expectation is 0.05 -- and that this is
  * near-identical across all four correction regimes  -> not the mean/SD fit, and
  * near-identical at n=31 sites vs n=352 pools         -> not pseudoreplication.

Remaining hypothesis: the correction standardizes only the FIRST TWO MOMENTS and then
reads p off a NORMAL reference (`1 - norm.cdf(Z, mean_pred, sd_pred)`, wza_script.py:125).
Within-block LD makes many variants carry the same signal, so the null block-Z is
HEAVY-TAILED. Standardizing mean and SD cannot change tail shape, and Bonferroni at ~57k
blocks probes p~9e-7 (~4.8 SD) -- exactly where a fat tail does maximum damage.

This script measures that directly, per variant class:
  * excess kurtosis of the null standardized Z (Normal = 0)
  * empirical vs Normal quantiles out to 0.9999
  * observed hits/replicate at the NOMINAL Bonferroni critical value
  * an HONEST critical value, via a generalized-Pareto fit to the permutation tail
    (the raw 1-0.05/nblocks quantile is NOT estimable from ~10 perms x nblocks samples,
    so the GPD extrapolation is the estimator -- and it doubles as a demo of the
    proposed empirical-null fix)

Run via sbatch (needs ~196G for the snp class). Writes null_tail_shape.csv.
"""
from __future__ import annotations
import os, sys, warnings
import numpy as np
import pandas as pd
from scipy.stats import norm, t as tdist, kurtosis, genpareto
from sklearn.isotonic import IsotonicRegression

warnings.filterwarnings("ignore")
HERE = os.path.dirname(os.path.abspath(__file__))
GEA = os.path.dirname(HERE)
sys.path.insert(0, GEA)
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), "blocks"))
import blocks_tiling as bt   # noqa: E402

CM = f"{GEA}/r2_gea_nonsnp/phase1_replication/results/class_matrices"
CLASSES = ["snp", "sv", "smallindel", "nonsnp"]
R2, NPERM, ROLLER, MINE, MAF_MIN = 0.9, 10, 50, 40, 0.05
TAIL_Q = 0.999          # GPD threshold: fit the top 0.1% of the null


def unitrank(v):
    r = pd.Series(v).rank().to_numpy().astype(np.float64)
    r -= r.mean()
    n = np.linalg.norm(r)
    return r / n if n else r


def main():
    pools = pd.read_csv(f"{CM}/gen9.pools.csv")
    site = pools["site"].to_numpy()
    us = np.unique(site)
    six = np.array([{s: i for i, s in enumerate(us)}[s] for s in site])
    sb = pools.groupby("site")["bio1"].first().reindex(us).to_numpy()
    npool = len(pools)
    rows = []

    for cls in CLASSES:
        rec = pd.read_csv(f"{CM}/{cls}_gen9.records.csv")
        rec["block"] = bt.merge_small_blocks(
            bt.assign_tiling(rec["chrom"].to_numpy(str), rec["pos"].to_numpy(np.int64), r2=R2))
        k = (rec["block"].to_numpy() != "") & (rec["maf"].to_numpy() > MAF_MIN)
        sub = rec[k]
        blocks = sub["block"].to_numpy().astype(str)
        pq = (sub["maf"] * (1 - sub["maf"])).to_numpy()
        af = np.load(f"{CM}/{cls}_gen9_af.npy", mmap_mode="r")
        A = np.asarray(af[:, sub.index.to_numpy()], dtype=np.float32)
        del af
        rk = A.argsort(0).argsort(0).astype(np.float32) + 1.0
        del A
        rk -= rk.mean(0)
        nr = np.sqrt((rk ** 2).sum(0)); nr[nr == 0] = 1.0; rk /= nr
        print(f"\n=== {cls}: {len(sub):,} records, {pd.Series(blocks).nunique():,} blocks ===",
              flush=True)

        allz = []
        for rep in range(NPERM):
            rng = np.random.default_rng(1000 + rep)
            clim = sb[rng.permutation(len(us))][six]
            r = np.clip(unitrank(clim).astype(np.float32) @ rk, -0.999999, 0.999999).astype(np.float64)
            tt = r * np.sqrt((npool - 2) / (1 - r ** 2))
            pv = 2 * tdist.sf(np.abs(tt), npool - 2)
            pV = np.clip(pd.Series(pv).rank(method="first").to_numpy() / len(pv), 1e-15, 1 - 1e-3)
            z = norm.ppf(1 - pV)
            d = pd.DataFrame({"b": blocks, "num": pq * z, "den": pq ** 2})
            g = d.groupby("b")
            res = pd.DataFrame({"n": g.size(), "Z": g["num"].sum() / np.sqrt(g["den"].sum())})
            res = res[res.n >= 2]
            s = res.sort_values("n")
            v = s.Z.rolling(ROLLER, min_periods=MINE).var()
            m = ~v.isnull()
            X = np.asarray(s.n.rolling(ROLLER, min_periods=MINE).mean()[m], float)
            sd = np.asarray(np.sqrt(v)[m], float)
            mn = np.asarray(s.Z.rolling(ROLLER, min_periods=MINE).mean()[m], float)
            t = res.n.to_numpy(float)
            sdp = IsotonicRegression(increasing=True, out_of_bounds="clip").fit(X, sd).predict(t)
            mnp = IsotonicRegression(increasing="auto", out_of_bounds="clip").fit(X, mn).predict(t)
            allz.append(((res.Z.to_numpy() - mnp) / sdp))
            nb = len(res)
        del rk

        Zs = np.concatenate(allz)
        alpha = 0.05 / nb
        zcrit = norm.ppf(1 - alpha)
        obs = float((Zs > zcrit).mean() * nb)

        # generalized-Pareto extrapolation of the permutation tail -> honest critical value
        u = np.quantile(Zs, TAIL_Q)
        exc = Zs[Zs > u] - u
        shape, loc, scale = genpareto.fit(exc, floc=0)
        # P(Z > u + y) = (1 - TAIL_Q) * SF_gpd(y)  ->  solve for the alpha quantile
        z_hon = u + genpareto.ppf(1 - alpha / (1 - TAIL_Q), shape, loc=0, scale=scale)

        print(f"  null standardized Z: n={len(Zs):,}  mean {Zs.mean():+.3f}  sd {Zs.std():.3f}")
        print(f"  EXCESS KURTOSIS {kurtosis(Zs):+.2f}   (Normal = 0)")
        for q in (0.99, 0.999, 0.9999):
            print(f"    q={q:<7}: empirical {np.quantile(Zs, q):6.3f}  vs Normal {norm.ppf(q):6.3f}")
        print(f"  nominal Bonferroni z_crit {zcrit:.3f} -> {obs:.2f} null hits/rep (nominal 0.05)")
        # NB heavy-tailedness RELATIVE TO NORMAL is established by the excess kurtosis and
        # the empirical-vs-Normal quantiles above -- NOT by the sign of xi. xi is the shape
        # of the GPD fitted to exceedances over u, so xi<0 means the tail decays FASTER than
        # exponential above u (bounded upper tail). That is good news: the extrapolation is
        # finite and well-behaved rather than power-law.
        print(f"  GPD-calibrated honest z_crit {z_hon:.3f}  (+{z_hon - zcrit:.3f} z; "
              f"GPD shape xi={shape:+.3f}, xi<0 => tail decays faster than exponential "
              f"above the {TAIL_Q} threshold, so this extrapolation is stable)", flush=True)

        rows.append(dict(cls=cls, blocks=nb, n_null=len(Zs), sd=Zs.std(),
                         excess_kurtosis=kurtosis(Zs),
                         q99=np.quantile(Zs, .99), q999=np.quantile(Zs, .999),
                         q9999=np.quantile(Zs, .9999),
                         z_nominal=zcrit, null_hits_at_nominal=obs,
                         gpd_xi=shape, z_honest=z_hon, z_shift=z_hon - zcrit))

    D = pd.DataFrame(rows)
    D.to_csv(f"{HERE}/null_tail_shape.csv", index=False)
    print("\n" + D.round(3).to_string(index=False))
    print(f"\nwrote {HERE}/null_tail_shape.csv")


if __name__ == "__main__":
    main()
