#!/usr/bin/env python
"""Permutation-null CALIBRATION test of the settled WZA regime (2026-07-28).

Question this answers
---------------------
`notebooks/cap_poly_decision.ipynb` settled the correction from fit quality:
tiling blocks, isotonic SD, isotonic_auto mean, NO CAP. Two residual risks could not
be closed by fit diagnostics alone:

  (1) With no cap, isotonic holds the SD FLAT above the rolling-window support edge.
      For the 8-15 largest blocks per class the flat clip is an extrapolation. If the
      true SD kept rising there, predicted SD would be too small -> p too small ->
      FALSE POSITIVES. It was directly measurable in only 2 of 8 blockdef x class
      cases (elsewhere <15 blocks sit past support, too few to estimate an SD).
  (2) The mean is now smoothed and clipped flat too.

Under a SIGNAL-FREE null there is nothing to find, so ANY block reaching significance
is a correction artifact. That is the decisive calibration check, and it is the same
argument that convicted deg-2 (see perm_null_sd.py, which established the null
SD-vs-n shape is monotone-up-then-plateau, isotonic-R2~=0.98).

Null construction: SITE-permute bio1 across the 31 GrENE-Net sites. This breaks the
genotype-climate association while preserving LD, pool/site structure and the
block-size distribution -- i.e. everything the SNP-number correction is estimated on.

`--level` selects the per-record test, which turned out to matter far more than the fit:
  pool  the production test -- rank-correlate AF across all 352 POOLS. But pools within
        a site share one climate value, so there are only 31 independent climates and
        n=352 df is pseudoreplicated. `_kendall_site_collapsed_test.py` measured this at
        6-7 orders of magnitude of p inflation for two example blocks.
  site  flower-weighted mean AF per site -> 31 points, n=31 df. Honest, and under
        site-permutation it is exactly exchangeable, so ANY residual excess of
        significant blocks here is the WZA aggregation's own doing, not
        pseudoreplication. Running both isolates the two causes.

Regimes compared (all on TILING blocks, so only the correction varies)
----------------------------------------------------------------------
  upstream_nocap   deg2 SD + deg2 mean,          no cap   <- upstream default (sample_snps=0)
  oldprod_cap      isotonic SD + interp mean,    cap      <- kMate production pre-2026-07-28
  proposed_nocap   isotonic SD + isotonic_auto,  no cap   <- SETTLED regime under test
  proposed_cap     isotonic SD + isotonic_auto,  cap      <- isolates the cap's contribution

Reported per regime per class, per replicate: fabricated p==0, NaN, Bonferroni hits and
BH-FDR hits (all should be ~0 under the null), plus how many of those hits are blocks
sitting PAST the support edge -- the targeted test of risk (1).

Run via sbatch (see perm_null_regime.sbatch). Writes perm_null_regime.csv + .npz.
"""
from __future__ import annotations
import os, sys
import numpy as np
import pandas as pd
from scipy.stats import norm, t as tdist
from sklearn.isotonic import IsotonicRegression

# self-locating: do NOT hardcode an absolute kmate path (the repo has moved once already)
HERE = os.path.dirname(os.path.abspath(__file__))
GEA = os.path.dirname(HERE)
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), "blocks"))
import lib                      # noqa: E402
import blocks_tiling as bt      # noqa: E402

CM = f"{GEA}/r2_gea_nonsnp/phase1_replication/results/class_matrices"
OUT = f"{HERE}/perm_null_regime"

CLASSES = ["snp", "sv", "smallindel", "nonsnp"]
# the caps the pre-2026-07-28 tiling arms used, for the cap-vs-nocap contrast only
CAP = {"snp": 2000, "sv": 700, "smallindel": 700, "nonsnp": 700}
R2 = 0.9
NPERM = 20
MAF_MIN, ROLLER, MINE, MIN_SNPS = 0.05, 50, 40, 2
RESAMPLES = 20      # wza_script uses 100; 20 is plenty for a null and 5x cheaper

REGIMES = {
    "upstream_nocap": dict(sd="deg2", mean="deg2", cap=False),
    "oldprod_cap":    dict(sd="isotonic", mean="interp", cap=True),
    "proposed_nocap": dict(sd="isotonic", mean="isotonic_auto", cap=False),
    "proposed_cap":   dict(sd="isotonic", mean="isotonic_auto", cap=True),
}


def unitrank(v):
    r = pd.Series(v).rank().to_numpy().astype(np.float64)
    r = r - r.mean()
    n = np.linalg.norm(r)
    return r / n if n else r


def fit_curve(Xs, y, target, kind):
    if kind == "isotonic":
        return IsotonicRegression(increasing=True, out_of_bounds="clip").fit(Xs, y).predict(target)
    if kind == "isotonic_auto":
        return IsotonicRegression(increasing="auto", out_of_bounds="clip").fit(Xs, y).predict(target)
    if kind == "interp":
        return np.interp(target, Xs, y)
    return np.poly1d(np.polyfit(Xs, y, int(kind[3:])))(target)


def bh_fdr(p):
    p = np.asarray(p, float)
    ok = np.isfinite(p)
    q = np.full(len(p), np.nan)
    v = p[ok]
    o = np.argsort(v)
    n = len(v)
    adj = np.minimum.accumulate((v[o] * n / np.arange(1, n + 1))[::-1])[::-1]
    out = np.empty(n)
    out[o] = np.clip(adj, 0, 1)
    q[ok] = out
    return q


def block_stat(pval, pq, blocks, cap, seed):
    """Production block weighted-Z: genome-wide rank-transform of p (wza_script.py:167),
    z = ppf(1-p), MAF-weighted Stouffer per block. With `cap`, oversized blocks are
    replaced by the resample-mean and recorded as SNPs=cap (wza_script.py:204-207)."""
    N = len(pval)
    ranks = pd.Series(pval).rank(method="first").to_numpy()
    pV = np.clip(ranks / N, 1e-15, 1 - 1e-3)
    z = norm.ppf(1 - pV)
    d = pd.DataFrame({"b": blocks, "num": pq * z, "den": pq ** 2})
    g = d.groupby("b")
    res = pd.DataFrame({"SNPs_raw": g.size(), "Z": g["num"].sum() / np.sqrt(g["den"].sum())})
    res = res[res.SNPs_raw >= MIN_SNPS].copy()
    res["SNPs"] = res.SNPs_raw.astype(float)
    if cap:
        big = res.index[res.SNPs_raw > cap]
        if len(big):
            gi = pd.Series(np.arange(len(blocks))).groupby(blocks).indices
            rng = np.random.default_rng(seed)
            for b in big:
                ix = gi[b]
                vals = [(lambda s: pq[s] @ z[s] / np.sqrt(pq[s] @ pq[s]))(
                            rng.choice(ix, cap, replace=False)) for _ in range(RESAMPLES)]
                res.loc[b, ["Z", "SNPs"]] = [np.mean(vals), float(cap)]
    return res


def correct(res, sd_kind, mean_kind):
    s = res.sort_values("SNPs")
    v = s.Z.rolling(ROLLER, min_periods=MINE).var()
    m = ~v.isnull()
    Xs = np.asarray(s.SNPs.rolling(ROLLER, min_periods=MINE).mean()[m], float)
    sd = np.asarray(np.sqrt(v)[m], float)
    mn = np.asarray(s.Z.rolling(ROLLER, min_periods=MINE).mean()[m], float)
    tgt = res.SNPs.to_numpy(float)
    with np.errstate(invalid="ignore"):
        p = 1 - norm.cdf(res.Z.to_numpy(float),
                         loc=fit_curve(Xs, mn, tgt, mean_kind),
                         scale=fit_curve(Xs, sd, tgt, sd_kind))
    return p, float(Xs.max())


def main():
    global CLASSES, NPERM, OUT
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--classes", nargs="+", default=CLASSES, choices=CLASSES)
    ap.add_argument("--nperm", type=int, default=NPERM)
    ap.add_argument("--level", default="pool", choices=["pool", "site"],
                    help="per-record test level: 'pool' = production (n=352, "
                         "pseudoreplicated) or 'site' = flower-weighted site means (n=31, honest)")
    ap.add_argument("--tag", default="", help="suffix for the output files (e.g. _smoke)")
    a = ap.parse_args()
    CLASSES, NPERM, OUT = a.classes, a.nperm, OUT + a.tag
    LEVEL = a.level

    pools = pd.read_csv(f"{CM}/gen9.pools.csv")
    site = pools["site"].to_numpy()
    usites = np.unique(site)
    site_ix = np.array([{s: i for i, s in enumerate(usites)}[s] for s in site])
    site_bio1 = pools.groupby("site")["bio1"].first().reindex(usites).to_numpy()
    # flower weights for the site collapse (same aggregation as _kendall_site_collapsed_test.py)
    w = pools["total_flowers"].to_numpy(float)
    w = np.where(np.isfinite(w) & (w > 0), w, 1.0)
    # [sites x pools] weight matrix, rows summing to 1 -> AF_site = W @ AF_pool
    W = np.zeros((len(usites), len(pools)))
    W[site_ix, np.arange(len(pools))] = w
    W /= W.sum(1, keepdims=True)
    ntest = len(pools) if LEVEL == "pool" else len(usites)
    print(f"{len(pools)} pools, {len(usites)} sites, {NPERM} permutations, r2={R2} TILING "
          f"blocks | LEVEL={LEVEL} -> per-record test n={ntest}", flush=True)

    rows = []
    for cls in CLASSES:
        print(f"\n=== {cls} ===", flush=True)
        af = np.load(f"{CM}/{cls}_gen9_af.npy")
        rec = pd.read_csv(f"{CM}/{cls}_gen9.records.csv")
        # TILING assignment -- identical to reblock_blockdef.py --how tiling (0% dropped)
        blk = bt.merge_small_blocks(
            bt.assign_tiling(rec["chrom"].to_numpy(str), rec["pos"].to_numpy(np.int64), r2=R2))
        rec["block"] = blk
        keep = (rec["block"].to_numpy() != "") & (rec["maf"].to_numpy() > MAF_MIN)
        sub = rec[keep]
        blocks = sub["block"].to_numpy().astype(str)
        pq = (sub["maf"] * (1 - sub["maf"])).to_numpy()
        afx = af[:, sub.index.to_numpy()].astype(np.float64)
        print(f"  records kept {len(sub):,}/{len(rec):,} ({100*len(sub)/len(rec):.1f}%) | "
              f"{pd.Series(blocks).nunique():,} blocks | max {pd.Series(blocks).value_counts().max():,}",
              flush=True)

        if LEVEL == "site":
            afx = W @ afx                      # [sites x recs] flower-weighted site mean AF
        rk = afx.argsort(0).argsort(0).astype(np.float64) + 1.0
        rk -= rk.mean(0)
        nrm = np.sqrt((rk ** 2).sum(0)); nrm[nrm == 0] = 1.0; rk /= nrm
        del afx

        def spearman_p(clim):
            r = np.clip(unitrank(clim) @ rk, -0.999999, 0.999999)
            tt = r * np.sqrt((ntest - 2) / (1 - r ** 2))
            return 2 * tdist.sf(np.abs(tt), ntest - 2)

        for rep in range(NPERM):
            rng = np.random.default_rng(1000 + rep)
            perm = site_bio1[rng.permutation(len(usites))]            # climate shuffled ACROSS SITES
            clim = perm if LEVEL == "site" else perm[site_ix]
            pv = spearman_p(clim)
            for name, cfg in REGIMES.items():
                res = block_stat(pv, pq, blocks, CAP[cls] if cfg["cap"] else None, seed=rep)
                p, sup = correct(res, cfg["sd"], cfg["mean"])
                n = len(p)
                q = bh_fdr(p)
                bonf = np.isfinite(p) & (p < 0.05 / n)
                past = res.SNPs_raw.to_numpy() > sup
                rows.append(dict(cls=cls, level=LEVEL, regime=name, rep=rep,
                                 blocks=n, support_to=sup,
                                 p0=int((p == 0).sum()), nan=int(np.isnan(p).sum()),
                                 bonf=int(bonf.sum()), fdr=int((q < 0.05).sum()),
                                 bonf_past_support=int((bonf & past).sum()),
                                 n_past_support=int(past.sum())))
            if rep == 0 or (rep + 1) % 5 == 0:
                d = pd.DataFrame([r for r in rows if r["cls"] == cls])
                print(f"  rep {rep+1}/{NPERM} cumulative means: " + " | ".join(
                    f"{nm}: bonf={d[d.regime==nm].bonf.mean():.1f} p0={d[d.regime==nm].p0.mean():.1f}"
                    for nm in REGIMES), flush=True)

    D = pd.DataFrame(rows)
    D.to_csv(f"{OUT}.csv", index=False)
    print(f"\nwrote {OUT}.csv ({len(D)} rows)\n")
    print("PER-REPLICATE MEANS UNDER THE SIGNAL-FREE NULL (all should be ~0):\n")
    piv = D.groupby(["cls", "regime"])[["p0", "nan", "bonf", "fdr", "bonf_past_support"]].mean()
    print(piv.round(2).to_string())
    print("\nWORST replicate per cls x regime (bonf):\n")
    print(D.groupby(["cls", "regime"]).bonf.max().unstack().to_string())


if __name__ == "__main__":
    main()
