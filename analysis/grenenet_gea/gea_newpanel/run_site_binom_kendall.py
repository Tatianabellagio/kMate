#!/usr/bin/env python
"""Site-level (honest 31-unit) binomial + Kendall climate-GEA, SNP vs non-SNP.

The plot-level (355-pool) binomial/Kendall are massively inflated: 355 pools carry
only 31 independent climate values (pseudoreplication) and the binomial treats
AF x flowers x 2 as that many independent draws (over-precision). The honest unit
is the 31 flower-weighted SITE means — one climate value = one observation — exactly
the unit the LFMM site-PC1 panel uses. This re-runs BOTH raw phase-1 tests on it.

For each class in {snp, nonsnp}, aggregate the 355-pool class matrix to 31 sites:
  site_af[s] = sum_p (flowers_p / sum flowers_s) * AF_p          flower-weighted mean
Then per record, across the 31 sites:
  Kendall : tau-b(site_af, climate)                              (levels == Delta p ranks)
  Binomial: GLM [succ, fail] ~ const + z(climate), with
            N_site  = (sum flowers in site) * 2                  diploid genomes / site
            succ    = round(site_af * N_site);  fail = round((1-site_af) * N_site)
Site-level MAF = min(mean_s site_af, 1-mean_s site_af); records filtered MAF>=--maf
in the notebook (kept here so the same file serves both classes).

Reads build_class_matrices.py outputs (same as run_kendall/run_binomial):
  {class}_gen{g}_af.npy   [355 x rec]   {class}_gen{g}.records.csv   gen{g}.pools.csv

Output (--out): {test}_{class}_site_{climate}.csv  per record:
  chrom,pos,ref_len,alt_len,MAF,site_maf,block,{tau|slope},pval

Usage (basic env has statsmodels):
  PY=.../envs/basic/bin/python
  $PY run_site_binom_kendall.py --class snp --gen 9 --climate bio1 --threads 16
"""
from __future__ import annotations
import argparse, os, sys, warnings
import numpy as np
import pandas as pd
from multiprocessing import Pool
import statsmodels.api as sm
from scipy.stats import kendalltau
from sklearn.preprocessing import StandardScaler
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import lib

warnings.filterwarnings("ignore")

_CLIM = None      # (n_site,) climate per site
_X = None         # (n_site, 2) const + z(climate)
_AFT = None       # (rec, n_site) site-mean AF
_SUCC = None      # (rec, n_site) alt counts
_FAIL = None      # (rec, n_site) ref counts


def _chunk(rng):
    a, b = rng
    tau = np.full(b - a, np.nan); pk = np.full(b - a, np.nan)
    slope = np.full(b - a, np.nan); pb = np.full(b - a, np.nan)
    for i in range(a, b):
        y = _AFT[i]
        ok = np.isfinite(y)
        if ok.sum() < 3 or np.ptp(y[ok]) == 0:          # need variation + >=3 sites
            continue
        # Kendall tau-b (climate vs site-mean AF)
        t, p = kendalltau(_CLIM[ok], y[ok])
        tau[i - a] = t; pk[i - a] = p
        # Binomial GLM on site-aggregated counts
        s = _SUCC[i][ok]; f = _FAIL[i][ok]
        if s.sum() == 0 or f.sum() == 0:
            continue
        X = _X[ok]
        if np.ptp(X[:, 1]) == 0:
            continue
        try:
            res = sm.GLM(np.column_stack([s, f]), X, family=sm.families.Binomial()).fit()
            slope[i - a] = res.params[1]; pb[i - a] = res.pvalues[1]
        except Exception:
            pass
    return a, tau, pk, slope, pb


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--class", dest="cls", required=True,
                    choices=["snp", "sv", "smallindel", "nonsnp"])
    ap.add_argument("--gen", type=int, default=9)
    ap.add_argument("--climate", default="bio1")
    ap.add_argument("--cmdir", default=f"{lib.GEA}/phase1_replication/class_matrices")
    ap.add_argument("--out", default=f"{lib.GEA}/gea_newpanel")
    ap.add_argument("--threads", type=int, default=8)
    args = ap.parse_args()

    recs = pd.read_csv(f"{args.cmdir}/{args.cls}_gen{args.gen}.records.csv")
    pools = pd.read_csv(f"{args.cmdir}/gen{args.gen}.pools.csv")
    af = np.load(f"{args.cmdir}/{args.cls}_gen{args.gen}_af.npy")    # [pools x rec]
    assert af.shape[0] == len(pools) and af.shape[1] == len(recs), (af.shape, len(pools), len(recs))

    sites = np.sort(pools.site.unique())
    site_of = pools.site.to_numpy()
    w = pools.total_flowers.to_numpy(float)
    # flower-weighted site-mean AF [n_site x rec] and per-site diploid genome count
    n_site = len(sites)
    site_af = np.empty((n_site, af.shape[1]), dtype=np.float64)
    N_site = np.empty(n_site, dtype=np.float64)
    for i, s in enumerate(sites):
        m = site_of == s
        ws = w[m]; N_site[i] = ws.sum() * 2.0
        wn = ws / ws.sum()
        site_af[i] = wn @ af[m].astype(np.float64)
    del af

    clim = (pools.drop_duplicates("site").set_index("site").loc[sites, args.climate]
            .to_numpy(float))
    z = StandardScaler().fit_transform(clim.reshape(-1, 1)).ravel()
    X = sm.add_constant(z, has_constant="add")                      # [n_site x 2]

    aft = np.ascontiguousarray(site_af.T)                          # [rec x n_site]
    finite = np.isfinite(aft)
    succ = np.where(finite, np.round(aft * N_site[None, :]), 0.0)
    fail = np.where(finite, np.round((1.0 - aft) * N_site[None, :]), 0.0)

    site_maf = np.minimum(np.nanmean(aft, axis=1), 1.0 - np.nanmean(aft, axis=1))

    print(f"{args.cls} gen{args.gen}: {n_site} sites x {len(recs):,} records vs "
          f"{args.climate} (range {clim.min():.1f}-{clim.max():.1f}); "
          f"N(genomes)/site {int(N_site.min())}-{int(N_site.max())}", flush=True)

    global _CLIM, _X, _AFT, _SUCC, _FAIL
    _CLIM = clim; _X = X; _AFT = aft; _SUCC = succ; _FAIL = fail

    n_rec = len(recs); step = 5000
    ranges = [(a, min(a + step, n_rec)) for a in range(0, n_rec, step)]
    tau = np.full(n_rec, np.nan); pk = np.full(n_rec, np.nan)
    slope = np.full(n_rec, np.nan); pb = np.full(n_rec, np.nan)
    with Pool(args.threads) as pool:
        for a, t, p, sl, q in pool.imap_unordered(_chunk, ranges):
            tau[a:a + len(t)] = t; pk[a:a + len(p)] = p
            slope[a:a + len(sl)] = sl; pb[a:a + len(q)] = q

    os.makedirs(f"{args.out}/kendall_site", exist_ok=True)
    os.makedirs(f"{args.out}/binomial_site", exist_ok=True)
    base = recs.rename(columns={"maf": "MAF"}).assign(site_maf=site_maf)
    kcols = ["chrom", "pos", "ref_len", "alt_len", "MAF", "site_maf", "block"]
    kout = f"{args.out}/kendall_site/kendall_{args.cls}_site_{args.climate}.csv"
    base.assign(tau=tau, pval=pk)[kcols + ["tau", "pval"]].to_csv(kout, index=False)
    bout = f"{args.out}/binomial_site/binomial_{args.cls}_site_{args.climate}.csv"
    base.assign(slope=slope, pval=pb)[kcols + ["slope", "pval"]].to_csv(bout, index=False)

    for nm, pv in [("kendall", pk), ("binomial", pb)]:
        fin = np.isfinite(pv)
        mm = np.isfinite(pv) & (site_maf >= 0.05)
        print(f"  {nm}: tested {int(fin.sum()):,}/{n_rec:,} | "
              f"min p={np.nanmin(pv) if fin.any() else float('nan'):.2e} | "
              f"site-MAF>=.05 p<1e-5: {int((pv[mm] < 1e-5).sum()):,}", flush=True)
    print(f"  -> {kout}\n  -> {bout}", flush=True)


if __name__ == "__main__":
    main()
