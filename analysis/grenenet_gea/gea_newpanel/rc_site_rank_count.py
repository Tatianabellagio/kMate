#!/usr/bin/env python
"""Site-level (honest 31-unit) RANK (Kendall-tau) + COUNT (binomial GLM) climate-GEA.

Alternatives to LFMM at the honest 31-site unit. For one class in {snp,nonsnp,sv}
and one climate axis, across the 31 SITES:

  site_af[s, v] = Delta-p[s, v] + p0[col[v]]        (from lfmm_site Y + af_store p0)
  Kendall : tau-b(site_af_v, climate) across 31 sites            (rank-based)
  Binomial: GLM [alt, ref] ~ const + z(climate), N_site genomes  (count-based)
            N_site = (sum total_flowers in site) * 2

p0: snp -> p0_snp.npy, {nonsnp,sv} -> p0_nonsnp.npy (sv col indexes the nonsnp p0).
Site order = np.sort(pools.site.unique()) == lfmm_site Y rows == env_site rows.

Reports per test:
  lambda_GC (genomic inflation): Kendall from chi2=chi2.isf(p,1); binomial from
            Wald z^2 (avoids p==0 underflow); lambda = median(chi2)/0.4549364.
  #Bonferroni blocks / #FDR blocks: clq0.9 blocks (blocks_clq09) with >=1 variant
            passing Bonferroni (0.05/n_tested) or BH-FDR q<0.05.
Tested family = site-MAF>=0.05 variants with a finite statistic.

SV class: writes annotated (genes) hit tables for both tests.

Outputs: analysis/grenenet_gea/gea_newpanel/results/rank_count/
  summary_{cls}_{axis}.json          lambda + block/variant counts (both tests)
  svhits_{test}_{cls}_{axis}.csv      (sv only) FDR/Bonferroni hits + genes

Usage (basic env has statsmodels):
  PY=.../envs/basic/bin/python
  $PY rc_site_rank_count.py --class sv --axis bio5 --threads 16
"""
from __future__ import annotations
import argparse, os, sys, json, time, warnings
import numpy as np
import pandas as pd
from multiprocessing import Pool
import statsmodels.api as sm
from scipy.stats import kendalltau, chi2
from sklearn.preprocessing import StandardScaler

sys.path.insert(0, "/global/scratch/users/tbellg/kmate/analysis/grenenet_gea")
sys.path.insert(0, "/global/scratch/users/tbellg/kmate/analysis/grenenet_gea/gea_newpanel")
from blocks_clq09 import assign_clq09_blocks
from lib import annotate_svs, load_genes

warnings.filterwarnings("ignore")

ROOT = "/global/scratch/users/tbellg/kmate/analysis/grenenet_gea"
LFMM = f"{ROOT}/gea_newpanel/lfmm_site"
ENVD = f"{ROOT}/gea_newpanel/env_site"
RECD = f"{ROOT}/phase1_replication/class_matrices"
OUTD = f"{ROOT}/gea_newpanel/rank_count"
os.makedirs(OUTD, exist_ok=True)

NSITES = 31
MAF_MIN = 0.05
CHI2_MED = chi2.ppf(0.5, 1)          # 0.4549364 -- median of a 1-df chi-square
P0 = {"snp": "p0_snp.npy", "nonsnp": "p0_nonsnp.npy", "sv": "p0_nonsnp.npy"}

_CLIM = None      # (nsite,) climate per site
_X = None         # (nsite, 2) const + z(climate)
_AFT = None       # (rec, nsite) site AF
_SUCC = None      # (rec, nsite) alt counts
_FAIL = None      # (rec, nsite) ref counts


def _chunk(rng):
    a, b = rng
    tau = np.full(b - a, np.nan); pk = np.full(b - a, np.nan)
    slope = np.full(b - a, np.nan); zb = np.full(b - a, np.nan); pb = np.full(b - a, np.nan)
    for i in range(a, b):
        y = _AFT[i]
        ok = np.isfinite(y)
        if ok.sum() < 3 or np.ptp(y[ok]) == 0:
            continue
        t, p = kendalltau(_CLIM[ok], y[ok])
        tau[i - a] = t; pk[i - a] = p
        s = _SUCC[i][ok]; f = _FAIL[i][ok]
        if s.sum() == 0 or f.sum() == 0:
            continue
        X = _X[ok]
        if np.ptp(X[:, 1]) == 0:
            continue
        try:
            res = sm.GLM(np.column_stack([s, f]), X, family=sm.families.Binomial()).fit()
            slope[i - a] = res.params[1]
            zb[i - a] = res.tvalues[1]         # Wald z (slope/bse)
            pb[i - a] = res.pvalues[1]
        except Exception:
            pass
    return a, tau, pk, slope, zb, pb


def bh_fdr(p):
    p = np.asarray(p, float); n = p.size
    order = np.argsort(p)
    ranked = p[order] * n / (np.arange(n) + 1)
    q = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n); out[order] = np.clip(q, 0, 1)
    return out


def lam_from_p(p):
    """lambda_GC from p-values via 1-df chi-square inversion."""
    p = np.clip(np.asarray(p, float), 1e-300, 1.0)
    c = chi2.isf(p, 1)
    return float(np.median(c) / CHI2_MED)


def lam_from_z(z):
    """lambda_GC from Wald z (chi2 = z^2), robust to p underflow."""
    z = np.asarray(z, float)
    return float(np.median(z * z) / CHI2_MED)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--class", dest="cls", required=True, choices=["snp", "nonsnp", "sv"])
    ap.add_argument("--axis", required=True)
    ap.add_argument("--threads", type=int, default=8)
    args = ap.parse_args()
    t0 = time.time()

    dims = open(f"{LFMM}/lfmm_{args.cls}_site_dims.txt").read().split()
    ns, nv = int(dims[0]), int(dims[1])
    assert ns == NSITES, (ns, NSITES)
    Y = np.fromfile(f"{LFMM}/lfmm_{args.cls}_site_Y.f64", dtype=np.float64).reshape(ns, nv)
    recs = pd.read_csv(f"{RECD}/{args.cls}_gen9.records.csv")
    assert len(recs) == nv, (len(recs), nv)
    p0 = np.load(f"{ROOT}/af_store/{P0[args.cls]}").astype(np.float64)
    col = recs["col"].to_numpy()
    site_af = Y + p0[col][None, :]                       # [31 x nv]
    del Y

    # site order and per-site diploid genome count (same np.sort(site) order as Y/env)
    pools = pd.read_csv(f"{RECD}/gen9.pools.csv")
    sites = np.sort(pools.site.unique())
    fl = pools.groupby("site").total_flowers.sum().loc[sites].to_numpy(float)
    N_site = np.round(fl * 2.0)

    clim = pd.read_csv(f"{ENVD}/env_site_{args.axis}.csv").iloc[:, 0].to_numpy(float)
    assert clim.size == NSITES
    z = StandardScaler().fit_transform(clim.reshape(-1, 1)).ravel()
    X = sm.add_constant(z, has_constant="add")           # [31 x 2]

    aft = np.ascontiguousarray(site_af.T)                # [rec x 31]
    finite = np.isfinite(aft)
    succ = np.where(finite, np.round(aft * N_site[None, :]), 0.0)
    fail = np.where(finite, np.round((1.0 - aft) * N_site[None, :]), 0.0)
    mean_af = np.nanmean(aft, axis=1)
    site_maf = np.minimum(mean_af, 1.0 - mean_af)

    print(f"[{args.cls}/{args.axis}] {ns} sites x {nv:,} rec; climate "
          f"{clim.min():.2f}..{clim.max():.2f}; N/site {int(N_site.min())}-{int(N_site.max())}; "
          f"MAF>=.05 {(site_maf>=MAF_MIN).sum():,}", flush=True)

    global _CLIM, _X, _AFT, _SUCC, _FAIL
    _CLIM = clim; _X = X; _AFT = aft; _SUCC = succ; _FAIL = fail

    step = 5000
    ranges = [(a, min(a + step, nv)) for a in range(0, nv, step)]
    tau = np.full(nv, np.nan); pk = np.full(nv, np.nan)
    slope = np.full(nv, np.nan); zb = np.full(nv, np.nan); pb = np.full(nv, np.nan)
    with Pool(args.threads) as pool:
        for a, t, p, sl, zz, q in pool.imap_unordered(_chunk, ranges):
            tau[a:a+len(t)] = t; pk[a:a+len(p)] = p
            slope[a:a+len(sl)] = sl; zb[a:a+len(zz)] = zz; pb[a:a+len(q)] = q
    print(f"  fits done ({time.time()-t0:.0f}s)", flush=True)

    blocks = assign_clq09_blocks(recs["chrom"].to_numpy(), recs["pos"].to_numpy())
    recs = recs.assign(site_maf=site_maf, block_clq09=blocks)
    genes = load_genes() if args.cls == "sv" else None

    maf_ok = site_maf >= MAF_MIN
    summary = {"cls": args.cls, "axis": args.axis, "nv": int(nv),
               "n_maf_ok": int(maf_ok.sum())}

    for test, stat, pv, zstat in [("kendall", tau, pk, None),
                                  ("binomial", slope, pb, zb)]:
        fin = maf_ok & np.isfinite(pv)
        n_test = int(fin.sum())
        pv_t = pv[fin]
        if zstat is not None:
            lam = lam_from_z(zstat[fin])
        else:
            lam = lam_from_p(pv_t)
        bonf_thr = 0.05 / max(n_test, 1)
        pass_bonf = fin & (pv < bonf_thr)
        q = np.full(nv, np.nan)
        q[fin] = bh_fdr(pv_t)
        pass_fdr = fin & (q < 0.05)
        nblk_bonf = recs.loc[pass_bonf, "block_clq09"].nunique()
        nblk_fdr = recs.loc[pass_fdr, "block_clq09"].nunique()
        summary[test] = dict(
            lambda_gc=round(lam, 4), n_tested=n_test,
            min_p=float(np.nanmin(pv_t)) if n_test else float("nan"),
            n_var_bonf=int(pass_bonf.sum()), n_var_fdr=int(pass_fdr.sum()),
            n_blocks_bonf=int(nblk_bonf), n_blocks_fdr=int(nblk_fdr),
            bonf_thr=float(bonf_thr))
        print(f"  {test:8s}: lambda={lam:.3f} nTest={n_test:,} minP={summary[test]['min_p']:.2e} "
              f"Bonf {pass_bonf.sum()}v/{nblk_bonf}blk  FDR {pass_fdr.sum()}v/{nblk_fdr}blk",
              flush=True)

        if args.cls == "sv":
            hit = recs[pass_fdr | pass_bonf].copy()
            hit["stat"] = stat[hit.index.to_numpy()]
            hit["pval"] = pv[hit.index.to_numpy()]
            hit["q"] = q[hit.index.to_numpy()]
            hit["pass_bonf"] = pass_bonf[hit.index.to_numpy()]
            hit["pass_fdr"] = pass_fdr[hit.index.to_numpy()]
            if len(hit):
                hit = annotate_svs(hit, flank=2000, genes=genes)
            hit.sort_values("pval").to_csv(
                f"{OUTD}/svhits_{test}_{args.cls}_{args.axis}.csv", index=False)

    with open(f"{OUTD}/summary_{args.cls}_{args.axis}.json", "w") as f:
        json.dump(summary, f, indent=2)
    print(f"  -> summary_{args.cls}_{args.axis}.json ({time.time()-t0:.0f}s)", flush=True)


if __name__ == "__main__":
    main()
