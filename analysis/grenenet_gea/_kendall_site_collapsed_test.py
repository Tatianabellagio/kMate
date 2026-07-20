#!/usr/bin/env python
"""Prototype/validation: SITE-COLLAPSED Kendall-tau (fix for pool-in-site
pseudoreplication in the phase1 Kendall model).

WHY: `wza_investigation/RESULTS.md` "BACK TO BASICS" section showed, for 2 example
blocks (CAM5 2_1265, Chr4 CRK 4_2519), that per-SNP Kendall computed across ALL
POOLS (326, gen1) gives p as extreme as 8e-9/1.2e-8, but the SAME correlation
computed across the 31 INDEPENDENT SITES (collapsing pools within a site to one
flower-weighted mean AF) gives an honest p of only ~0.015/0.012 -- a 6-7 order-of-
magnitude inflation from pseudoreplication, on top of (and separate from) the
WZA SNP-count issue. This script generalizes that fix genome-wide (gen9, bio1,
snp + nonsnp) and validates it against the production pool-level Kendall.

Site-level AF = flower-weighted mean of pool AF across pools in the same site
(same aggregation as wza_investigation/raw_signal.py's `sdf.groupby("site").mean()`,
except that used an unweighted site mean of pool AFs; here we additionally weight
by total_flowers, the natural pool-seq weight -- both collapse pools to one point
per site so the pseudoreplication fix is identical; flower-weighting is the more
principled aggregation and used throughout this pipeline, e.g. run_betabinom_latent).

DOES NOT touch/overwrite anything under clq90/{kendall,wza,wza_in}/ or
clq90/compare_clq90.csv (frozen production). Writes only to:
  analysis/grenenet_gea/phase1_replication/results/clq90/kendall_fix_test/site_collapsed/

Usage:
  PY=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python
  $PY analysis/grenenet_gea/_kendall_site_collapsed_test.py
"""
from __future__ import annotations
import os, sys, time
import numpy as np
import pandas as pd
from multiprocessing import Pool
from scipy.stats import kendalltau, chi2

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lib

CMDIR = f"{lib.GEA}/phase1_replication/results/class_matrices"
PRODDIR = f"{lib.GEA}/phase1_replication/results/clq90/wza_in"       # frozen production, READ-ONLY
OUT = f"{lib.GEA}/phase1_replication/results/clq90/kendall_fix_test/site_collapsed"
GEN = 9
CLIM = "bio1"
CLASSES = ["snp", "nonsnp"]
THREADS = 4
MAF_FLOOR = 0.05

# spot-check windows (old hapFIRE-block coords, from wza_investigation/RESULTS.md
# "BACK TO BASICS" section; re-identified against kendall_snp_gen1_bio1.csv):
#   CAM5    2_1265: Chr2:11,533,904-11,534,244
#   Chr4CRK 4_2519: Chr4:12,133,539-12,168,827
SPOTCHECK = {
    "CAM5 (2_1265)": ("Chr2", 11_530_000, 11_540_000),
    "Chr4 CRK (4_2519)": ("Chr4", 12_130_000, 12_172_000),
}

_X = None      # site-level climate vector [n_sites]
_AFT = None    # site-level AF, transposed [records x n_sites]


def _chunk(rng):
    a, b = rng
    tau = np.full(b - a, np.nan); pv = np.full(b - a, np.nan)
    for i in range(a, b):
        y = _AFT[i]
        ok = np.isfinite(y)
        if ok.sum() >= 3 and np.ptp(y[ok]) > 0:
            t, p = kendalltau(_X[ok], y[ok])
            tau[i - a] = t; pv[i - a] = p
    return a, tau, pv


def gif_proxy(p):
    p = p[np.isfinite(p) & (p > 0)]
    if len(p) == 0:
        return np.nan
    return np.median(chi2.isf(p, 1)) / chi2.ppf(0.5, 1)


def site_collapse(af, site_codes, sites, weights):
    """af: [pools x records] float32. Returns site-level AF [n_sites x records]."""
    n_sites = len(sites)
    n_pools = af.shape[0]
    W = np.zeros((n_sites, n_pools), dtype=np.float32)
    for s_i, s in enumerate(sites):
        m = site_codes == s
        W[s_i, m] = weights[m]

    M = np.isfinite(af)
    AFf = np.where(M, af, 0.0).astype(np.float32)
    Mf = M.astype(np.float32)
    num = W @ AFf          # [n_sites x records]
    den = W @ Mf
    with np.errstate(invalid="ignore", divide="ignore"):
        site_af = num / den
    site_af[den == 0] = np.nan
    return site_af


def main():
    os.makedirs(OUT, exist_ok=True)
    pools = pd.read_csv(f"{CMDIR}/gen{GEN}.pools.csv")
    site_codes = pools["site"].to_numpy()
    weights = pools["total_flowers"].to_numpy(float)
    sites = np.sort(pools["site"].unique())
    n_sites = len(sites)

    # site-level climate: bio1 is constant within a site -> take the (unique) mean
    site_bio1 = pools.groupby("site")[CLIM].mean().reindex(sites).to_numpy(float)
    nuniq = pools.groupby("site")[CLIM].nunique()
    assert (nuniq <= 1).all(), "bio1 not constant within a site -- unexpected"
    print(f"n_pools={len(pools)}  n_sites={n_sites}  climate={CLIM}", flush=True)

    global _X, _AFT
    summary = []
    for cls in CLASSES:
        t0 = time.time()
        recs = pd.read_csv(f"{CMDIR}/{cls}_gen{GEN}.records.csv")
        af = np.load(f"{CMDIR}/{cls}_gen{GEN}_af.npy")   # [pools x records]
        assert af.shape[0] == len(pools) and af.shape[1] == len(recs)
        n_rec = af.shape[1]
        print(f"[{cls}] {af.shape[0]} pools x {n_rec:,} records -> collapsing to "
              f"{n_sites} sites ...", flush=True)

        site_af = site_collapse(af, site_codes, sites, weights)   # [sites x records]
        del af
        n_used = np.isfinite(site_af).sum(axis=0)
        print(f"[{cls}] collapse done ({time.time()-t0:.0f}s); "
              f"median sites-used/record = {np.median(n_used):.0f}", flush=True)

        _X = site_bio1
        _AFT = np.ascontiguousarray(site_af.T)   # [records x sites]
        del site_af

        step = 50000
        ranges = [(a, min(a + step, n_rec)) for a in range(0, n_rec, step)]
        tau = np.empty(n_rec, np.float64); pv = np.empty(n_rec, np.float64)
        t1 = time.time()
        with Pool(THREADS) as pool:
            for a, t, p in pool.imap_unordered(_chunk, ranges):
                tau[a:a + len(t)] = t; pv[a:a + len(p)] = p
        print(f"[{cls}] kendall done ({time.time()-t1:.0f}s)", flush=True)

        recs = recs.assign(tau=tau, pval=pv)
        recs = recs.rename(columns={"maf": "MAF"})
        out = f"{OUT}/kendall_site_{cls}_gen{GEN}_{CLIM}.csv"
        cols = (["hap_id"] if "hap_id" in recs.columns else []) + \
               ["chrom", "pos", "ref_len", "alt_len", "MAF", "block", "tau", "pval"]
        recs[cols].to_csv(out, index=False)
        fin = np.isfinite(pv)
        print(f"[{cls}] tested {int(fin.sum()):,}/{n_rec:,} | p<0.05: "
              f"{int((pv[fin] < 0.05).sum()):,} | min p={np.nanmin(pv):.2e} -> {out}",
              flush=True)

        # ---------------- diagnostics vs production pool-level kendall -------
        # NOTE: wza_in/kendall_{cls}_gen9_bio1.csv is ALREADY MAF>=0.05-floored
        # (1.2M/1.99M snp rows) -- not the same row set/order as the full
        # class_matrices catalog, so align by (chrom,pos,ref_len,alt_len).
        mine = recs.assign(pval_site=pv, MAF_site=recs["MAF"].to_numpy(float))
        prod = pd.read_csv(f"{PRODDIR}/kendall_{cls}_gen{GEN}_{CLIM}.csv",
                            usecols=["chrom", "pos", "ref_len", "alt_len", "MAF", "pval"])
        prod = prod.rename(columns={"MAF": "MAF_pool", "pval": "pval_pool"})
        merged = mine.merge(prod, on=["chrom", "pos", "ref_len", "alt_len"], how="inner")
        print(f"[{cls}] matched {len(merged):,}/{len(prod):,} production rows "
              f"(of {n_rec:,} in full catalog)", flush=True)

        keep = merged["MAF_site"].to_numpy(float) >= MAF_FLOOR
        p_site = np.where(keep, merged["pval_site"].to_numpy(float), np.nan)
        p_pool = np.where(keep, merged["pval_pool"].to_numpy(float), np.nan)

        gif_site = gif_proxy(p_site)
        gif_pool = gif_proxy(p_pool)
        n_ok_site = np.isfinite(p_site).sum()
        n_ok_pool = np.isfinite(p_pool).sum()
        frac05_site = np.nanmean(p_site < 0.05)
        frac05_pool = np.nanmean(p_pool < 0.05)
        minp_site = np.nanmin(p_site)
        minp_pool = np.nanmin(p_pool)

        print(f"\n[{cls}] MAF>={MAF_FLOOR} diagnostics (n_site={n_ok_site:,}, "
              f"n_pool={n_ok_pool:,}):")
        print(f"  GIF proxy:      site-collapsed={gif_site:.3f}   "
              f"pool-level(prod)={gif_pool:.3f}")
        print(f"  frac p<0.05:    site-collapsed={frac05_site:.4f}   "
              f"pool-level(prod)={frac05_pool:.4f}")
        print(f"  min p:          site-collapsed={minp_site:.2e}   "
              f"pool-level(prod)={minp_pool:.2e}\n", flush=True)

        summary.append(dict(cls=cls, n_rec=n_rec, n_ok_site=int(n_ok_site),
                             n_ok_pool=int(n_ok_pool), gif_site=gif_site,
                             gif_pool=gif_pool, frac_p05_site=frac05_site,
                             frac_p05_pool=frac05_pool, minp_site=minp_site,
                             minp_pool=minp_pool))

        # ---------------- spot-check CAM5 / Chr4 CRK --------------------------
        mchrom = merged["chrom"].to_numpy()
        mpos = merged["pos"].to_numpy()
        for name, (c, lo, hi) in SPOTCHECK.items():
            sel = (mchrom == c) & (mpos >= lo) & (mpos <= hi)
            if sel.sum() == 0:
                continue
            sub_site_p = merged["pval_site"].to_numpy(float)[sel]
            sub_pool_p = merged["pval_pool"].to_numpy(float)[sel]
            fin_s = np.isfinite(sub_site_p)
            fin_p = np.isfinite(sub_pool_p)
            print(f"[{cls}] spot-check {name} ({c}:{lo}-{hi}, n={sel.sum()}): "
                  f"site-collapsed min p={np.nanmin(sub_site_p) if fin_s.any() else np.nan:.2e}, "
                  f"median p={np.nanmedian(sub_site_p[fin_s]) if fin_s.any() else np.nan:.2e}  |  "
                  f"production pool-level min p={np.nanmin(sub_pool_p) if fin_p.any() else np.nan:.2e}",
                  flush=True)

        _X = None; _AFT = None

    summ = pd.DataFrame(summary)
    summ.to_csv(f"{OUT}/diagnostics_summary.csv", index=False)
    print("\n=== SUMMARY ===")
    print(summ.to_string(index=False))
    print(f"\nwrote {OUT}/diagnostics_summary.csv")


if __name__ == "__main__":
    main()
