#!/usr/bin/env python
"""Phase-1 binomial-regression GEA on a filtered class matrix (one class, one gen, one bioclim).

Mirrors the phase-1 `binomial_regression/.../run_partition_binomial_reg_last_gen.py`:
for each kept record, fit a per-SNP binomial GLM

    [successes, failures] ~ const + scaled(climate)

with response = focal-allele (alt) counts vs reference counts, predictor = the
standardized per-pool bioclim value across the generation's site_gen_plot pools.
The env-coefficient slope and its Wald p-value are the per-record stats -> WZA.

Allele counts follow phase-1 exactly (create_partitions_allele_freq.generate_allele_counts):
    N         = total_flowers * 2          (diploid genomes sampled in the pool)
    successes = round(AF       * N)        ("minor" in phase-1, here the alt allele)
    failures  = round((1 - AF) * N)        ("major")
Pools within a site share a climate value (ties) and are pseudoreplicated -- exactly
as in phase-1; structure correction is the job of LFMM / WZA, not this step.

Reads the outputs of build_class_matrices.py (same as run_kendall.py):
  {class}_gen{g}_af.npy        [pools x records]
  {class}_gen{g}.records.csv   record meta (incl. maf, block)
  gen{g}.pools.csv             pool rowmeta with total_flowers + bio1..bio19

Output (--out):
  binomial_{class}_gen{g}_{climate}.csv   per-record: chrom,pos,ref_len,alt_len,
                                          MAF,block,slope,pval  (-> WZA input)
Columns named so wza_script.py consumes it directly: needs `MAF`, a p-value
summary stat, and a window column (`block`).

Usage:
  PY=.../kmate/bin/python
  $PY run_binomial.py --class sv --gen 1 --climate bio1
"""
from __future__ import annotations
import argparse, os, sys, warnings
import numpy as np
import pandas as pd
from multiprocessing import Pool
import statsmodels.api as sm
from sklearn.preprocessing import StandardScaler
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import lib

_X = None       # design matrix [pools x 2]: const + scaled climate
_SUCC = None    # successes (alt counts) [records x pools]
_FAIL = None    # failures  (ref counts) [records x pools]
_W = None       # per-pool var_weights (None = ordinary binomial)

# WZA's z-transform (norm.ppf(1-p)) breaks for p<~1e-16 (float64 1-p precision) and
# NaN/0 poison the aggregation, so floor here at a tiny positive value.
P_FLOOR = 1e-300

warnings.filterwarnings("ignore")          # statsmodels PerfectSeparation / convergence noise


def _chunk(rng):
    a, b = rng
    slope = np.full(b - a, np.nan)
    pv = np.full(b - a, np.nan)
    for i in range(a, b):
        succ = _SUCC[i]
        fail = _FAIL[i]
        n = succ + fail
        ok = n > 0                                   # pool observed (finite AF -> count)
        # need >=3 pools, allele segregating across them, and climate variation
        if ok.sum() < 3:
            continue
        s = succ[ok]; f = fail[ok]
        if (s.sum() == 0) or (f.sum() == 0):         # monomorphic across observed pools
            continue
        X = _X[ok]
        if np.ptp(X[:, 1]) == 0:                     # no climate gradient among observed pools
            continue
        y = np.column_stack([s, f])
        # var_weights is a GLM *constructor* arg (NOT a .fit() arg — .fit silently
        # ignores unknown kwargs, leaving the fit unweighted).
        kw = {} if _W is None else {"var_weights": _W[ok]}
        try:
            res = sm.GLM(y, X, family=sm.families.Binomial(), **kw).fit()
            slope[i - a] = res.params[1]
            p = res.pvalues[1]
            pv[i - a] = max(p, P_FLOOR) if np.isfinite(p) else np.nan
        except Exception:
            pass
    return a, slope, pv


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--class", dest="cls", required=True,
                    choices=["snp", "sv", "smallindel", "nonsnp"])
    ap.add_argument("--gen", type=int, required=True)
    ap.add_argument("--climate", default="bio1")
    ap.add_argument("--cmdir", default=f"{lib.GEA}/r2_gea_nonsnp/phase1_replication/results/class_matrices")
    ap.add_argument("--out", default=f"{lib.GEA}/r2_gea_nonsnp/phase1_replication/results/binomial")
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--effective-n", dest="eff_n", default="none",
                    choices=["none", "site"],
                    help="deflate the pool binomial's over-precision at source. "
                         "'none'=phase-1 (N=flowers*2, pools independent); "
                         "'site'=var_weights 1/pools-per-site so the information "
                         "scales to the ~31 independent climates, not the 355 "
                         "pseudoreplicated pools (fixes the p-underflow: p<0.05 "
                         "84%->48%, p<1e-15 41%->1.4%, 0 exact-0/NaN).")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    recs = pd.read_csv(f"{args.cmdir}/{args.cls}_gen{args.gen}.records.csv")
    pools = pd.read_csv(f"{args.cmdir}/gen{args.gen}.pools.csv")
    clim = pools[args.climate].to_numpy(float)
    nflowers = pools["total_flowers"].to_numpy(float)
    N = np.round(nflowers * 2).astype(np.float64)            # diploid genomes per pool
    af = np.load(f"{args.cmdir}/{args.cls}_gen{args.gen}_af.npy")   # [pools x rec]
    assert af.shape[0] == len(clim), (af.shape, len(clim))
    assert af.shape[1] == len(recs), (af.shape, len(recs))
    n_rec = af.shape[1]

    # standardized climate design matrix [pools x 2] (const + z(climate)) -- phase-1 StandardScaler
    z = StandardScaler().fit_transform(clim.reshape(-1, 1)).ravel()
    X = sm.add_constant(z, has_constant="add")               # [pools x 2]

    # per-pool integer allele counts (phase-1: AF * flowers * 2, round). NaN AF -> count 0 (masked).
    aft = np.ascontiguousarray(af.T)                          # [records x pools]
    del af
    finite = np.isfinite(aft)
    succ = np.where(finite, np.round(aft * N[None, :]), 0.0)
    fail = np.where(finite, np.round((1.0 - aft) * N[None, :]), 0.0)
    succ[~finite] = 0.0; fail[~finite] = 0.0

    # effective-N: down-weight each pool by 1/(pools in its site) so the total
    # binomial information equals the number of independent climates (~31 sites),
    # not the 355 pools (pools within a site share one climate -> pseudoreplicated).
    # Within-site flower-weighting is retained via the per-pool N. See STATUS_clq90.
    w = None
    if args.eff_n == "site":
        site = pools["site"].to_numpy()
        ppsite = pd.Series(site).map(pd.Series(site).value_counts()).to_numpy(float)
        w = 1.0 / ppsite
        print(f"  effective-N=site: {len(np.unique(site))} sites, "
              f"pools/site {int(ppsite.min())}-{int(ppsite.max())}, "
              f"sum(var_weights)={w.sum():.1f} (~n_sites)", flush=True)

    print(f"{args.cls} gen{args.gen}: {len(clim)} pools x {n_rec:,} records vs "
          f"{args.climate} (range {np.nanmin(clim):.1f}-{np.nanmax(clim):.1f}); "
          f"N(genomes) {int(N.min())}-{int(N.max())}", flush=True)

    global _X, _SUCC, _FAIL, _W
    _X = X; _SUCC = succ; _FAIL = fail; _W = w

    step = 5000
    ranges = [(a, min(a + step, n_rec)) for a in range(0, n_rec, step)]
    slope = np.full(n_rec, np.nan); pv = np.full(n_rec, np.nan)
    with Pool(args.threads) as pool:
        for a, sl, p in pool.imap_unordered(_chunk, ranges):
            slope[a:a + len(sl)] = sl; pv[a:a + len(p)] = p

    recs = recs.assign(slope=slope, pval=pv)
    recs = recs.rename(columns={"maf": "MAF"})
    out = f"{args.out}/binomial_{args.cls}_gen{args.gen}_{args.climate}.csv"
    recs[["chrom", "pos", "ref_len", "alt_len", "MAF", "block",
          "slope", "pval"]].to_csv(out, index=False)
    fin = np.isfinite(pv)
    print(f"  tested {int(fin.sum()):,}/{n_rec:,} | p<0.05: {int((pv[fin]<0.05).sum()):,} "
          f"| min p={np.nanmin(pv) if fin.any() else float('nan'):.2e}\n  -> {out}", flush=True)


if __name__ == "__main__":
    main()
