#!/usr/bin/env python
"""Latent-factor binomial GEA — the LFMM structure correction inside the binomial.

LFMM's engine is not the GIF, it is the K LATENT FACTORS: it fits
    Y = X*beta + U*V + E          (U*V = K factors absorbing confounding structure)
and the GIF only touches up residual inflation afterward. That U*V slot maps cleanly
onto the per-record binomial GLM: add K latent factors (PCs of the pool AF matrix,
estimated like LFMM's U from the class's own response) as fixed covariates:

    [alt, ref] ~ const + z(climate) + LF_1 + ... + LF_K

so the climate slope is tested *net of* the K structure axes. (Kendall has no
covariate slot, so this applies only to the parametric test; a *linear* model of
Delta p with K factors would simply reproduce LFMM.)

Latent factors: center each variant column of the [pools x rec] AF matrix, take the
top-K eigenvectors of the pool-pool gram matrix (== left singular vectors == the
per-pool structure loadings), standardized. K=16 matches the LFMM run.

Reads build_class_matrices.py outputs; writes (--out):
  binomial_lf{K}_{class}_gen{g}_{climate}.csv  chrom,pos,ref_len,alt_len,MAF,block,slope,pval

Usage (basic env):
  $PY run_binomial_latent.py --class snp --gen 9 --climate bio1 --k 16 --threads 16
"""
from __future__ import annotations
import argparse, os, sys, warnings
import numpy as np
import pandas as pd
from multiprocessing import Pool
import statsmodels.api as sm
from sklearn.preprocessing import StandardScaler
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import lib

warnings.filterwarnings("ignore")

_X = None       # [pools x (2+K)]: const + z(climate) + K latent factors
_SUCC = None    # [rec x pools] alt counts
_FAIL = None    # [rec x pools] ref counts


def _chunk(rng):
    a, b = rng
    slope = np.full(b - a, np.nan); pv = np.full(b - a, np.nan)
    for i in range(a, b):
        succ = _SUCC[i]; fail = _FAIL[i]; n = succ + fail
        ok = n > 0
        if ok.sum() < _X.shape[1] + 1:                # need > n_params pools
            continue
        s = succ[ok]; f = fail[ok]
        if s.sum() == 0 or f.sum() == 0:
            continue
        X = _X[ok]
        if np.ptp(X[:, 1]) == 0:
            continue
        try:
            res = sm.GLM(np.column_stack([s, f]), X, family=sm.families.Binomial()).fit()
            slope[i - a] = res.params[1]; pv[i - a] = res.pvalues[1]   # climate = col 1
        except Exception:
            pass
    return a, slope, pv


def latent_factors(af: np.ndarray, k: int) -> np.ndarray:
    """Top-k per-pool structure loadings from the [pools x rec] AF matrix (LFMM's U)."""
    Xc = af.astype(np.float64)
    Xc -= Xc.mean(axis=0, keepdims=True)              # center each variant (== center Delta p)
    G = Xc @ Xc.T                                      # [pools x pools] gram
    w, V = np.linalg.eigh(G)                           # ascending eigenvalues
    U = V[:, ::-1][:, :k]                              # top-k eigenvectors [pools x k]
    return StandardScaler().fit_transform(U)          # standardize factor columns


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--class", dest="cls", required=True,
                    choices=["snp", "sv", "smallindel", "nonsnp"])
    ap.add_argument("--gen", type=int, default=9)
    ap.add_argument("--climate", default="bio1")
    ap.add_argument("--k", type=int, default=16)
    ap.add_argument("--cmdir", default=f"{lib.GEA}/phase1_replication/class_matrices")
    ap.add_argument("--out", default=f"{lib.GEA}/gea_newpanel/binomial_latent")
    ap.add_argument("--threads", type=int, default=8)
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    recs = pd.read_csv(f"{args.cmdir}/{args.cls}_gen{args.gen}.records.csv")
    pools = pd.read_csv(f"{args.cmdir}/gen{args.gen}.pools.csv")
    clim = pools[args.climate].to_numpy(float)
    N = np.round(pools["total_flowers"].to_numpy(float) * 2).astype(np.float64)
    af = np.load(f"{args.cmdir}/{args.cls}_gen{args.gen}_af.npy")     # [pools x rec]
    assert af.shape[0] == len(clim) and af.shape[1] == len(recs)
    n_rec = af.shape[1]

    U = latent_factors(af, args.k)                                   # [pools x K]
    z = StandardScaler().fit_transform(clim.reshape(-1, 1))          # [pools x 1]
    X = np.column_stack([np.ones(len(clim)), z.ravel(), U])          # const, clim, factors

    aft = np.ascontiguousarray(af.T); del af
    finite = np.isfinite(aft)
    succ = np.where(finite, np.round(aft * N[None, :]), 0.0)
    fail = np.where(finite, np.round((1.0 - aft) * N[None, :]), 0.0)

    print(f"{args.cls} gen{args.gen}: {len(clim)} pools x {n_rec:,} records vs "
          f"{args.climate}, K={args.k} latent factors ({X.shape[1]} design cols)", flush=True)

    global _X, _SUCC, _FAIL
    _X = X; _SUCC = succ; _FAIL = fail

    step = 5000
    ranges = [(a, min(a + step, n_rec)) for a in range(0, n_rec, step)]
    slope = np.full(n_rec, np.nan); pv = np.full(n_rec, np.nan)
    with Pool(args.threads) as pool:
        for a, sl, p in pool.imap_unordered(_chunk, ranges):
            slope[a:a + len(sl)] = sl; pv[a:a + len(p)] = p

    recs = recs.rename(columns={"maf": "MAF"}).assign(slope=slope, pval=pv)
    out = f"{args.out}/binomial_lf{args.k}_{args.cls}_gen{args.gen}_{args.climate}.csv"
    recs[["chrom", "pos", "ref_len", "alt_len", "MAF", "block", "slope", "pval"]].to_csv(out, index=False)
    fin = np.isfinite(pv)
    print(f"  tested {int(fin.sum()):,}/{n_rec:,} | p<1e-5: {int((pv[fin]<1e-5).sum()):,} "
          f"| p<0.05: {int((pv[fin]<0.05).sum()):,} | min p={np.nanmin(pv):.2e}\n  -> {out}", flush=True)


if __name__ == "__main__":
    main()
