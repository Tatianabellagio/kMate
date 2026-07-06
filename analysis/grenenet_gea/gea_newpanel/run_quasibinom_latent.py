#!/usr/bin/env python
"""Overdispersion-corrected binomial GEA (two routes), WITH K latent factors.

The plain per-variant binomial GLM `[alt,ref] ~ const + z(clim) + LF_1..LF_K`
treats N = AF*flowers*2 (up to ~680 genomes/pool) as independent Bernoulli draws.
The pool AF is really estimated from ~6x reads over tens of flowers, so that count
is ~10x too precise; the Wald variance is ~10x too small and GIF stays ~12 even
with K=16 factors (factors fix STRUCTURE, not the wrong count variance).

Two fixes, each layered ON TOP of the K latent factors:

(A) PER-VARIANT QUASI-BINOMIAL. Fit the binomial+LF GLM at full N, estimate the
    variant's dispersion  phi_i = PearsonChi2 / df_resid, clip phi_i>=1, and
    scale the climate Wald statistic  z_i -> z_i/sqrt(phi_i); recompute p from N(0,1).
    Per-variant phi (not one global scaling == genomic control, which over-corrects
    to 0 hits) keeps a real-signal tail while deflating the null bulk.

(B) EFFECTIVE SAMPLE SIZE (ACER-style pool-seq effective coverage). The counts
    overstate independent info; the pool AF is a two-stage sample (n_chrom = 2*flowers
    individuals, then ~mean_coverage reads), so the effective number of independent
    observations is
        n_eff = 1 / (1/n_chrom + 1/coverage)
    (harmonic combination; Kolaczkowski/ACER). Median here n_eff~5 vs N~56 => ~11x
    down-weight. Refit the binomial+LF GLM with N_eff counts; the likelihood itself
    now carries the right precision, so the Wald p is directly deflated.

Reuses latent_factors() (LFMM's U) and the multiprocessing-chunk GLM machinery of
run_binomial_latent.py.

Out (--out): quasibinom_lf{K}_{class}_gen{g}_{climate}.csv
  chrom,pos,ref_len,alt_len,MAF,block,slope,pval_binom,pval_quasi,phi,pval_effN

Usage (basic env): $PY run_quasibinom_latent.py --class snp --k 16 --threads 16
"""
from __future__ import annotations
import argparse, os, sys, warnings
import numpy as np
import pandas as pd
from multiprocessing import Pool
import statsmodels.api as sm
from scipy import stats
from sklearn.preprocessing import StandardScaler
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import lib

warnings.filterwarnings("ignore")

_X = None            # [pools x (2+K)] const + z(clim) + K latent factors
_SUCC = _FAIL = None      # full-N alt/ref counts   [rec x pools]
_SUCC_E = _FAIL_E = None  # effective-N alt/ref counts
_NP = None           # n design params (for df_resid)


def _chunk(rng):
    a, b = rng
    m = b - a
    slope = np.full(m, np.nan)
    pbin = np.full(m, np.nan)
    pqua = np.full(m, np.nan)
    phi = np.full(m, np.nan)
    peff = np.full(m, np.nan)
    Fam = sm.families.Binomial()
    for i in range(a, b):
        j = i - a
        succ = _SUCC[i]; fail = _FAIL[i]; n = succ + fail
        ok = n > 0
        if ok.sum() < _NP + 2:                 # need df_resid >= 2
            continue
        s = succ[ok]; f = fail[ok]
        if s.sum() == 0 or f.sum() == 0:
            continue
        X = _X[ok]
        if np.ptp(X[:, 1]) == 0:
            continue
        # ---- (A) full-N binomial + quasi ----
        try:
            res = sm.GLM(np.column_stack([s, f]), X, family=Fam).fit()
            b1 = res.params[1]; z = res.tvalues[1]
            slope[j] = b1; pbin[j] = res.pvalues[1]
            dfres = res.df_resid
            ph = res.pearson_chi2 / dfres if dfres > 0 else np.nan
            if not np.isfinite(ph) or ph < 1.0:
                ph = 1.0
            phi[j] = ph
            zq = z / np.sqrt(ph)
            pqua[j] = 2.0 * stats.norm.sf(abs(zq))
        except Exception:
            pass
        # ---- (B) effective-N binomial ----
        try:
            se = _SUCC_E[i][ok]; fe = _FAIL_E[i][ok]
            if se.sum() > 0 and fe.sum() > 0:
                res2 = sm.GLM(np.column_stack([se, fe]), X, family=Fam).fit()
                peff[j] = res2.pvalues[1]
        except Exception:
            pass
    return a, slope, pbin, pqua, phi, peff


def latent_factors(af: np.ndarray, k: int) -> np.ndarray:
    """Top-k per-pool structure loadings from the [pools x rec] AF matrix (LFMM's U)."""
    if k == 0:
        return np.zeros((af.shape[0], 0))
    Xc = af.astype(np.float64)
    Xc -= Xc.mean(axis=0, keepdims=True)
    G = Xc @ Xc.T
    w, V = np.linalg.eigh(G)
    U = V[:, ::-1][:, :k]
    return StandardScaler().fit_transform(U)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--class", dest="cls", required=True,
                    choices=["snp", "sv", "smallindel", "nonsnp"])
    ap.add_argument("--gen", type=int, default=9)
    ap.add_argument("--climate", default="bio1")
    ap.add_argument("--k", type=int, default=16)
    ap.add_argument("--cmdir", default=f"{lib.GEA}/phase1_replication/class_matrices")
    ap.add_argument("--out", default=f"{lib.GEA}/gea_newpanel/quasibinom")
    ap.add_argument("--threads", type=int, default=16)
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    recs = pd.read_csv(f"{args.cmdir}/{args.cls}_gen{args.gen}.records.csv")
    pools = pd.read_csv(f"{args.cmdir}/gen{args.gen}.pools.csv")
    clim = pools[args.climate].to_numpy(float)
    N = np.round(pools["total_flowers"].to_numpy(float) * 2).astype(np.float64)
    cov = pools["mean_coverage"].to_numpy(float)
    Neff = 1.0 / (1.0 / N + 1.0 / cov)          # ACER effective pool-seq sample size
    af = np.load(f"{args.cmdir}/{args.cls}_gen{args.gen}_af.npy")     # [pools x rec]
    assert af.shape[0] == len(clim) and af.shape[1] == len(recs)
    n_rec = af.shape[1]

    U = latent_factors(af, args.k)
    z = StandardScaler().fit_transform(clim.reshape(-1, 1))
    X = np.column_stack([np.ones(len(clim)), z.ravel(), U])

    aft = np.ascontiguousarray(af.T); del af
    finite = np.isfinite(aft)
    succ = np.where(finite, np.round(aft * N[None, :]), 0.0).astype(np.float32)
    fail = np.where(finite, np.round((1.0 - aft) * N[None, :]), 0.0).astype(np.float32)
    succ_e = np.where(finite, np.round(aft * Neff[None, :]), 0.0).astype(np.float32)
    fail_e = np.where(finite, np.round((1.0 - aft) * Neff[None, :]), 0.0).astype(np.float32)
    del aft, finite

    print(f"{args.cls} gen{args.gen}: {len(clim)} pools x {n_rec:,} recs vs {args.climate}, "
          f"K={args.k} ({X.shape[1]} cols)\n"
          f"  N med={np.median(N):.0f}  cov med={np.median(cov):.1f}  "
          f"Neff med={np.median(Neff):.1f}  downweight med={np.median(N/Neff):.1f}x",
          flush=True)

    global _X, _SUCC, _FAIL, _SUCC_E, _FAIL_E, _NP
    _X = X; _SUCC = succ; _FAIL = fail; _SUCC_E = succ_e; _FAIL_E = fail_e
    _NP = X.shape[1]

    step = 4000
    ranges = [(a, min(a + step, n_rec)) for a in range(0, n_rec, step)]
    slope = np.full(n_rec, np.nan); pbin = np.full(n_rec, np.nan)
    pqua = np.full(n_rec, np.nan); phi = np.full(n_rec, np.nan); peff = np.full(n_rec, np.nan)
    with Pool(args.threads) as pool:
        for a, sl, pb, pq, ph, pe in pool.imap_unordered(_chunk, ranges):
            L = len(sl)
            slope[a:a+L] = sl; pbin[a:a+L] = pb; pqua[a:a+L] = pq
            phi[a:a+L] = ph; peff[a:a+L] = pe

    recs = recs.rename(columns={"maf": "MAF"}).assign(
        slope=slope, pval_binom=pbin, pval_quasi=pqua, phi=phi, pval_effN=peff)
    out = f"{args.out}/quasibinom_lf{args.k}_{args.cls}_gen{args.gen}_{args.climate}.csv"
    recs[["chrom", "pos", "ref_len", "alt_len", "MAF", "block",
          "slope", "pval_binom", "pval_quasi", "phi", "pval_effN"]].to_csv(out, index=False)

    def gif(p):
        p = p[np.isfinite(p)]
        return np.median(stats.chi2.isf(p, 1)) / stats.chi2.ppf(0.5, 1)
    fin = np.isfinite(pbin)
    print(f"  tested {int(fin.sum()):,}/{n_rec:,} | phi med={np.nanmedian(phi):.2f}\n"
          f"  GIF binom={gif(pbin):.2f}  quasi={gif(pqua):.2f}  effN={gif(peff):.2f}\n"
          f"  minp binom={np.nanmin(pbin):.1e} quasi={np.nanmin(pqua):.1e} effN={np.nanmin(peff):.1e}\n"
          f"  -> {out}", flush=True)


if __name__ == "__main__":
    main()
