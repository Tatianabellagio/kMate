#!/usr/bin/env python
"""Per-SV linear mixed model: log(p3/p0) ~ bio1 + (1|site), gen-3 SV climate GEA.

The principled fix for the pseudoreplication that inflated the naive Kendall and
forced LFMM's blunt GIF deflation: model the within-site correlation explicitly
with a site random intercept (partial pooling). If it works, the p-values are
self-calibrated (GIF ~ 1) WITHOUT genomic control, recovering signal LFMM's
uniform deflation removes.

For each of the 23,344 filtered SVs (MAF>=0.05 + min-count), across the 193 gen-3
pools:
  y = log(p3_clip / p0_clip)   (freqs clipped to [EPS,1-EPS] so lost/fixed alleles
                                stay finite),  x = standardized bio1,  groups = site.
  statsmodels MixedLM(y ~ bio1, groups=site, random intercept) -> bio1 coef + p.

Output (--out): mixedmodel_gen3_bio1.npz (beta, pval, se, converged + locus meta)
and the GIF of the resulting p-values (should be ~1, vs LFMM's ~2.2).
"""
from __future__ import annotations
import argparse, os, sys, warnings
import numpy as np
import pandas as pd
from multiprocessing import Pool
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import lib
warnings.filterwarnings("ignore")

EPS = 1e-3
_Y = None; _X = None; _GRP = None      # shared via fork


def _fit_chunk(rng):
    import statsmodels.formula.api as smf
    a, b = rng
    beta = np.full(b - a, np.nan); pval = np.full(b - a, np.nan)
    se = np.full(b - a, np.nan); conv = np.zeros(b - a, bool)
    for i in range(a, b):
        y = _Y[i]
        if not np.all(np.isfinite(y)) or np.ptp(y) == 0:
            continue
        df = pd.DataFrame({"y": y, "bio1": _X, "site": _GRP})
        try:
            m = smf.mixedlm("y ~ bio1", df, groups=df["site"]).fit(reml=True, method="lbfgs")
            beta[i - a] = m.params["bio1"]; pval[i - a] = m.pvalues["bio1"]
            se[i - a] = m.bse["bio1"]; conv[i - a] = m.converged
        except Exception:
            pass
    return a, beta, pval, se, conv


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gen", type=int, default=3)
    ap.add_argument("--climate", default="bio1")
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--limit", type=int, default=0, help="benchmark: only first N SVs")
    ap.add_argument("--out", default=f"{lib.GEA}/gea")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    D = f"{lib.GEA}/pool_matrices"
    mt = pd.read_csv(f"{D}/pool_gen{args.gen}_nonsnp.meta.csv")
    af = np.load(f"{D}/pool_gen{args.gen}_nonsnp_af.npy", mmap_mode="r")
    p0 = np.load(f"{lib.AF_STORE}/p0_nonsnp.npy")
    loci = pd.read_csv(f"{lib.GEA}/lfmm/locus_index_gen{args.gen}_sv.csv")
    rec = loci.rec_index.to_numpy()
    if args.limit:
        rec = rec[:args.limit]; loci = loci.iloc[:args.limit]

    p3 = np.asarray(af[:, rec], dtype=np.float64)               # [pools x nSV]
    p3c = np.clip(p3, EPS, 1 - EPS); p0c = np.clip(p0[rec], EPS, 1 - EPS)
    global _Y, _X, _GRP
    _Y = np.log(p3c / p0c[None, :]).T                           # [nSV x pools]
    x = mt[args.climate].to_numpy(float); _X = (x - x.mean()) / x.std()
    _GRP = mt.site.to_numpy()
    n = _Y.shape[0]
    print(f"gen{args.gen} {args.climate}: {n:,} SVs x {len(mt)} pools, "
          f"{mt.site.nunique()} sites, threads={args.threads}", flush=True)

    step = max(50, n // (args.threads * 8))
    ranges = [(a, min(a + step, n)) for a in range(0, n, step)]
    beta = np.empty(n); pval = np.empty(n); se = np.empty(n); conv = np.empty(n, bool)
    with Pool(args.threads) as pool:
        for a, bt, pv, s, cv in pool.imap_unordered(_fit_chunk, ranges):
            beta[a:a+len(bt)] = bt; pval[a:a+len(pv)] = pv
            se[a:a+len(s)] = s; conv[a:a+len(cv)] = cv

    ok = np.isfinite(pval)
    # GIF of the mixed-model p-values (chi-sq 1 df)
    from scipy.stats import chi2, norm
    z2 = norm.isf(pval[ok] / 2) ** 2
    gif = np.median(z2) / chi2.ppf(0.5, 1)
    np.savez(f"{args.out}/mixedmodel_gen{args.gen}_{args.climate}.npz",
             beta=beta, pval=pval, se=se, converged=conv,
             chrom=loci.chrom.to_numpy().astype("U5"), pos=loci.pos.to_numpy(),
             sv_size=loci.sv_size.to_numpy(), p0=loci.p0.to_numpy(),
             rec_index=rec, gif=gif)
    print(f"converged: {conv.sum():,}/{n:,} | finite p: {ok.sum():,}")
    print(f"GIF (mixed model) = {gif:.3f}   (LFMM was ~2.2; target ~1)")
    print(f"p<1e-4: {(pval[ok]<1e-4).sum()} | min p: {np.nanmin(pval):.2e}")
    print(f"-> {args.out}/mixedmodel_gen{args.gen}_{args.climate}.npz")


if __name__ == "__main__":
    main()
