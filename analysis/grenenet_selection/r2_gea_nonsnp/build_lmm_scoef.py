#!/usr/bin/env python
"""Per-SV random-intercept LMM on the selection coefficients (pipeline B, Stage 2').

Replaces the two-stage average+meta-regression with the plain mixed model the user
asked for — NO manual weighting, just a site random intercept to absorb the
non-independence of replicate plots within a site:

    s_ij  ~  climate_site  +  (1 | site)

where s_ij = per-plot-lineage logit-slope selection coefficient (build_two_stage_gea
stage1_scoef), fit by REML (statsmodels MixedLM) for each SV across its 142
plot-lineages / 20 sites. Output: bio1 fixed-effect beta, p, SE, convergence, GIF.

Caveat: climate is a site-level (Level-2) predictor, so the (1|site) intercept
competes with it; with ~20 sites the parametric p can be conservative/degenerate
(this is what collapsed the earlier pool-level model). GIF reported so we can see.

  PY=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python
  $PY analysis/grenenet_selection/r2_gea_nonsnp/build_lmm_scoef.py --climate bio1 --threads 16
"""
from __future__ import annotations
import argparse, os, sys, warnings
import numpy as np
import pandas as pd
from multiprocessing import Pool
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import lib
from build_two_stage_gea import sv_mask_and_meta, stage1_scoef
warnings.filterwarnings("ignore")

_Y = None; _X = None; _GRP = None     # shared via fork: s [nSV x nL], climate, site


def _fit_chunk(rng):
    import statsmodels.formula.api as smf
    a, b = rng
    n = b - a
    beta = np.full(n, np.nan); pval = np.full(n, np.nan)
    se = np.full(n, np.nan); conv = np.zeros(n, bool)
    for i in range(a, b):
        y = _Y[i]
        if not np.all(np.isfinite(y)) or np.ptp(y) == 0:
            continue
        df = pd.DataFrame({"y": y, "climate": _X, "site": _GRP})
        try:
            m = smf.mixedlm("y ~ climate", df, groups=df["site"]).fit(
                reml=True, method="lbfgs")
            beta[i - a] = m.params["climate"]; pval[i - a] = m.pvalues["climate"]
            se[i - a] = m.bse["climate"]; conv[i - a] = m.converged
        except Exception:
            pass
    return a, beta, pval, se, conv


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--climate", default="bio1")
    ap.add_argument("--threads", type=int, default=16)
    ap.add_argument("--limit", type=int, default=0, help="debug: first N SVs only")
    ap.add_argument("--out", default=f"{lib.GEA}/gea")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    mask, sv_meta = sv_mask_and_meta()
    sv_idx = np.where(mask)[0]
    p0 = np.load(f"{lib.AF_STORE}/p0_nonsnp.npy")[mask].astype(np.float64)
    s, site, _, unit = stage1_scoef(sv_idx, p0)      # s [nL x nSV], site [nL]
    # climate per lineage (its site's bio1), standardized
    m3 = pd.read_csv(f"{lib.GEA}/pool_matrices/pool_gen3_nonsnp.meta.csv")
    bio = m3.groupby("site")[args.climate].first()
    x = bio.loc[site].to_numpy(float); x = (x - x.mean()) / x.std()

    global _Y, _X, _GRP
    _Y = np.ascontiguousarray(s.T)                   # [nSV x nL]
    _X = x; _GRP = np.asarray(site)
    if args.limit:
        _Y = _Y[:args.limit]; sv_meta = sv_meta.iloc[:args.limit]
    nSV = _Y.shape[0]
    print(f"LMM s~{args.climate}+(1|site): {nSV:,} SVs x {len(site)} lineages, "
          f"{len(np.unique(site))} sites, threads={args.threads}", flush=True)

    step = max(50, nSV // (args.threads * 8))
    ranges = [(a, min(a + step, nSV)) for a in range(0, nSV, step)]
    beta = np.empty(nSV); pval = np.empty(nSV); se = np.empty(nSV); conv = np.empty(nSV, bool)
    with Pool(args.threads) as pool:
        for a, bt, pv, s_, cv in pool.imap_unordered(_fit_chunk, ranges):
            k = len(bt); beta[a:a+k] = bt; pval[a:a+k] = pv; se[a:a+k] = s_; conv[a:a+k] = cv

    from scipy.stats import chi2, norm
    ok = np.isfinite(pval) & (pval > 0)
    gif = float(np.median(norm.isf(pval[ok] / 2) ** 2) / chi2.ppf(0.5, 1))
    out = f"{args.out}/lmm_scoef_{args.climate}.npz"
    np.savez(out, beta=beta, pval=pval, se=se, converged=conv,
             chrom=sv_meta.chrom.to_numpy().astype("U5"), pos=sv_meta.pos.to_numpy(),
             ref_len=sv_meta.ref_len.to_numpy(), alt_len=sv_meta.alt_len.to_numpy(),
             sv_size=sv_meta.sv_size.to_numpy(), p0=p0[:nSV], gif=gif)
    print(f"converged {int(conv.sum()):,}/{nSV:,} | finite p {int(ok.sum()):,}")
    print(f"GIF = {gif:.3f} | p<1e-4: {int((pval[ok]<1e-4).sum())} | min p {np.nanmin(pval):.2e}")
    print(f"-> {out}")


if __name__ == "__main__":
    main()
