#!/usr/bin/env python
"""Pipeline B v2 — pooled-trajectory selection coefficients + IV-weighted climate GEA.

Implements analysis/grenenet_gea/PIPELINE_B_POOLED_MODEL.md exactly:
  Step 1  pool persistent plots -> flower-weighted site frequency p̄_{t,g} (gen0 = p0).
  Step 2  per-point log-odds variance V_{t,g} = max(drift σ², binomial v_samp)/n_g
          (variance-components, on logit scale; no double-count). Drift var shrunk
          toward the per-SV pool (d0=4). Gen0 var from the 8 SEEDMIX reps.
  Step 3  weighted GLS slope of logit(p̄) on t -> s_g and SE(s_g) (one fit, ω=1/V).
  Step 4  random-effects meta-regression s_g ~ β0 + β1·climate, w=1/(SE²+τ̂²) (DL);
          β1 calibrated by site-label permutation (two-sided + one-sided up-in-warm).
          β1 also reported under 1/SE², n_g, equal weights (robustness).

Output: twostage_scoefpooled_bio1.npz (same schema as build_two_stage_gea -> WZA &
significant-blocks work unchanged) + printed diagnostics.

  PY=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python
  $PY analysis/grenenet_gea/build_two_stage_pooled.py --climate bio1
"""
from __future__ import annotations
import argparse, glob, os, sys
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lib

PM = f"{lib.GEA}/pool_matrices"; STORE = lib.AF_STORE
NC_MIN, SV_MIN_BP, EPS, D0, N_PERM = 150, 50, 1e-3, 4.0, 2000
PLOIDY = 2


def logit(p):
    p = np.clip(p, EPS, 1 - EPS)
    return np.log(p / (1 - p))


def sv_mask_and_meta():
    idx = np.load(f"{STORE}/index_nonsnp.npz")
    rl = idx["ref_len"].astype(np.int64); al = idx["alt_len"].astype(np.int64)
    size = np.abs(al - rl)
    nc = np.asarray(np.load(sorted(glob.glob(f"{STORE}/nc_nonsnp/*.npy"))[0]))
    mask = (size > SV_MIN_BP) & (nc >= NC_MIN)
    meta = pd.DataFrame(dict(chrom=idx["chrom"].astype("U5")[mask], pos=idx["pos"][mask],
                             ref_len=rl[mask], alt_len=al[mask], sv_size=size[mask]))
    return mask, meta


def gen0_logit_var(mask):
    """Among-SEEDMIX-rep variance of logit(p0), /n_reps. Handles old-panel SEEDMIX TSVs
    (10.33M) via the old2new mask. Falls back to a small floor if unavailable."""
    snp = np.load(f"{STORE}/snp_mask.npy")
    full_sv = np.zeros(len(snp), bool); full_sv[~snp] = mask        # SV rows in 8.49M panel
    N_OLD = 10_325_364
    o2n = np.load(f"{STORE}/old2new_mask.npy") if os.path.exists(f"{STORE}/old2new_mask.npy") else None
    fs = [f for f in glob.glob(f"{lib.SEEDMIX}/SEEDMIX_S*.tsv") if "_Chr" not in os.path.basename(f)]
    L = []
    for f in fs:
        a = pd.read_csv(f, sep="\t", usecols=["alt_freq"]).alt_freq.to_numpy()
        if len(a) == N_OLD and o2n is not None:
            a = a[o2n]                                             # 10.33M -> 8.49M subset
        if len(a) != len(full_sv):
            print(f"  [gen0] {os.path.basename(f)} len {len(a):,} != panel — skipping"); continue
        L.append(logit(a[full_sv]))
    if len(L) < 2:
        return None
    Lm = np.vstack(L)                                              # [n_reps x nSV]
    print(f"  [gen0] using {len(L)} SEEDMIX reps for founding variance", flush=True)
    return np.nanvar(Lm, axis=0, ddof=1) / Lm.shape[0]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--climate", default="bio1")
    ap.add_argument("--n-perm", type=int, default=N_PERM)
    ap.add_argument("--out", default=f"{lib.GEA}/gea")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)
    clim = args.climate

    mask, sv_meta = sv_mask_and_meta()
    sv_idx = np.where(mask)[0]; nSV = len(sv_idx)
    p0 = np.load(f"{STORE}/p0_nonsnp.npy")[mask].astype(np.float64)
    V0 = gen0_logit_var(mask)
    if V0 is None:
        V0 = np.full(nSV, 1e-4)
    print(f"SVs: {nSV:,} | gen0 logit-var median {np.nanmedian(V0):.2e}", flush=True)

    # ---- load gens, find persistent plots (present in gen1,2,3) ----
    P, M, key = {}, {}, {}
    for g in (1, 2, 3):
        P[g] = np.load(f"{PM}/pool_gen{g}_nonsnp_af.npy")[:, sv_idx].astype(np.float64)
        M[g] = pd.read_csv(f"{PM}/pool_gen{g}_nonsnp.meta.csv")
        key[g] = (M[g].site.astype(str) + "_" + M[g]["plot"].astype(str)).to_numpy()
    persistent = set(key[1]) & set(key[2]) & set(key[3])
    sites = np.array(sorted(M[3].site.unique()))
    bio1 = M[3].groupby("site")[clim].first().loc[sites].to_numpy(float)
    print(f"persistent plots: {len(persistent)} | sites: {len(sites)}", flush=True)

    # ---- Step 1+2: per (site, gen) pooled freq, drift var, sampling var ----
    nS = len(sites)
    ybar = {1: {}, 2: {}, 3: {}}; sig2 = {1: {}, 2: {}, 3: {}}
    vsamp = {1: {}, 2: {}, 3: {}}; ng = {1: {}, 2: {}, 3: {}}
    for g in (1, 2, 3):
        row_of = {k: i for i, k in enumerate(key[g])}
        flw = M[g]["total_flowers"].to_numpy(float)
        cov = M[g]["mean_coverage"].to_numpy(float)
        for si, s in enumerate(sites):
            plots = [k for k in np.unique(key[g][M[g].site.to_numpy() == s]) if k in persistent]
            rows = [row_of[k] for k in plots]
            if not rows:
                ng[g][si] = 0; continue
            psub = P[g][rows]                                      # [n x SV]
            f = flw[rows]; c = cov[rows]; n = len(rows)
            pbar = (f @ psub) / f.sum()                           # flower-weighted mean [SV]
            ybar[g][si] = logit(pbar)
            pq = np.clip(pbar * (1 - pbar), 1e-6, None)
            neff = 1.0 / (1.0 / (PLOIDY * f) + 1.0 / c)            # per-plot effective N
            vsamp[g][si] = (np.mean(1.0 / neff)) / pq             # v̄_samp on logit scale [SV]
            sig2[g][si] = (np.var(logit(psub), axis=0, ddof=1) if n >= 2
                           else np.full(nSV, np.nan))             # among-plot logit var
            ng[g][si] = n

    # per-SV pooled drift variance (df-weighted across cells with n>=2)
    num = np.zeros(nSV); den = 0.0
    for g in (1, 2, 3):
        for si in range(nS):
            if ng[g][si] >= 2:
                num += (ng[g][si] - 1) * np.nan_to_num(sig2[g][si]); den += (ng[g][si] - 1)
    sig2_pool = num / max(den, 1.0)

    # ---- Step 3: weighted GLS slope per site (vectorized over SVs) ----
    VFLOOR = 1e-6                                                  # avoid inf weights (V=0)
    V0 = np.maximum(V0, VFLOOR)
    s_g = np.full((nS, nSV), np.nan); SE_g = np.full((nS, nSV), np.nan)
    site_ok = np.zeros(nS, bool)
    tvec = np.array([0., 1., 2., 3.])
    for si in range(nS):
        if not all(si in ybar[g] for g in (1, 2, 3)):
            continue
        site_ok[si] = True
        Y = [logit(p0)] + [ybar[g][si] for g in (1, 2, 3)]        # 4 x [SV]
        # per-point variance V_{t,g}
        Vs = [V0]
        for g in (1, 2, 3):
            df = ng[g][si] - 1
            s2 = sig2[g][si]
            s2_shr = (np.nan_to_num(s2) * df + D0 * sig2_pool) / (df + D0) if df >= 1 \
                else sig2_pool                                    # n=1 -> pool
            Vs.append(np.maximum.reduce([s2_shr, vsamp[g][si], np.full(nSV, VFLOOR)]) / ng[g][si])
        Y = np.vstack(Y); Vw = np.maximum(np.vstack(Vs), VFLOOR)  # [4 x SV]
        w = 1.0 / Vw
        Sw = w.sum(0); Swt = (w * tvec[:, None]).sum(0)
        Swtt = (w * tvec[:, None] ** 2).sum(0)
        Swy = (w * Y).sum(0); Swty = (w * tvec[:, None] * Y).sum(0)
        Delta = Sw * Swtt - Swt ** 2
        s_g[si] = (Sw * Swty - Swt * Swy) / Delta
        SE_g[si] = np.sqrt(Sw / Delta)

    # ---- Step 4: random-effects meta-regression + permutation ----
    use = np.where(site_ok)[0]                                    # sites with persistent plots (site-level)
    sg = s_g[use]; seg = SE_g[use]; x = bio1[use]
    x = (x - x.mean()) / x.std()
    k = len(use)
    print(f"sites used in Stage 2: {k}", flush=True)

    def wls_beta(w, xx):
        Sw = w.sum(0); Swx = (w * xx[:, None]).sum(0)
        Swxx = (w * xx[:, None] ** 2).sum(0)
        Sws = (w * sg).sum(0); Swxs = (w * xx[:, None] * sg).sum(0)
        D = Sw * Swxx - Swx ** 2
        b1 = (Sw * Swxs - Swx * Sws) / D
        b0 = (Sws - b1 * Swx) / Sw
        return b0, b1, D

    # DL tau^2 from 1/SE^2-weighted residuals
    w0 = 1.0 / (seg ** 2)
    b0, b1, _ = wls_beta(w0, x)
    resid = sg - (b0[None, :] + b1[None, :] * x[:, None])
    Q = (w0 * resid ** 2).sum(0)
    C = w0.sum(0) - (w0 ** 2).sum(0) / w0.sum(0)
    tau2 = np.maximum(0.0, (Q - (k - 2)) / np.where(C > 0, C, np.nan))
    w = 1.0 / (seg ** 2 + tau2[None, :])
    _, beta, _ = wls_beta(w, x)

    # weight-concentration diagnostic: Kish effective # sites = (Σw)²/Σw² per SV
    eff_n = (w.sum(0) ** 2) / (w ** 2).sum(0)
    share = (w / w.sum(0)[None, :]).mean(1)            # mean weight share per site
    print(f"WEIGHT CONCENTRATION: effective #sites (Kish) median {np.nanmedian(eff_n):.1f} "
          f"of {k} | p10 {np.nanpercentile(eff_n,10):.1f} p90 {np.nanpercentile(eff_n,90):.1f}")
    print(f"  top-1/top-3 site mean weight share: {np.sort(share)[::-1][0]:.2f} / "
          f"{np.sort(share)[::-1][:3].sum():.2f}", flush=True)

    # robustness: other weightings
    ng3 = np.array([ng[3][si] for si in use], float)
    _, beta_se2, _ = wls_beta(w0, x)
    _, beta_ng, _ = wls_beta(np.repeat(ng3[:, None], nSV, 1), x)
    _, beta_eq, _ = wls_beta(np.ones((k, nSV)), x)

    # site-label permutation null (per-SV weights w)
    rng = np.random.default_rng(0)
    beta_perm = np.empty((args.n_perm, nSV))
    for b in range(args.n_perm):
        xp = x[rng.permutation(k)]
        _, bp, _ = wls_beta(w, xp)
        beta_perm[b] = bp
    sd = beta_perm.std(0); sd[sd == 0] = np.nan
    z_emp = beta / sd
    znull = (beta_perm / sd[None, :]).ravel(); znull = znull[np.isfinite(znull)]
    nN = len(znull); sn = np.sort(znull); asn = np.sort(np.abs(znull))
    p_perm = np.clip(1.0 - np.searchsorted(asn, np.abs(z_emp)) / nN, 1.0 / nN, 1.0)
    p_perm_up = np.clip((nN - np.searchsorted(sn, z_emp)) / nN, 1.0 / nN, 1.0)
    p_perm[~np.isfinite(z_emp)] = np.nan; p_perm_up[~np.isfinite(z_emp)] = np.nan

    from scipy.stats import chi2, norm
    g = np.isfinite(p_perm) & (p_perm > 0)
    gif = float(np.median(norm.isf(p_perm[g] / 2) ** 2) / chi2.ppf(0.5, 1))

    out = f"{args.out}/twostage_scoefpooled_{clim}.npz"
    np.savez(out, beta=beta.astype(np.float32), z_emp=z_emp.astype(np.float32),
             p_perm=p_perm.astype(np.float32), p_perm_up=p_perm_up.astype(np.float32),
             tau2=tau2.astype(np.float32),
             chrom=sv_meta.chrom.to_numpy().astype("U5"), pos=sv_meta.pos.to_numpy(),
             ref_len=sv_meta.ref_len.to_numpy(), alt_len=sv_meta.alt_len.to_numpy(),
             sv_size=sv_meta.sv_size.to_numpy(), p0=p0.astype(np.float32),
             n_sites=k, gif_emp=gif)

    # diagnostics
    fin = np.isfinite(beta) & np.isfinite(beta_se2)
    print(f"\nGIF (perm) = {gif:.2f} | p_perm<0.05: {(p_perm<0.05).sum():,} "
          f"(exp {int(0.05*np.isfinite(p_perm).sum()):,}) | min {np.nanmin(p_perm):.1e}")
    print(f"tau2 median {np.nanmedian(tau2):.2e} | frac tau2>0 {np.mean(tau2>0):.2f}")
    print("β1 weighting robustness (corr to chosen 1/(SE²+τ²)):")
    for nm, bb in [("1/SE²", beta_se2), ("n_g", beta_ng), ("equal", beta_eq)]:
        print(f"   {nm:7s}: r={np.corrcoef(beta[fin], bb[fin])[0,1]:.3f}")
    print(f"-> {out}")


if __name__ == "__main__":
    main()
