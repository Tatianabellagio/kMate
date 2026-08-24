#!/usr/bin/env python
"""Sign-flip permutation null + FWER threshold for the sampling-variance-weighted block LD-LMM.

The parametric residual test (z_b ~ N(0,1)) is calibrated in the bulk after the sampling-variance
floor, but a mild correlated-residual / EM-projection inflation remains in the upper-middle tail
(site 4: lambda_0.75-0.9 ~ 1.1-1.18) that NO parametric null can capture. This builds the EMPIRICAL
null by SIGN-FLIPPING the per-plot temporal slopes:

  for a sign vector eps in {+1,-1}^n (one entry per replicate plot, the SAME vector applied to all
  blocks):  slopes_perm[j,:] = eps_j * slopes[j,:]   ->   s_b^perm = mean_j eps_j slopes[j,b]
  then re-run the EXACT Stage-2 pipeline (recompute replicate-SE -> D0, BLUP linked background,
  residual z) at the FITTED variance components (tau, sigma_e^2).

This preserves: (a) boundary curvature (the frequency values, hence per-block slope magnitudes, are
reused); (b) cross-block residual correlation (one eps for all blocks); (c) founder/LD structure
(A and C_LD unchanged; the random founder loadings that survive each flip are re-absorbed by the
re-fit BLUP). It DESTROYS the across-plot-consistent block-specific signal -> a valid H0 for
"selection beyond linkage". With n=11 plots we ENUMERATE all 2^n = 2048 sign patterns exactly
(|z| is invariant to a global sign flip, so 1024 are distinct; we keep all for simplicity).

Threshold: maxT FWER -- per pattern record max_b |z_b^perm|; the 95th percentile is the genome-wide
5% threshold. Per-block FWER p_b = (1 + #{patterns: maxabsz >= |z_b^obs|}) / (1 + n_pattern).
Also a pooled-null sample for a marginal calibrated tail. Writes site<ID>_block_ld_lmm_permnull.csv
(observed + p_perm/p_fwer) + _meta.json. Env: kmate. SITE via env.
"""
import os, sys, json, itertools
import numpy as np
import pandas as pd
from scipy import stats, optimize
from scipy.linalg import cho_factor, cho_solve

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lib
from founder_genotype import build_genotype

H = "results/grenenet_gea/hapfreq"
SITE = int(os.environ.get("SITE", 4))
EPS = 1e-3
Tg = np.array([0.0, 1.0, 2.0, 3.0])


def logit(p):
    p = np.clip(p, EPS, 1 - EPS)
    return np.log(p / (1 - p))


def logit_var(f, Ngam):
    f = np.clip(f, EPS, 1 - EPS)
    return 1.0 / (np.clip(Ngam, 1.0, None) * f * (1 - f))


def main():
    G, founders, reg = build_genotype()
    mat = np.load(f"{H}/hapfreq_matrix.npy")
    samples = [l.strip() for l in open(f"{H}/hapfreq_samples.txt") if l.strip()]
    sidx = {s: i for i, s in enumerate(samples)}
    p0 = np.load(f"{H}/hapfreq_p0_seedmix.npy")
    fr = pd.read_csv(f"{H}/hapfreq_registry.csv")
    nhap = mat.shape[1]

    # ---- STAGE 1 (identical to block_ld_lmm_temporal_sampvar.py) ----
    pt = lib.pool_table()
    s = pt[(pt.site == SITE) & pt.sampleid.astype(str).isin(sidx)].copy()
    cell, cellN = {}, {}
    for (gen, plot), g in s.groupby(["generation", "plot"]):
        rows = g.sampleid.astype(str).map(sidx).to_numpy()
        w = g.flowerscollected.to_numpy(float); w = np.where(np.isfinite(w) & (w > 0), w, 1.0)
        cell[(int(gen), int(plot))] = (mat[rows] * w[:, None]).sum(0) / w.sum()
        cellN[(int(gen), int(plot))] = float(w.sum())
    present = {}
    for (gen, plot) in cell:
        present.setdefault(plot, set()).add(gen)
    plots = sorted(pl for pl, gs in present.items() if {1, 2, 3} <= gs)
    n = len(plots)
    coef = (Tg - 1.5)
    Nmed = np.median([cellN[(g, pl)] for (g, pl) in cellN])
    v0 = logit_var(p0, 2 * Nmed)
    v0_term = (coef[0] / 5.0) ** 2 * v0

    slopes = np.zeros((n, nhap))
    evolved_samp = np.zeros(nhap)
    for j, pl in enumerate(plots):
        cells = [cell[(1, pl)], cell[(2, pl)], cell[(3, pl)]]
        Ns = [cellN[(1, pl)], cellN[(2, pl)], cellN[(3, pl)]]
        y = np.vstack([logit(p0), logit(cells[0]), logit(cells[1]), logit(cells[2])])
        slopes[j] = (coef[:, None] * y).sum(0) / 5.0
        for ti, (c, Nc) in enumerate(zip(cells, Ns), start=1):
            evolved_samp += (coef[ti] / 5.0) ** 2 * logit_var(c, 2 * Nc)
    evolved_samp /= n ** 2

    fr["unit"] = fr.chrom + ":" + fr.unit_start.astype(str) + "-" + fr.unit_end.astype(str)
    testable = (fr.covered & fr.panel_freq.between(0.05, 0.95)).to_numpy()
    keep = np.zeros(nhap, bool)
    for u, grp in fr[testable].groupby("unit"):
        kept = grp.index.drop(grp.panel_freq.idxmax()) if len(grp) > 1 else grp.index
        keep[kept.to_numpy()] = True
    ki = np.where(keep)[0]; M = len(ki)

    SL = slopes[:, ki]                                  # (n, M) per-plot slopes, testable blocks
    v0k = v0_term[ki]; sampk = evolved_samp[ki]
    A = G[:, ki].astype(np.float64)
    pf = A.mean(0); A = (A - pf) / np.sqrt(pf * (1 - pf) + 1e-9)
    nF = A.shape[0]
    one = np.ones(M)

    def D0_of(SLp):
        rep = SLp.std(0, ddof=1) ** 2 / n                # replicate var of the mean slope
        D0 = v0k + np.maximum(rep, sampk)
        return np.clip(D0, np.median(v0k + np.maximum(rep, sampk)) * 1e-3, None)

    def pieces(tau, se2, D0):
        Winv = 1.0 / (se2 + D0)
        AWA = (A * Winv[None, :]) @ A.T / nF
        cf = cho_factor(np.eye(nF) / tau + AWA, lower=True)

        def Vinv(x):
            Wx = Winv * x
            u = cho_solve(cf, (A @ Wx) / np.sqrt(nF))
            return Wx - Winv * ((A.T @ u) / np.sqrt(nF))
        logdetV = np.sum(np.log(se2 + D0)) + np.linalg.slogdet(np.eye(nF) + tau * AWA)[1]
        return Vinv, logdetV

    def resid_z(SLp, tau, se2):
        """residual z for a (possibly sign-flipped) per-plot slope matrix, at fixed (tau,se2)."""
        D0 = D0_of(SLp)
        sb = SLp.mean(0); sb = sb - sb.mean()
        Vinv, _ = pieces(tau, se2, D0)
        Vi1 = Vinv(one); mu = (one @ Vinv(sb)) / (one @ Vi1)
        r = sb - mu
        g_b = (tau / nF) * (A.T @ (A @ Vinv(r)))
        resid = r - g_b
        return resid / np.sqrt(se2 + D0)

    # ---- observed fit: REML for (tau, se2) on the real data ----
    D0 = D0_of(SL)
    sb0 = SL.mean(0); sb0 = sb0 - sb0.mean()

    def neg_reml(par):
        tau, se2 = np.exp(par)
        Vinv, logdetV = pieces(tau, se2, D0)
        Vi1 = Vinv(one); mu = (one @ Vinv(sb0)) / (one @ Vi1)
        r = sb0 - mu
        return 0.5 * (logdetV + r @ Vinv(r) + np.log(one @ Vi1))

    res = optimize.minimize(neg_reml, x0=np.log([np.var(sb0), np.var(sb0) / 2]),
                            method="Nelder-Mead", options={"xatol": 1e-3, "fatol": 1e-3})
    tau, se2 = np.exp(res.x)
    z_obs = resid_z(SL, tau, se2)
    absz_obs = np.abs(z_obs)

    # ---- sign-flip permutation null: enumerate all 2^n patterns ----
    patterns = list(itertools.product([1.0, -1.0], repeat=n))
    np_ = len(patterns)
    maxabs = np.empty(np_)
    pool = []                                            # subsample of null |z| for marginal tail
    ge_obs = np.zeros(M, dtype=np.int64)                 # per-block: #patterns with maxabsz >= |z_obs_b|
    rng = np.random.default_rng(0)
    for k, eps in enumerate(patterns):
        eps = np.array(eps)
        if np.all(eps == eps[0]):                        # global flip == observed (|z| invariant)
            zk = absz_obs
        else:
            zk = np.abs(resid_z(eps[:, None] * SL, tau, se2))
        mk = float(zk.max())
        maxabs[k] = mk
        ge_obs += (mk >= absz_obs)
        pool.append(zk[rng.integers(0, M, size=64)])     # 64 random null z per pattern
        if (k + 1) % 256 == 0:
            print(f"  ...{k+1}/{np_} patterns")
    pool = np.concatenate(pool)

    thr05 = float(np.quantile(maxabs, 0.95))             # FWER 5% maxT threshold
    thr10 = float(np.quantile(maxabs, 0.90))
    p_fwer = (1 + ge_obs) / (1 + np_)                    # per-block FWER p-value (maxT)
    # marginal calibrated p from the pooled null (two-sided already, |z|)
    pool_sorted = np.sort(pool)
    p_marg = 1.0 - np.searchsorted(pool_sorted, absz_obs, side="right") / len(pool_sorted)
    p_marg = np.clip(p_marg, 1.0 / len(pool_sorted), 1.0)

    d = fr.loc[ki, ["chrom", "unit_start", "unit_end", "panel_freq"]].copy()
    d["s"] = sb0; d["z"] = z_obs; d["p_param"] = 2 * stats.norm.sf(absz_obs)
    d["p_marg_perm"] = p_marg; d["p_fwer"] = p_fwer
    m = len(d); o = d.p_marg_perm.values.argsort()
    q = np.empty(m); q[o] = np.minimum.accumulate(
        (d.p_marg_perm.values[o] * m / (np.arange(m) + 1))[::-1])[::-1]
    d["q_perm"] = np.clip(q, 0, 1)
    d.sort_values("p_marg_perm").to_csv(f"{H}/site{SITE}_block_ld_lmm_permnull.csv", index=False)

    n_fwer = int((d.p_fwer < 0.05).sum())
    n_qperm = int((d.q_perm < 0.05).sum())
    lam_param = np.median(stats.chi2.isf(np.clip(d.p_param, 1e-300, 1), 1)) / stats.chi2.ppf(0.5, 1)
    meta = dict(site=SITE, n_plot=n, M=int(M), n_patterns=np_, tau=float(tau), se2=float(se2),
                maxT_thr_fwer05=thr05, maxT_thr_fwer10=thr10,
                obs_max_absz=float(absz_obs.max()), lambda_param_median=float(lam_param),
                n_fwer05=n_fwer, n_qperm05=n_qperm,
                pool_z_q=[float(np.quantile(pool, qq)) for qq in [0.5, 0.95, 0.99, 0.999]])
    json.dump(meta, open(f"{H}/site{SITE}_block_ld_lmm_permnull_meta.json", "w"), indent=2)

    print(f"\nsite {SITE}: {M:,} blocks, {n} plots, {np_} sign-flip patterns (exact)")
    print(f"  fitted tau={tau:.4f}, se2={se2:.4f}")
    print(f"  observed max |z| = {absz_obs.max():.2f}")
    print(f"  maxT FWER threshold: 5% |z| > {thr05:.2f}   10% |z| > {thr10:.2f}")
    print(f"  pooled null |z| quantiles 50/95/99/99.9 = "
          + "/".join(f"{np.quantile(pool, qq):.2f}" for qq in [0.5, 0.95, 0.99, 0.999])
          + "   (param N(0,1): 0.67/1.96/2.58/3.29)")
    print(f"  hits: FWER p<0.05 = {n_fwer} | permutation-BH q<0.05 = {n_qperm}")
    print("\n  TOP blocks by permutation marginal p:")
    cols = ["chrom", "unit_start", "unit_end", "panel_freq", "s", "z", "p_param", "p_marg_perm", "p_fwer", "q_perm"]
    print(d.sort_values("p_marg_perm").head(12)[cols].to_string(index=False))
    print(f"\n[done] {H}/site{SITE}_block_ld_lmm_permnull.csv + _meta.json")


if __name__ == "__main__":
    main()
