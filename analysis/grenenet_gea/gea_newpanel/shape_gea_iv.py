#!/usr/bin/env python
"""Two-stage, replicate-aware response-function GEA (fixes the audit).

The first version collapsed the 7-12 plot replicates per site into ONE flower-
weighted mean and then ran OLS (or a broken WLS) across 31 equally-trusted means
-> it used only the sqrt(N) averaging benefit of the replicates, and it had two
bugs (mis-calibrated WLS permutation; logit-clip leverage at near-fixed sites).

This version uses the replicates for their real power:

  STAGE 1 (per site s, per locus): from the n_s=7..12 plot AFs, flower-weighted
    site mean Δp_s AND its variance V_s = among-plot SEM^2 + binomial floor
    (the 2-component point variance: among-plot drift + finite-sample sampling).
  STAGE 2 (per locus, across the 31 sites): inverse-variance-weighted (GLS/WLS)
    regression of Δp_s on the climate-shape basis, so a site whose 12 plots agree
    counts more than one with 7 noisy plots. Near-fixed / high-drift sites are
    down-weighted by their SE (no logit needed -> no clip artifact). Weights are
    capped at 10x the per-locus median so no single site dominates.

Response = raw Δp (site_af - shared p0); intercept absorbs p0. Bases per env
{bio1,bio12}: linear, orthogonal quadratic (symmetric hump = intermediate optimum),
and hinge max(z,0) at the climate mean (asymmetric ramp = conditional neutrality).
Sequential F: F_lin (intercept->+lin), F_quad (+lin->+quad), F_hinge (+lin->+hinge).

Null = CORRECT climate permutation: shuffle the 31 climate rows of the design,
keep each site's (Δp, weight) fixed, rebuild the weighted fit. weight in {ols, iv}
(ols = unit weights = clean unweighted baseline; iv = replicate precision).

Classes snp|indel|sv (indel = nonsnp |Δlen|<=50, sv >=51). Vectorised per-locus
WLS via batched 5x5 normal equations.
"""
from __future__ import annotations
import argparse, os, sys, time
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import lib
import shape_gea as S   # reuse load_class / site_frame / build_design / ENVS

OUT = f"{lib.GEA}/gea_newpanel/shape_gea"
CMDIR = S.CMDIR
WCAP = 10.0            # cap a site's weight at WCAP x per-locus median


# --------------------------------------------------- stage 1: per-site Δp + measurement var
def stage1(cls):
    """Return dp [31 x L] (=site_af-p0) and m [31 x L] (per-site MEASUREMENT variance
    of the site mean = among-plot SEM^2 + binomial floor). NOTE: this is only the
    within-site (plot-level) uncertainty; the site~climate regression also carries a
    shared site-level drift component tau^2, added later via DerSimonian-Laird."""
    base = "snp" if cls == "snp" else "nonsnp"
    recs = pd.read_csv(f"{CMDIR}/{base}_gen9.records.csv")
    af = np.load(f"{CMDIR}/{base}_gen9_af.npy")                 # [355 x L] float32
    p0 = np.load(f"{lib.AF_STORE}/{S.P0[base]}")[recs["col"].to_numpy()].astype(np.float64)
    pools = pd.read_csv(f"{CMDIR}/gen9.pools.csv")
    sites = np.sort(pools.site.unique())
    site_of = pools.site.to_numpy(); flw = pools.total_flowers.to_numpy(float)
    L = af.shape[1]
    site_af = np.empty((len(sites), L)); mvar = np.empty((len(sites), L))
    for i, s in enumerate(sites):
        m = site_of == s
        w = flw[m]; sw = w.sum(); w2 = (w * w).sum()
        P = af[m].astype(np.float64)                           # [n_s x L]
        mu = (w[:, None] * P).sum(0) / sw
        # unbiased weighted among-plot variance of a single plot, then SEM^2 of the mean
        denom = sw - w2 / sw
        v_among = (w[:, None] * (P - mu) ** 2).sum(0) / max(denom, 1e-9)
        sem2 = v_among * (w2 / sw ** 2)
        floor = (mu * (1 - mu) + 1e-3) / (2 * sw)              # binomial sampling floor
        site_af[i] = mu; mvar[i] = sem2 + floor
    if cls in ("indel", "sv"):
        sz = (recs.alt_len - recs.ref_len).abs().to_numpy()
        keep = (sz <= 50) if cls == "indel" else (sz >= 51)
        recs = recs.loc[keep].reset_index(drop=True)
        site_af = site_af[:, keep]; mvar = mvar[:, keep]; p0 = p0[keep]
    dp = site_af - p0[None, :]
    return dp, mvar, recs, p0, site_af


def iv_weights(dp, mvar):
    """Random-effects (DerSimonian-Laird) weights: 1/(mvar + tau^2), tau^2 from the
    climate-INDEPENDENT intercept-only heterogeneity so the permutation null stays
    exactly valid. Capped at WCAP x per-locus median."""
    w0 = 1.0 / mvar                                           # [n x L] within-site precision
    n = dp.shape[0]
    sw0 = w0.sum(0); mu = (w0 * dp).sum(0) / sw0              # weighted grand mean per locus
    Q = (w0 * (dp - mu) ** 2).sum(0)                          # heterogeneity statistic
    c = sw0 - (w0 * w0).sum(0) / sw0
    tau2 = np.maximum(0.0, (Q - (n - 1)) / np.maximum(c, 1e-12))   # DL among-site drift var
    V = 1.0 / (mvar + tau2[None, :])
    med = np.median(V, axis=0)
    return np.minimum(V, WCAP * med[None, :]), tau2


# ------------------------------------------------------------------ design (lin/quad/hinge)
def design(sp, envs):
    n = len(sp); lins, quads, hinges = [], [], []
    for e in envs:
        x = sp[e].to_numpy(float); z = (x - x.mean()) / x.std(ddof=0)
        q = z * z; A = np.c_[np.ones(n), z]
        q = q - A @ np.linalg.lstsq(A, q, rcond=None)[0]; q = q / q.std(ddof=0)
        h = np.maximum(z, 0.0)                                 # hinge at climate mean
        lins.append(z); quads.append(q); hinges.append(h)
    X = np.column_stack([np.ones(n)] + lins + quads + hinges)  # [1|lin|quad|hinge]
    ne = len(envs)
    idx = dict(intr=[0], lin=list(range(1, 1 + ne)),
               quad=list(range(1 + ne, 1 + 2 * ne)), hinge=list(range(1 + 2 * ne, 1 + 3 * ne)))
    return X, idx


# ----------------------------------------------------- vectorised per-locus weighted RSS
def wls_rss(Xc, Y, V):
    """Weighted RSS per locus for design columns Xc [n x k]; Y,V [n x L]."""
    XtWX = np.einsum('sp,sq,si->ipq', Xc, Xc, V, optimize=True)   # [L,k,k]
    XtWy = np.einsum('sp,si,si->ip', Xc, V, Y, optimize=True)     # [L,k]
    beta = np.linalg.solve(XtWX, XtWy[:, :, None])[:, :, 0]       # [L,k]
    resid = Y - np.einsum('sp,ip->si', Xc, beta, optimize=True)   # [n,L]
    return (V * resid * resid).sum(0)                            # [L]


def rss0(Y, V):
    mu = (V * Y).sum(0) / V.sum(0)
    return (V * (Y - mu) ** 2).sum(0)


def fset(X, idx, Y, V, n):
    """F_lin, F_quad, F_hinge per locus."""
    cL = idx['intr'] + idx['lin']
    cQ = cL + idx['quad']; cH = cL + idx['hinge']
    kL, kQ, kH = len(cL), len(cQ), len(cH)
    r0 = rss0(Y, V); r1 = wls_rss(X[:, cL], Y, V)
    r2 = wls_rss(X[:, cQ], Y, V); r3 = wls_rss(X[:, cH], Y, V)
    r1 = np.maximum(r1, 1e-12); r2 = np.maximum(r2, 1e-12); r3 = np.maximum(r3, 1e-12)
    F_lin = ((r0 - r1) / (kL - 1)) / (r1 / (n - kL))
    F_quad = ((r1 - r2) / (kQ - kL)) / (r2 / (n - kQ))
    F_hinge = ((r1 - r3) / (kH - kL)) / (r3 / (n - kH))
    return F_lin, F_quad, F_hinge


def quad_sign(X, idx, Y, V):
    cQ = idx['intr'] + idx['lin'] + idx['quad']; Xc = X[:, cQ]
    XtWX = np.einsum('sp,sq,si->ipq', Xc, Xc, V, optimize=True)
    XtWy = np.einsum('sp,si,si->ip', Xc, V, Y, optimize=True)
    beta = np.linalg.solve(XtWX, XtWy[:, :, None])[:, :, 0]       # [L,k]
    return beta[:, len(idx['intr']) + len(idx['lin'])]            # coef on first quad col (bio1)


# --------------------------------------------------------------------------------- run
def run(cls, weight, n_perm, seed, smoke):
    t0 = time.time()
    dp, mvar, recs, p0, site_af = stage1(cls)
    _, sp, _ = S.site_frame()
    n = dp.shape[0]
    if smoke:
        dp = dp[:, :smoke]; mvar = mvar[:, :smoke]; recs = recs.iloc[:smoke].reset_index(drop=True)
        p0 = p0[:smoke]
    if weight == "ols":
        V = np.ones_like(dp)
    else:
        V, tau2 = iv_weights(dp, mvar)
    Y = dp; L = Y.shape[1]
    X, idx = design(sp, S.ENVS)
    clim = X[:, 1:]                                              # non-intercept cols

    Fl, Fq, Fh = fset(X, idx, Y, V, n)
    qs = quad_sign(X, idx, Y, V)

    rng = np.random.default_rng(seed)
    cl = np.zeros(L, np.int32); cq = np.zeros(L, np.int32); ch = np.zeros(L, np.int32)
    for _ in range(n_perm):
        pm = rng.permutation(n)                                  # permute climate rows
        Xp = np.column_stack([np.ones(n), clim[pm, :]])
        fl, fq, fh = fset(Xp, idx, Y, V, n)                     # Y,V fixed to sites
        cl += (fl >= Fl); cq += (fq >= Fq); ch += (fh >= Fh)
    p_lin = (cl + 1) / (n_perm + 1); p_quad = (cq + 1) / (n_perm + 1); p_hinge = (ch + 1) / (n_perm + 1)

    sz = (recs.alt_len - recs.ref_len).abs().to_numpy()
    os.makedirs(OUT, exist_ok=True)
    fn = f"{OUT}/shapeiv_{cls}_{weight}.npz"
    np.savez(fn,
        chrom=recs.chrom.to_numpy().astype("U5"), pos=recs.pos.to_numpy(np.int64),
        sv_size=sz.astype(np.int32), maf=recs.maf.to_numpy(np.float32), p0=p0.astype(np.float32),
        F_lin=Fl.astype(np.float32), F_quad=Fq.astype(np.float32), F_hinge=Fh.astype(np.float32),
        p_lin=p_lin.astype(np.float32), p_quad=p_quad.astype(np.float32), p_hinge=p_hinge.astype(np.float32),
        qsign_bio1=qs.astype(np.float32))
    print(f"[{cls}/{weight}] L={L:,} | lin p<.05 {(p_lin<.05).mean()*100:.2f}% | "
          f"quad {(p_quad<.05).mean()*100:.2f}% | hinge {(p_hinge<.05).mean()*100:.2f}% | "
          f"concave-of-quad {((p_quad<.05)&(qs<0)).sum()/max((p_quad<.05).sum(),1)*100:.0f}% | "
          f"{time.time()-t0:.0f}s -> {os.path.basename(fn)}", flush=True)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--cls", required=True, choices=["snp", "indel", "sv", "nonsnp"])
    ap.add_argument("--weight", default="iv", choices=["iv", "ols"])
    ap.add_argument("--n-perm", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--smoke", type=int, default=0)
    a = ap.parse_args()
    run(a.cls, a.weight, a.n_perm, a.seed, a.smoke)
