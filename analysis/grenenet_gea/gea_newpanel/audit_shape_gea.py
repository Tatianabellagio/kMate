#!/usr/bin/env python
"""Adversarial audit of shape_gea.py. Three concrete checks:
  (1) F-test math vs statsmodels on random loci (is the vectorized nested F right?)
  (2) WLS/flowers permutation-null calibration: our shortcut permutes rows of the
      *whitened* Y with a fixed hat matrix — valid for OLS, but in WLS the weight is
      tied to the SITE, so this scrambles the weight<->response pairing. Recompute the
      CORRECT WLS permutation (permute climate in X, rebuild H, weights fixed to site)
      on a subset and compare the genome-wide quad rate.
  (3) logit-clipping leverage: on 31 points a fixed/lost site (p9~0/1) becomes a
      logit outlier that can manufacture curvature. Quantify, and re-do the SV-vs-SNP
      curvature ratio excluding loci with any near-fixed site.
"""
from __future__ import annotations
import os, sys, time
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import lib
import shape_gea as S

rng = np.random.default_rng(1)


def build_XY(cls):
    p9, recs, p0 = S.load_class(cls)
    sites, sp, flw = S.site_frame()
    Y = S.lib_logit(np.clip(p9, S.EPS, 1 - S.EPS))
    X, lin, quad = S.build_design(sp, S.ENVS)
    return Y, X, lin, quad, flw, recs, p9, p0, sp


# ---------------------------------------------------------------- (1) F math
def check_f_math():
    import scipy.stats as st
    Y, X, lin, quad, flw, recs, p9, p0, sp = build_XY("snp")
    n = X.shape[0]; k1, k2 = 3, 5
    H0, H1, H2 = S.hat(X[:, :1]), S.hat(X[:, :k1]), S.hat(X)
    Fl, Fq, Ff = S.f_stats(Y, H0, H1, H2, 1, k1, k2, n)
    idx = rng.integers(0, Y.shape[1], 6)
    print("(1) F_quad vs manual nested-F (should match to ~1e-6):")
    for j in idx:
        y = Y[:, j]
        # manual: reduced [1,lin] vs full [1,lin,quad]
        def rss(Xm):
            b, *_ = np.linalg.lstsq(Xm, y, rcond=None); r = y - Xm @ b; return r @ r
        r1 = rss(X[:, :k1]); r2 = rss(X)
        Fman = ((r1 - r2) / (k2 - k1)) / (r2 / (n - k2))
        p_par = st.f.sf(Fman, k2 - k1, n - k2)
        print(f"   locus {j}: F_vec={Fq[j]:.6f}  F_manual={Fman:.6f}  diff={abs(Fq[j]-Fman):.2e}  parp={p_par:.3f}")


# ------------------------------------------------- (2) WLS null calibration
def correct_wls_quad_rate(cls, nsub=40000, nperm=500):
    """Genome-wide quad p<.05 rate under the CORRECT WLS permutation on a subset."""
    Y, X, lin, quad, flw, recs, p9, p0, sp = build_XY(cls)
    L = Y.shape[1]; sub = rng.choice(L, min(nsub, L), replace=False)
    Y = Y[:, sub]; n = X.shape[0]; k1, k2 = 3, 5
    w = flw / flw.mean(); sw = np.sqrt(w)
    Yw = Y * sw[:, None]
    Xbase = X.copy()                                   # [1 | z1 z12 | q1 q12]
    # observed (weighted)
    Xw = Xbase * sw[:, None]
    H0, H1, H2 = S.hat(Xw[:, :1]), S.hat(Xw[:, :k1]), S.hat(Xw)
    Fq_obs = S.f_stats(Yw, H0, H1, H2, 1, k1, k2, n)[1]
    # ---- our SHORTCUT null: permute rows of Yw, fixed H (what shape_gea does) ----
    cnt_short = np.zeros(len(sub), np.int32)
    for _ in range(nperm):
        pm = rng.permutation(n)
        Fq = S.f_stats(Yw[pm, :], H0, H1, H2, 1, k1, k2, n)[1]
        cnt_short += (Fq >= Fq_obs)
    rate_short = ((cnt_short + 1) / (nperm + 1) < 0.05).mean() * 100
    # ---- CORRECT null: permute climate rows in X (weights tied to site), rebuild H ----
    cnt_ok = np.zeros(len(sub), np.int32)
    clim = Xbase[:, 1:]                                # linear+quad climate cols
    for _ in range(nperm):
        pm = rng.permutation(n)
        Xp = np.column_stack([np.ones(n), clim[pm, :]]) * sw[:, None]
        H0p, H1p, H2p = S.hat(Xp[:, :1]), S.hat(Xp[:, :k1]), S.hat(Xp)
        Fq = S.f_stats(Yw, H0p, H1p, H2p, 1, k1, k2, n)[1]
        cnt_ok += (Fq >= Fq_obs)
    rate_ok = ((cnt_ok + 1) / (nperm + 1) < 0.05).mean() * 100
    print(f"(2) {cls} WLS quad p<.05 rate on {len(sub):,} loci: "
          f"shortcut(shape_gea)={rate_short:.2f}%  CORRECT-perm={rate_ok:.2f}%  (target ~5%)")


# ------------------------------------------------- (3) logit-clip leverage
def check_clip_leverage(tol=0.02):
    for cls in ["sv", "snp"]:
        p9, recs, p0 = S.load_class(cls)
        z = np.load(f"{S.OUT}/shape_{cls}_none.npz", allow_pickle=True)
        pq = z["p_quad"]
        extreme = ((p9 < tol) | (p9 > 1 - tol)).any(0)   # any near-fixed/lost site
        q = pq < 0.05
        print(f"(3) {cls}: {extreme.mean()*100:.1f}% loci have >=1 near-fixed site "
              f"| quad p<.05 among extreme={q[extreme].mean()*100:.2f}% vs clean={q[~extreme].mean()*100:.2f}%")
    # SV-vs-matched-SNP quad ratio, clean loci only
    p9s, rs, p0s = S.load_class("snp"); zs = np.load(f"{S.OUT}/shape_snp_none.npz")["p_quad"]
    p9v, rv, p0v = S.load_class("sv"); zv = np.load(f"{S.OUT}/shape_sv_none.npz")["p_quad"]
    ext_s = ((p9s < tol) | (p9s > 1 - tol)).any(0); ext_v = ((p9v < tol) | (p9v > 1 - tol)).any(0)
    order = np.argsort(p0s); ps = p0s[order]
    mi = order[np.clip(np.searchsorted(ps, p0v), 0, len(ps) - 1)]
    for lab, mask_v, mask_s in [("ALL", np.ones(len(zv), bool), np.ones(len(zs), bool))]:
        pass
    def ratio(keep_v):
        rv_ = (zv[keep_v] < 0.05).mean(); rs_ = (zs[mi[keep_v]] < 0.05).mean()
        return rv_ * 100, rs_ * 100, rv_ / max(rs_, 1e-9)
    a = ratio(np.ones(len(zv), bool)); b = ratio(~ext_v)
    print(f"(3) SV/matchedSNP quad ratio  ALL: SV={a[0]:.2f}% SNP={a[1]:.2f}% ratio={a[2]:.3f}"
          f"  | CLEAN(no near-fixed site): SV={b[0]:.2f}% SNP={b[1]:.2f}% ratio={b[2]:.3f}")


if __name__ == "__main__":
    t = time.time()
    check_f_math()
    for c in ["snp", "sv"]:
        correct_wls_quad_rate(c)
    check_clip_leverage()
    print(f"[audit done {time.time()-t:.0f}s]")
