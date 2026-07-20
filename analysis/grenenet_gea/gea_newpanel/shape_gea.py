#!/usr/bin/env python
"""Response-function ("shape") GEA — beyond the linear LFMM.

Standard LFMM tests ONE thing per locus: a monotone-linear allele-frequency /
climate relationship (Y = beta*X + latent). In an experimental-evolution frame
(Y = per-site selection response) that captures only **antagonistic pleiotropy /
clinal** adaptation (up in warm gardens, down in cold). It is blind to:

  - **intermediate optimum** ("adaptive only in the middle") -> a concave *hump*
    in Δp vs climate; the linear term integrates to ~0, LFMM has ~no power.
  - **conditional neutrality** (adaptive one end, neutral the other) -> a
    one-sided *ramp/hinge*; linear has only partial power and mis-specifies shape.

This module builds "our own LFMM" for arbitrary response shapes by BASIS-
EXPANDING the environment and doing nested tests, at the honest **site level**
(31 gardens = 31 independent climate values), calibrated by a **climate-
permutation null** (shuffle the 31 climate labels; preserves genome covariance
and whatever inflation the parametric F carries -> no GIF fudge, unlike the
K=16 pool-level LFMM whose GIF never ->1).

Per locus, on the variance-stabilised **logit frequency** scale (kills the
bounded-Δp fake-curvature artifact; p0 is a single shared founder frequency so
it drops into the intercept), fit for envs {bio1, bio12}:

    logit(p9_site)  ~  1  +  [lin: z1, z12]  +  [quad: q1, q12]

with orthogonal within-env quads. Sequential F tests:

    F_lin  : intercept        vs +linear     (2 df) -> clinal / antag. pleiotropy
    F_quad : +linear          vs +quad        (2 df) -> curvature / interm. optimum
    F_full : intercept        vs +lin+quad   (4 df) -> any climate shape

Two response weightings (the user asked to compare):
    --weight none    : OLS on logit(p9)             (raw logit-Δp shape)
    --weight flowers : WLS, w_site ∝ site flower census (drift-aware ~ s-coef)

Classes: snp | indel | sv  (indel = nonsnp |Δlen|<=50, sv = >=50; split from the
existing nonsnp site matrix — no rebuild). Output: per-class npz of per-locus
F/empirical-p/quad-signs for the SNP-vs-SV shape comparison notebook.
"""
from __future__ import annotations
import argparse, os, sys, time
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import lib  # noqa: E402

CMDIR = f"{lib.GEA}/phase1_replication/results/class_matrices"
SITEDIR = f"{lib.GEA}/gea_newpanel/results/lfmm_site"
OUT = f"{lib.GEA}/gea_newpanel/results/shape_gea"
P0 = {"snp": "p0_snp.npy", "nonsnp": "p0_nonsnp.npy"}
EPS = 1e-4
ENVS = ["bio1", "bio12"]


# ----------------------------------------------------------------------------- data
def load_class(cls):
    """Return (p9 [31 x L], records DataFrame, p0 [L]) for snp|indel|sv|nonsnp."""
    base = "snp" if cls == "snp" else "nonsnp"
    recs = pd.read_csv(f"{CMDIR}/{base}_gen9.records.csv")
    n_site, n_rec = map(int, open(f"{SITEDIR}/lfmm_{base}_site_dims.txt").read().split())
    assert n_rec == len(recs), f"{base}: dims {n_rec} != records {len(recs)}"
    dp = np.fromfile(f"{SITEDIR}/lfmm_{base}_site_Y.f64", dtype=np.float64).reshape(n_site, n_rec)
    p0 = np.load(f"{lib.AF_STORE}/{P0[base]}")[recs["col"].to_numpy()].astype(np.float64)
    p9 = dp + p0[None, :]                                   # site_af = Δp + shared p0
    if cls in ("indel", "sv"):                              # size-split the non-SNP set
        sz = (recs.alt_len - recs.ref_len).abs().to_numpy()
        keep = (sz <= 50) if cls == "indel" else (sz >= 51)
        recs = recs.loc[keep].reset_index(drop=True)
        p9 = p9[:, keep]; p0 = p0[keep]
    return p9, recs, p0


def site_frame():
    """31 sites (sorted, matching build_site_matrices) with env + flower weight."""
    pools = pd.read_csv(f"{CMDIR}/gen9.pools.csv")
    sites = np.sort(pools.site.unique())
    sp = pools.drop_duplicates("site").set_index("site").loc[sites]
    flw = pools.groupby("site").total_flowers.sum().loc[sites].to_numpy(float)
    return sites, sp, flw


# ----------------------------------------------------------------------- design
def build_design(sp, envs):
    """[31 x k] columns: intercept, linear(env)..., quad(env)... (within-env orth).
    Returns X, and the slices for the linear and quad blocks."""
    n = len(sp)
    lins, quads = [], []
    for e in envs:
        x = sp[e].to_numpy(float)
        z = (x - x.mean()) / x.std(ddof=0)
        q = z * z
        # residualise q against [1, z] then standardise -> pure curvature, sign(+z^2) kept
        A = np.c_[np.ones(n), z]
        q = q - A @ np.linalg.lstsq(A, q, rcond=None)[0]
        q = q / q.std(ddof=0)
        lins.append(z); quads.append(q)
    # blocks contiguous: [intercept | linear... | quad...]
    X = np.column_stack([np.ones(n)] + lins + quads)
    lin = list(range(1, 1 + len(envs)))
    quad = list(range(1 + len(envs), 1 + 2 * len(envs)))
    return X, lin, quad


def hat(X):
    """Projection (hat) matrix H = X (X'X)^-1 X'  [n x n]."""
    return X @ np.linalg.solve(X.T @ X, X.T)


# ------------------------------------------------------------------- statistics
def rss_via_H(Y, H):
    """Residual SS per column for design with hat matrix H. Y [n x L]."""
    return (Y * Y).sum(0) - (Y * (H @ Y)).sum(0)


def f_stats(Y, H0, H1, H2, k0, k1, k2, n):
    """Sequential F per locus: F_lin (H0->H1), F_quad (H1->H2), F_full (H0->H2)."""
    r0 = rss_via_H(Y, H0); r1 = rss_via_H(Y, H1); r2 = rss_via_H(Y, H2)
    r1 = np.maximum(r1, 1e-12); r2 = np.maximum(r2, 1e-12)
    F_lin = ((r0 - r1) / (k1 - k0)) / (r1 / (n - k1))
    F_quad = ((r1 - r2) / (k2 - k1)) / (r2 / (n - k2))
    F_full = ((r0 - r2) / (k2 - k0)) / (r2 / (n - k2))
    return F_lin, F_quad, F_full


def quad_sign(Y, X, quad_idx):
    """Full-model coefficient sign on each quad column (<0 = concave = hump)."""
    B = np.linalg.solve(X.T @ X, X.T @ Y)              # [k x L]
    return B[quad_idx, :]                               # [n_env x L]


# ------------------------------------------------------------------------- main
def run(cls, weight, n_perm, seed, smoke):
    t0 = time.time()
    p9, recs, p0 = load_class(cls)
    sites, sp, flw = site_frame()
    n = len(sites)
    Y = lib_logit(np.clip(p9, EPS, 1 - EPS))            # [31 x L] response
    if smoke:
        Y = Y[:, :smoke]; recs = recs.iloc[:smoke].reset_index(drop=True); p0 = p0[:smoke]
    L = Y.shape[1]

    X, lin, quad = build_design(sp, ENVS)
    k0, k1, k2 = 1, 1 + len(lin), 1 + len(lin) + len(quad)
    if weight == "flowers":                             # WLS: whiten rows by sqrt(w)
        w = flw / flw.mean()
        sw = np.sqrt(w)[:, None]
        Yw = Y * sw; Xw = X * np.sqrt(w)[:, None]
    else:
        Yw, Xw = Y, X

    H0 = hat(Xw[:, :k0]); H1 = hat(Xw[:, :k1]); H2 = hat(Xw)
    Fl, Fq, Ff = f_stats(Yw, H0, H1, H2, k0, k1, k2, n)
    qsign = quad_sign(Yw, Xw, quad)                     # [n_env x L]

    # ---- climate-permutation null: shuffle the 31 rows of Y (design fixed) ----
    # RSS_perm(design) with fixed Y  ==  RSS(design) with row-permuted Y  (SST invariant)
    rng = np.random.default_rng(seed)
    cnt_l = np.zeros(L, np.int32); cnt_q = np.zeros(L, np.int32); cnt_f = np.zeros(L, np.int32)
    for _ in range(n_perm):
        perm = rng.permutation(n)
        Yp = Yw[perm, :]
        fl, fq, ff = f_stats(Yp, H0, H1, H2, k0, k1, k2, n)
        cnt_l += (fl >= Fl); cnt_q += (fq >= Fq); cnt_f += (ff >= Ff)
    p_lin = (cnt_l + 1) / (n_perm + 1)
    p_quad = (cnt_q + 1) / (n_perm + 1)
    p_full = (cnt_f + 1) / (n_perm + 1)

    sz = (recs.alt_len - recs.ref_len).abs().to_numpy()
    out = dict(
        chrom=recs.chrom.to_numpy().astype("U5"), pos=recs.pos.to_numpy(np.int64),
        sv_size=sz.astype(np.int32), maf=recs.maf.to_numpy(np.float32), p0=p0.astype(np.float32),
        F_lin=Fl.astype(np.float32), F_quad=Fq.astype(np.float32), F_full=Ff.astype(np.float32),
        p_lin=p_lin.astype(np.float32), p_quad=p_quad.astype(np.float32), p_full=p_full.astype(np.float32),
        qsign_bio1=qsign[0].astype(np.float32), qsign_bio12=qsign[1].astype(np.float32),
    )
    os.makedirs(OUT, exist_ok=True)
    fn = f"{OUT}/shape_{cls}_{weight}.npz"
    np.savez(fn, **out)
    # quick genome-wide summary
    sig_q = (p_quad < 0.05); sig_l = (p_lin < 0.05)
    hump = sig_q & (qsign[0] < 0)
    print(f"[{cls}/{weight}] L={L:,} | lin p<.05 {sig_l.mean()*100:.2f}% | "
          f"quad p<.05 {sig_q.mean()*100:.2f}% ({sig_q.sum():,}) | of those concave(bio1) "
          f"{ (hump.sum()/max(sig_q.sum(),1))*100:.0f}% | quad-only(no lin) "
          f"{((sig_q&~sig_l).sum()):,} | wrote {os.path.basename(fn)} | {time.time()-t0:.0f}s",
          flush=True)


def lib_logit(p):
    return np.log(p) - np.log1p(-p)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--cls", required=True, choices=["snp", "indel", "sv", "nonsnp"])
    ap.add_argument("--weight", default="none", choices=["none", "flowers"])
    ap.add_argument("--n-perm", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--smoke", type=int, default=0, help="only first N loci (debug)")
    a = ap.parse_args()
    run(a.cls, a.weight, a.n_perm, a.seed, a.smoke)
