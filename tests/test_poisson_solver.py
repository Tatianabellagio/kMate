"""
Try Poisson NLL solver instead of WLS.

Poisson MLE objective: maximize Σ_k [c_k log μ_k − μ_k]
Equivalent to minimize Σ_k [μ_k - c_k log μ_k]
CVXPY: cp.sum(μ) - cp.sum(cp.multiply(c, cp.log(μ)))
This is convex in μ.

Test on uniform-82 simulation.
"""
from __future__ import annotations
import sys, os, time
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import numpy as np
import pandas as pd
import cvxpy as cp
from scipy.sparse import load_npz


DATA = os.path.join(os.path.dirname(__file__), "..", "data")


def solve_poisson(counts, cn, coverage, omega=None, eps=1e-3):
    K = counts.shape[0]
    F = cn.shape[0]
    if omega is None: omega = np.zeros(K)
    h = cp.Variable(F, nonneg=True)
    mu = coverage * (cn.T @ h) + coverage * omega + eps
    # Poisson NLL: minimize Σ [μ - c log μ]
    obj = cp.Minimize(cp.sum(mu) - cp.sum(cp.multiply(counts, cp.log(mu))))
    prob = cp.Problem(obj, [cp.sum(h) == 1])
    prob.solve(solver="SCS", verbose=False)
    return np.asarray(h.value)


def solve_l2_normalized(counts, cn, coverage, omega=None):
    """L2 on the relative-rate scale: c_k/(λ × AC_k) ≈ Σ_f h[f] cn[f,k]/AC_k"""
    K = counts.shape[0]
    F = cn.shape[0]
    if omega is None: omega = np.zeros(K)
    ac = cn.sum(axis=0)
    # design matrix: cn[f,k] / AC[k]  → predicts h-weighted "presence fraction"
    rate_obs = counts / np.maximum(coverage * ac, 1e-6)  # observed rate per k-mer
    cn_norm = cn.astype(float) / np.maximum(ac[None, :], 1)  # F × K, sums to 1 over founders
    h = cp.Variable(F, nonneg=True)
    pred = cn_norm.T @ h  # K-vector, sums to 1 only if h sums to 1 and cn_norm columns sum to 1
    residual = rate_obs - pred
    obj = cp.Minimize(cp.sum_squares(residual))
    prob = cp.Problem(obj, [cp.sum(h) == 1])
    prob.solve(solver="SCS", verbose=False)
    return np.asarray(h.value)


def solve_kl(counts, cn, coverage, omega=None, eps=1e-3):
    """KL-divergence loss (= Poisson NLL up to constants).
    cp.kl_div(c, μ) = c log(c/μ) - c + μ ≥ 0, convex.
    For c=0: kl_div = μ. For c>0: penalizes μ near 0 strongly.
    """
    K = counts.shape[0]
    F = cn.shape[0]
    if omega is None: omega = np.zeros(K)
    h = cp.Variable(F, nonneg=True)
    mu = coverage * (cn.T @ h) + coverage * omega + eps
    obj = cp.Minimize(cp.sum(cp.kl_div(counts, mu)))
    prob = cp.Problem(obj, [cp.sum(h) == 1])
    prob.solve(solver="SCS", verbose=False)
    return np.asarray(h.value)


def main():
    cn_sparse = load_npz(os.path.join(DATA, "test_chr1_first200.cn.npz"))
    meta = np.load(os.path.join(DATA, "test_chr1_first200.meta.npz"), allow_pickle=True)
    founders = meta["founders"]
    F, K = cn_sparse.shape
    cn = np.asarray(cn_sparse.todense()).astype(np.int8)
    ac = cn.sum(axis=0)
    counts = np.load(os.path.join(DATA, "sim_chr1", "uniform82_counts.npz"))["counts"]
    cov = counts.sum() * F / ac.sum()
    h_true = np.full(F, 1.0/F)

    print(f"cn: {F} × {K:,}, cov={cov:.1f}×")
    print(f"truth: uniform 1/{F} = {1/F:.4f}")

    def report(label, h_hat):
        h_hat = np.asarray(h_hat) / np.asarray(h_hat).sum()
        err = np.linalg.norm(h_hat - h_true)
        cv = h_hat.std() / h_hat.mean()
        nz = (h_hat < 1e-4).sum()
        print(f"  {label:55s}  ||Δh||={err:.4f}  CV={cv:.3f}  n_zero={nz}/{F}")

    print(f"\n{'-'*72}")
    # P1 Poisson NLL via kl_div
    t = time.time()
    h_p = solve_kl(counts, cn, cov)
    report(f"P1 Poisson NLL via kl_div  [{time.time()-t:.0f}s]", h_p)

    # P2 KL with proper coverage scan
    for c_try in [cov*0.5, cov, cov*2.0]:
        h = solve_kl(counts, cn, c_try)
        report(f"P2 kl_div, cov={c_try:.1f}×", h)

    # P3 L2 on normalized rates
    h_n = solve_l2_normalized(counts, cn, cov)
    report("P3 L2 on rate=c/(cov×AC)", h_n)

    # P4 Poisson on AC≥2 only
    keep = ac >= 2
    h_p2 = solve_kl(counts[keep], cn[:, keep], cov)
    report(f"P4 KL, AC≥2 only (n={keep.sum():,})", h_p2)

    # P5 Poisson on AC≥5 only
    keep5 = ac >= 5
    h_p5 = solve_kl(counts[keep5], cn[:, keep5], cov)
    report(f"P5 KL, AC≥5 only (n={keep5.sum():,})", h_p5)

    # show top founders by inferred h for the best variant
    print(f"\n  Top 10 founders by Poisson-h (P1):")
    for i in np.argsort(-h_p)[:10]:
        print(f"    {founders[i]:>10s}: h={h_p[i]/h_p.sum():.4f}")
    print(f"  Bottom 10 founders by Poisson-h (P1):")
    for i in np.argsort(h_p)[:10]:
        print(f"    {founders[i]:>10s}: h={h_p[i]/h_p.sum():.4f}")


if __name__ == "__main__":
    main()
