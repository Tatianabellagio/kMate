"""
Per-block, per-sample CVXPY solver for pool-seq founder frequencies.

The model (per block b, sample i):
    E[c_{ik}] = λ_i · α_k · h_b · cn_k  +  ω_k · λ_i
    c_{ik} ~ Poisson(E[c_{ik}])

with the simplex constraint h_b ∈ Δ^F (h ≥ 0, sum h = 1).

Two solvers:
- solve_block_wls: one-shot weighted least squares with sensible default weights
- solve_block_irls: iteratively reweighted (Newton on Poisson NLL)

Both are DCP-compliant: μ_k(h) is affine in h, residual squared is convex.

(Earlier attempt: Anscombe transform 2·sqrt(c+3/8) − 2·sqrt(μ+3/8) which is
convex in h but its square is not DCP. Switched to weighted L2 on raw scale.)

See MODEL_SPEC.md for full math.
"""
from __future__ import annotations
import numpy as np
import cvxpy as cp


def solve_block_wls(
    counts: np.ndarray,          # K-vector of observed k-mer counts in this block
    cn: np.ndarray,              # F × K matrix: copy number of each k-mer in each founder (0/1)
    coverage: float,             # λ_i, expected k-mer coverage
    alpha: np.ndarray | float = 1.0,   # K-vector (or scalar) of per-k-mer reach
    omega: np.ndarray | None = None,   # K-vector of per-k-mer contamination rate (or None = 0)
    weights: np.ndarray | None = None, # K-vector of weights (or None = 1/sqrt(c+1))
    solver: str = "SCS",
    verbose: bool = False,
) -> tuple[np.ndarray, float]:
    """One-shot weighted least squares.

    Default weights w_k = 1 / sqrt(c_k + 1). This approximates the Poisson
    weighting w = 1/sqrt(μ) using observed counts as a proxy. Good for low-noise.
    """
    K = counts.shape[0]
    F = cn.shape[0]
    assert cn.shape[1] == K, f"cn shape {cn.shape} mismatch counts shape {counts.shape}"

    if omega is None:
        omega = np.zeros(K)
    if np.isscalar(alpha):
        alpha = np.full(K, alpha)
    if weights is None:
        weights = 1.0 / np.sqrt(counts + 1.0)

    h = cp.Variable(F, nonneg=True)
    mu = coverage * cp.multiply(alpha, cn.T @ h) + coverage * omega    # affine in h
    residual = cp.multiply(weights, counts - mu)                        # affine in h
    obj = cp.Minimize(cp.sum_squares(residual))                         # convex
    cons = [cp.sum(h) == 1]
    prob = cp.Problem(obj, cons)
    prob.solve(solver=solver, verbose=verbose)
    if h.value is None:
        raise RuntimeError(f"CVXPY solver failed: status={prob.status}")
    return h.value, prob.value


def solve_block_irls(
    counts: np.ndarray,
    cn: np.ndarray,
    coverage: float,
    alpha: np.ndarray | float = 1.0,
    omega: np.ndarray | None = None,
    max_iter: int = 8,
    tol: float = 1e-4,
    solver: str = "SCS",
    verbose: bool = False,
) -> tuple[np.ndarray, float]:
    """IRLS / Newton on Poisson NLL.

    Start uniform, predict μ, weight 1/sqrt(μ+ε), solve weighted-L2, refit.
    """
    K = counts.shape[0]
    F = cn.shape[0]
    if omega is None:
        omega = np.zeros(K)
    if np.isscalar(alpha):
        alpha = np.full(K, alpha)

    h_curr = np.full(F, 1.0 / F)
    obj_val = np.inf

    for it in range(max_iter):
        mu_curr = coverage * alpha * (cn.T @ h_curr) + coverage * omega
        weights = 1.0 / np.sqrt(np.maximum(mu_curr, 1e-3))

        h_new, obj_val = solve_block_wls(
            counts, cn, coverage, alpha=alpha, omega=omega,
            weights=weights, solver=solver, verbose=False
        )
        diff = np.linalg.norm(h_new - h_curr)
        h_curr = h_new
        if verbose:
            print(f"  IRLS iter {it}: ||Δh||={diff:.2e}, obj={obj_val:.4f}")
        if diff < tol:
            break

    return h_curr, obj_val


def update_omega(
    counts_all: np.ndarray,        # N × K matrix of counts (one row per sample)
    cn: np.ndarray,                # F × K matrix
    h_all: np.ndarray,             # N × F matrix of fitted h per sample (within this block)
    coverage: np.ndarray,          # N-vector of λ_i per sample
    alpha: np.ndarray | float = 1.0,
) -> np.ndarray:
    """Estimate per-k-mer contamination rate ω_k from across-sample residuals.

    For each k-mer, ω_k is the average residual count not explained by founder
    evidence, normalized by per-sample coverage. Clipped at 0 from below.
    """
    K = counts_all.shape[1]
    if np.isscalar(alpha):
        alpha = np.full(K, alpha)

    pred_signal = (h_all @ cn) * coverage[:, None] * alpha[None, :]       # N × K
    residual = counts_all - pred_signal
    omega = np.maximum(0, residual.mean(axis=0) / coverage.mean())
    return omega
