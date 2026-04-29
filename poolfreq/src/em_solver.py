"""
EM solver for pool-seq founder frequencies.

Mixture-model formulation:
    For each k-mer k, the observed count c_k is the sum of contributions from
    each founder f, weighted by h[f] · cn[f,k] · λ. The expected total rate is
    μ_k = λ · Σ_f h[f] · cn[f,k].

EM treats this as a multinomial mixture:
    Each "k-mer event" (each unit of count) comes from one founder f, with
    probability  P(f | k) = h[f] · cn[f,k] / Σ_f' h[f'] · cn[f',k].

E-step:
    For each k-mer, distribute its count fractionally over founders carrying it:
        n_attributed[f, k] = c_k · h[f] · cn[f,k] / sum_f'(h[f'] · cn[f',k])

M-step:
    Update h[f] ∝ Σ_k n_attributed[f, k]
    (equivalent to weighted average of count-mass attributed to founder f)

This iteration:
  - never escapes the open interior of the simplex if started uniform
  - is monotone-improving in the Poisson likelihood
  - is HARP / PanGenie's same idea, generalized for k-mer counts vs reads

Reference: Kessner et al. 2013 (HARP); Dempster, Laird, Rubin 1977 (EM).
"""
from __future__ import annotations
import numpy as np


def solve_em(
    counts: np.ndarray,            # K-vector of observed counts
    cn: np.ndarray,                # F × K binary copy-number matrix
    coverage: float,               # not used directly in EM (it cancels in the M-step)
    h_init: np.ndarray | None = None,
    max_iter: int = 100,
    tol: float = 1e-6,
    verbose: bool = False,
) -> tuple[np.ndarray, dict]:
    """Run EM until h converges.

    Returns (h, info_dict).
    """
    K = counts.shape[0]
    F = cn.shape[0]
    # Keep everything float32 — mixing dtypes triggers an in-place upcast of
    # cn (F × K) to float64 inside numpy's matmul, which is the fast path
    # killer (>10× slowdown at K~10M+).
    if cn.dtype != np.float32:
        cn = cn.astype(np.float32)
    counts = counts.astype(np.float32)

    h = np.full(F, 1.0 / F, dtype=np.float32) if h_init is None \
        else h_init.astype(np.float32)
    total_c = counts.sum()

    history = []
    for it in range(max_iter):
        denom = np.maximum(h @ cn, np.float32(1e-7))
        cw = counts / denom
        h_new = h * (cn @ cw) / total_c
        h_new = h_new / h_new.sum()

        delta = float(np.linalg.norm(h_new - h))
        history.append(delta)
        if verbose and it % 10 == 0:
            print(f"  EM iter {it}: ||Δh||={delta:.2e}, h_min={h_new.min():.4f}, "
                  f"h_max={h_new.max():.4f}")
        h = h_new
        if delta < tol:
            if verbose:
                print(f"  EM converged at iter {it}: ||Δh||={delta:.2e}")
            break

    return h.astype(np.float64), {"iterations": it + 1,
                                  "delta_history": history,
                                  "converged": delta < tol}


def solve_em_with_omega(
    counts: np.ndarray,
    cn: np.ndarray,
    coverage: float,
    omega_init: np.ndarray | None = None,
    h_init: np.ndarray | None = None,
    max_iter: int = 100,
    tol: float = 1e-6,
):
    """EM with per-k-mer contamination ω. Adds a "ghost founder" representing
    background contamination that doesn't depend on h.

    Hidden mixture: each count is from one of (founders) or from "contamination"
    at rate ω_k. Update ω jointly with h.
    """
    K = counts.shape[0]
    F = cn.shape[0]
    if cn.dtype != np.float32:
        cn = cn.astype(np.float32)
    counts = counts.astype(np.float64)

    h = np.full(F, 1.0 / F) if h_init is None else h_init.copy()
    omega = np.zeros(K) if omega_init is None else omega_init.copy()

    for it in range(max_iter):
        # E-step: per-k-mer denominator includes ω
        # rate from founders: λ · h^T · cn[:, k]
        # rate from contam:  λ · ω[k]  (independent of h)
        # Posterior P(founder f | k) ∝ h[f] · cn[f,k]; P(contam | k) ∝ ω[k]
        # Total: h^T cn + ω
        rate_f = h @ cn  # K-vector
        denom = rate_f + omega + 1e-15

        # Attribute count to founders vs contamination
        c_to_founders = counts * rate_f / denom  # K-vector — total mass attributed to founders at k
        c_to_contam = counts * omega / denom  # K-vector

        # M-step for h (same as before, but with founder-attributed counts)
        cw = c_to_founders / np.maximum(rate_f, 1e-15)
        # ^ this is h-independent denominator; equivalent: cw = c_to_founders / rate_f
        # = counts × h^T cn / denom / rate_f = counts / denom (when rate_f > 0)
        # Use simpler form
        cw_simple = counts / denom
        h_new = h * (cn @ cw_simple) / max(c_to_founders.sum(), 1e-15)
        h_new = h_new / h_new.sum()

        # M-step for ω: ω_k = c_to_contam[k] / coverage
        # (per-k-mer contamination rate, normalized by coverage)
        omega_new = c_to_contam / max(coverage, 1e-9)

        delta_h = np.linalg.norm(h_new - h)
        delta_o = np.linalg.norm(omega_new - omega)
        h = h_new
        omega = omega_new
        if delta_h < tol and delta_o < tol * coverage:
            break

    return h, omega, {"iterations": it + 1}
