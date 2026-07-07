"""Uncertainty / confidence on the founder-mixture estimate ĥ.

Two estimators of the sampling covariance of the EM solution ĥ on the simplex:

1. fisher_cov_h  — analytic, O(F²·n_nz): the observed Fisher information of the
   (weighted) Poisson likelihood, inverted on the simplex tangent space.
       J = K diag(ω_k c_k / μ_k²) Kᵀ,   μ_k = Σ_f h_f K_{f,k}
   Cheap enough to run per sample (and per window). Captures BOTH
     - low coverage   → small c_k → small J → wide SE, and
     - complex/collinear regions → near-singular J → wide SE on the
       non-identifiable founder directions.

2. bootstrap_cov_h — gold-standard but B× the cost: a parametric Poisson
   bootstrap. Resimulate counts c* ~ Poisson(λ̂ μ_k(ĥ)) at the fitted ĥ
   (depth-matched to the observed total), refit the EM B times, take the
   empirical covariance of the refits. Use to validate that the analytic
   Fisher SEs are calibrated.

Both return a full F×F covariance Σ that is zero outside the estimated support
S = {f : ĥ_f > support_eps}: founders the pool does not contain sit on the
simplex boundary where the interior asymptotics do not apply, so their variance
is reported as 0 (they are flagged via `support`, not given a spurious SE).

Per-founder SE = sqrt(diag(Σ)). To propagate to a per-record AF SE, use the
delta method on the projection AF_r = (ĥᵀv_r)/(ĥᵀu_r):
    g_r = (v_r - AF_r·u_r) / (ĥᵀu_r),   Var(AF_r) = g_rᵀ Σ g_r
(see af_se_from_cov below).
"""
from __future__ import annotations
import numpy as np

from .em_solver import solve_em


def fisher_information_h(h, kmer_pa, counts, omega=None, eps_mu=1e-7, chunk=None):
    """Observed Fisher information J = K diag(ω_k c_k / μ_k²) Kᵀ at h (F×F).

    kmer_pa: F×K (dense). counts: K. omega: K or None (ω_k≡1).
    Only k-mers with c_k>0 contribute (zeros drop out), matching the EM.

    chunk: if set, accumulate J over column blocks of this width instead of
    materializing the full F×K float64 temporary `kmer_pa*w` — needed on large
    panels (e.g. 231 founders × millions of k-mers, where the temp is ~10 GB).
    """
    h = np.asarray(h, dtype=np.float64)
    kmer_pa = np.asarray(kmer_pa)
    counts = np.asarray(counts, dtype=np.float64)
    mu = np.maximum(h @ kmer_pa, eps_mu)          # K
    w = counts / (mu * mu)                         # c_k / μ_k²
    if omega is not None:
        w = w * np.asarray(omega, dtype=np.float64)
    # J = kmer_pa @ diag(w) @ kmer_paᵀ = (kmer_pa * w) @ kmer_paᵀ
    if chunk is None:
        Kw = kmer_pa.astype(np.float64) * w        # F×K (broadcast over columns)
        return Kw @ kmer_pa.astype(np.float64).T   # F×F
    F, Kn = kmer_pa.shape
    J = np.zeros((F, F), dtype=np.float64)
    for s0 in range(0, Kn, chunk):
        sl = slice(s0, min(s0 + chunk, Kn))
        Kc = kmer_pa[:, sl].astype(np.float64)
        J += (Kc * w[sl]) @ Kc.T
    return J


def _tangent_pinv_on_support(J, support, rcond=1e-10):
    """(F×F) covariance: pseudo-inverse of J restricted to `support`, computed
    on the support-simplex tangent space (Σ_{f∈S} h_f = 1). Zero elsewhere."""
    F = J.shape[0]
    Sigma = np.zeros((F, F), dtype=np.float64)
    s = np.asarray(support, dtype=int)
    if s.size <= 1:
        return Sigma                              # nothing resolvable
    m = s.size
    Jss = J[np.ix_(s, s)]
    P = np.eye(m) - np.ones((m, m)) / m            # tangent projector
    Jt = P @ Jss @ P
    Sig_ss = P @ np.linalg.pinv(Jt, rcond=rcond) @ P
    Sigma[np.ix_(s, s)] = Sig_ss
    return Sigma


def fisher_cov_h(h, kmer_pa, counts, omega=None, support_eps=1e-3, rcond=1e-2,
                 chunk=None):
    """Analytic covariance Σ (F×F) of ĥ from the observed Fisher information.

    Returns (Sigma, support). SE per founder = sqrt(diag(Sigma)).

    Defaults (support_eps=1e-3, rcond=1e-2) are tuned so the analytic SE tracks
    the parametric bootstrap: a loose support or small rcond lets near-boundary /
    collinear founder directions (tiny Fisher eigenvalues) blow the variance up
    10-20×, because the unconstrained Gaussian asymptotics break where the
    multiplicative EM actually pins those founders. rcond drops those
    unresolvable directions, matching the bootstrap. The right rcond is mildly
    data-dependent — see benchmarks/h_uncertainty.
    """
    J = fisher_information_h(h, kmer_pa, counts, omega=omega, chunk=chunk)
    support = np.flatnonzero(np.asarray(h) > support_eps)
    Sigma = _tangent_pinv_on_support(J, support, rcond=rcond)
    return Sigma, support


def identifiability(h, kmer_pa, counts, omega=None, support_eps=1e-6):
    """Scalar region-resolvability diagnostics from the support Fisher info:
    effective rank and condition number of J_SS. High condition number / low
    eff-rank ⇒ a complex / collinear region that cannot resolve the mixture."""
    J = fisher_information_h(h, kmer_pa, counts, omega=omega)
    s = np.flatnonzero(np.asarray(h) > support_eps)
    if s.size <= 1:
        return {"support_size": int(s.size), "eff_rank": float(s.size),
                "cond": np.inf}
    ev = np.linalg.eigvalsh(J[np.ix_(s, s)])
    ev = np.clip(ev, 0, None)
    pos = ev[ev > 0]
    # effective rank = exp(entropy of normalized eigenvalue spectrum)
    p = pos / pos.sum()
    eff_rank = float(np.exp(-(p * np.log(p)).sum())) if pos.size else 0.0
    cond = float(pos.max() / pos.min()) if pos.size else np.inf
    return {"support_size": int(s.size), "eff_rank": eff_rank, "cond": cond}


def bootstrap_cov_h(h_hat, kmer_pa, counts, omega=None, B=200,
                    coverage=None, seed=0, max_iter=200, tol=1e-7,
                    rate_h=None, return_samples=False, normalize="per_founder"):
    """Parametric Poisson bootstrap covariance of ĥ.

    Resimulate c* ~ Poisson(λ̂ μ_k(rate_h)) — depth-matched to Σcounts — and
    refit the EM B times. rate_h defaults to h_hat (standard bootstrap); pass
    the TRUE h to get the gold Monte-Carlo sampling distribution instead.
    Returns (Sigma, samples|None). SE = sqrt(diag(Sigma)).
    """
    rng = np.random.default_rng(seed)
    h_hat = np.asarray(h_hat, dtype=np.float64)
    kmer_pa = np.asarray(kmer_pa)
    counts = np.asarray(counts, dtype=np.float64)
    rh = h_hat if rate_h is None else np.asarray(rate_h, dtype=np.float64)
    mu = rh @ kmer_pa
    lam = counts.sum() / max(mu.sum(), 1e-12)      # depth-match expected total
    rate = lam * mu
    F = h_hat.size
    samples = np.empty((B, F), dtype=np.float64)
    for b in range(B):
        c_star = rng.poisson(rate).astype(np.float32)
        nz = c_star > 0
        om = None if omega is None else np.asarray(omega)[nz]
        h_b, _ = solve_em(c_star[nz], kmer_pa[:, nz], coverage or lam,
                          max_iter=max_iter, tol=tol, omega=om, normalize=normalize)
        samples[b] = h_b
    Sigma = np.cov(samples.T)
    return (Sigma, samples if return_samples else None)


def af_se_from_cov(h, Sigma, var_pa_col, var_called_col, support=None,
                   eps_den=1e-12):
    """Delta-method per-record AF SE from the founder covariance Σ.

    AF_r = (hᵀv_r)/(hᵀu_r);  g_r = (v_r - AF_r·u_r)/(hᵀu_r);  Var = g_rᵀ Σ g_r.

    var_pa_col, var_called_col: F×R dense (or array-like) columns v, u for the R
    records to score. Returns (af, se) length-R. If support given, Σ and the
    gradient are restricted to it (faster; Σ is zero off-support anyway).
    """
    h = np.asarray(h, dtype=np.float64)
    V = np.asarray(var_pa_col, dtype=np.float64)        # F×R
    U = np.asarray(var_called_col, dtype=np.float64)    # F×R
    den = np.maximum(h @ U, eps_den)                    # R
    af = (h @ V) / den                                  # R
    G = (V - af[None, :] * U) / den[None, :]            # F×R gradient
    if support is not None and len(support) and len(support) < h.size:
        s = np.asarray(support, dtype=int)
        Gs = G[s, :]
        var = np.einsum("ir,ij,jr->r", Gs, Sigma[np.ix_(s, s)], Gs)
    else:
        var = np.einsum("ir,ij,jr->r", G, Sigma, G)
    se = np.sqrt(np.clip(var, 0, None))
    return af, se
