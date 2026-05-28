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
    dirichlet_alpha: float = 0.0,  # symmetric Dirichlet prior strength on h
    prior_h: np.ndarray | None = None,    # NEW: anchor h to prior (per-window EM)
    prior_weight: float = 0.0,            # NEW: λ — strength of anchor toward prior_h
    omega: np.ndarray | None = None,      # NEW: per-k-mer weight ω_k (e.g. 1/m_b)
) -> tuple[np.ndarray, dict]:
    """Run EM until h converges.

    With dirichlet_alpha > 0, applies a symmetric Dirichlet(α) prior on h.
    The MAP update becomes:
        h_new[f] ∝ (h[f] · Σ_k cn[f,k] · counts[k] / μ_k) + (α - 1) / total_c
    For α = 1: no regularization (uniform prior, equivalent to MLE).
    For α > 1: pulls h toward uniform, prevents winner-takes-all collapse.
    Use small values like α = 1.01 to 1.5 for mild regularization.

    With prior_h + prior_weight > 0, applies a Dirichlet pseudocount centered
    on prior_h (instead of uniform):
        h_new[f] ∝ (h[f] · Σ_k cn[f,k] · counts[k] / μ_k) + λ · total_c · prior_h[f]
    This is the MAP update under a Dirichlet(α_f = 1 + λ·N·prior_h[f]) prior —
    pulls the per-window solution toward `prior_h` (typically the chrom-wide
    `h_global`). Local k-mer evidence has to overcome the prior to move h away
    from prior_h. λ=0 → pure MLE; λ=1 → prior is as influential as the data;
    sweet spot is typically λ ∈ [0.05, 0.5].

    With omega (K-vector of per-k-mer weights ω_k), every count is reweighted:
    the M-step uses ω_k·c_k everywhere it used c_k, i.e.
        h_new[f] ∝ h[f] · Σ_k cn[f,k] · (ω_k·c_k)/μ_k ,  normalized by Σ_k ω_k·c_k.
    ω_k = 1/m_b (m_b = #k-mers in k's bubble) is per-bubble de-replication; it
    turns "h ∝ k-mer count" into "h ∝ locus count" and removes the imbalanced-
    design over-credit. omega=None reproduces the unweighted MLE exactly.

    Returns (h, info_dict).
    """
    K = counts.shape[0]
    F = cn.shape[0]
    if cn.dtype != np.float32:
        cn = cn.astype(np.float32)
    counts = counts.astype(np.float32)

    h = np.full(F, 1.0 / F, dtype=np.float32) if h_init is None \
        else h_init.astype(np.float32)
    # ω-weighted counts (omega=None → unweighted, identical to MLE)
    wc = counts if omega is None else (omega.astype(np.float32) * counts)
    total_c = wc.sum()

    # Dirichlet prior pseudo-count (added to numerator before normalization)
    prior_pseudo = np.float32(max(0.0, dirichlet_alpha - 1.0))

    # Anchor pseudocount toward prior_h (per-window anchor toward h_global)
    use_anchor = prior_h is not None and prior_weight > 0 and total_c > 0
    if use_anchor:
        anchor_term = (np.float32(prior_weight) * np.float32(total_c)
                       * prior_h.astype(np.float32))
    else:
        anchor_term = None

    history = []
    for it in range(max_iter):
        denom = np.maximum(h @ cn, np.float32(1e-7))
        cw = wc / denom
        em_term = h * (cn @ cw)
        if anchor_term is not None:
            em_term = em_term + anchor_term
        if prior_pseudo > 0:
            h_new = (em_term + prior_pseudo) / (total_c + F * prior_pseudo
                                                 + (anchor_term.sum() if anchor_term is not None else 0))
        else:
            denom_norm = total_c
            if anchor_term is not None:
                denom_norm = denom_norm + anchor_term.sum()
            h_new = em_term / max(denom_norm, np.float32(1e-12))
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
