"""
Within-bubble correlation fix for the Poisson k-mer founder-frequency EM.

PROBLEM
-------
The pool model is  c_k ~ Poisson(lambda * (h^T cn)_k), and the EM treats each of
the K ~ 11M k-mers as an independent Poisson observation. But k-mers inside one
"bubble" (a local variant cluster, given by meta['bubble_id']) are physically
LINKED: if a founder haplotype is present in the pool, ALL of that haplotype's
k-mers in the bubble fire together. So a bubble with m_b k-mers is NOT m_b
independent observations -- it is closer to ~1 (or, more precisely, ~#distinct
carrier patterns) independent observation.

Because the cactus founders have DENSER bubbles (more k-mers per bubble), the
naive likelihood gives the cactus founders disproportionately more effective
sample size, inflating their estimated mass (real-SEEDMIX cactus 0.518 vs
truth ~0.35).

METHOD IMPLEMENTED: (c2) OVERDISPERSION / EFFECTIVE-n DOWN-WEIGHTING
--------------------------------------------------------------------
We DO NOT touch the design (mu_k = lambda * (h^T cn)_k is unchanged) and we DO
NOT delete any data. Instead we down-weight the LIKELIHOOD contribution of every
k-mer so that a size-m_b bubble contributes m_eff effective independent
observations rather than m_b near-duplicates, where

    m_eff = m_b / (1 + (m_b - 1) * rho)      (the classic "design effect" for
                                              an exchangeable/equicorrelated
                                              cluster with intra-cluster
                                              correlation rho)

    omega_k = m_eff / m_b = 1 / (1 + (m_b - 1) * rho)   (per-k-mer weight,
                                                          identical for all
                                                          k-mers in a bubble)

This is a weighted Poisson EM. With per-observation weights omega_k the M-step
multiplicative update is

    h_f <- h_f * ( sum_k omega_k * cn_fk * c_k / mu_k )
                 / ( sum_k omega_k * cn_fk )

then renormalize to the simplex. (This is the standard weighted-EM /
fractional-count generalization: every place the unweighted update sums c_k or
counts a carried k-mer, we instead sum omega_k * c_k and omega_k. mu_k is
computed from the UNWEIGHTED h^T cn, i.e. the model rate is unchanged -- omega_k
only scales how much each k-mer's residual c_k/mu_k pulls on h.)

WHY c2 (not c1): it is the safer first implementation -- no clustering of
carrier patterns, no aggregation bookkeeping, exact reproduction of the existing
EM at rho=0 (omega==1), and it lets us SWEEP rho in {0.5, 0.9, 0.99} to trace
how much the cactus inflation comes from the within-bubble pseudo-replication.
rho=1 would collapse each bubble to exactly 1 effective obs (omega_k = 1/m_b);
rho->0 recovers plain filt2.

INTUITION ON THE BIAS DIRECTION: a cactus-dense bubble (large m_b) gets a SMALL
omega_k, so its many linked k-mers stop counting as many independent votes for
the cactus founders -- exactly the over-credit we want to remove. PG/shared
bubbles tend to be sparser, so they are down-weighted less, relatively raising
their voice. No founder's k-mers are deleted, so identifiability on the
cactus-heavy sims is preserved (unlike subsampMedian).

This script mirrors poolfreq/src/sweep_shape_norm_h_only.py for cn/meta load,
nz-count filtering, densify, lambda/cov bookkeeping, and simplex output. It
reuses the precomputed filt2-indexed read counts (NO jellyfish).
"""
from __future__ import annotations
import argparse, gc, json, os, sys, time
import numpy as np
from scipy.sparse import load_npz


def weighted_solve_em(counts, cn_dense, omega, max_iter=200, tol=1e-7,
                      h_init=None, verbose=False):
    """Weighted multinomial-mixture (Poisson) EM on the simplex.

    counts:    K-vector of observed k-mer counts c_k  (float32)
    cn_dense:  F x K binary copy-number, UNWEIGHTED (float32, design unchanged)
    omega:     K-vector of per-k-mer likelihood weights (float32)
    Returns (h float64, info dict).

    This is the weighted generalization of poolfreq/src/em_solver.solve_em.
    The model rate mu_k = (h^T cn)_k is UNWEIGHTED (the design is untouched);
    omega only scales how much each k-mer's count mass contributes to the
    M-step. Concretely we replace the observed count mass c_k by the EFFECTIVE
    mass omega_k * c_k everywhere in the E/M step:

        mu_k  = max(h^T cn, eps)                  # responsibilities use raw rate
        cw_k  = omega_k * counts_k / mu_k         # weighted count residual
        num_f = h_f * (cn @ cw)                   # sum_k omega cn c/mu
        h_f  <- num_f / sum_f' num_f'             # == num_f / sum_k omega_k c_k
        h    <- h / h.sum()

    The normalizer sum_f num_f = sum_k omega_k c_k (the total effective count
    mass), because sum_f h_f cn_fk / mu_k == 1 for every observed k. At
    omega==1 this is EXACTLY em_solver.solve_em's update (numerator h*(cn@(c/mu)),
    normalized by total_c), so rho=0 reproduces the filt2 baseline bit-for-bit.

    A bubble of size m_b contributes total effective mass omega * (its count
    mass) = (1/(1+(m_b-1)rho)) * mass, i.e. m_eff/m_b of its raw weight -- the
    intended design-effect down-weighting, applied to the LIKELIHOOD only.
    """
    F, K = cn_dense.shape
    counts = counts.astype(np.float32)
    omega = omega.astype(np.float32)
    h = (np.full(F, 1.0 / F, dtype=np.float32)
         if h_init is None else h_init.astype(np.float32))

    weighted_counts = omega * counts            # effective count mass per k-mer
    total_eff = float(weighted_counts.sum())    # = sum_f num_f at every iter

    history = []
    delta = np.inf
    for it in range(max_iter):
        mu = np.maximum(h @ cn_dense, np.float32(1e-7))     # (K,) UNWEIGHTED rate
        cw = weighted_counts / mu                           # (K,)
        num = h * (cn_dense @ cw)                           # (F,)
        h_new = num / max(total_eff, 1e-12)
        h_new = h_new / h_new.sum()
        delta = float(np.linalg.norm(h_new - h))
        history.append(delta)
        if verbose and it % 20 == 0:
            print(f"    EM iter {it}: ||dh||={delta:.2e} "
                  f"h_max={h_new.max():.4f} eff_n={1/(h_new**2).sum():.1f}",
                  flush=True)
        h = h_new
        if delta < tol:
            break
    return h.astype(np.float64), {"iterations": it + 1,
                                  "converged": delta < tol,
                                  "delta_history": history}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cn-prefix", required=True,
                    help="expects {prefix}.cn.npz + {prefix}.meta.npz")
    ap.add_argument("--counts", required=True,
                    help="precomputed filt2-indexed counts .npy (K-vector)")
    ap.add_argument("--sample", required=True)
    ap.add_argument("--out", required=True, help="output .npz path")
    ap.add_argument("--rhos", default="0.0,0.5,0.9,0.99",
                    help="comma-separated intra-bubble correlation rho values")
    ap.add_argument("--em-max-iter", type=int, default=300)
    args = ap.parse_args()

    rhos = [float(r) for r in args.rhos.split(",")]
    print(f"=== correlation (c2 overdispersion) sample={args.sample} "
          f"rhos={rhos} ===", flush=True)
    t0 = time.time()

    cn = load_npz(args.cn_prefix + ".cn.npz").tocsr()
    meta = np.load(args.cn_prefix + ".meta.npz", allow_pickle=True)
    founders = np.asarray(meta["founders"]).astype(str)
    bubble_id = np.asarray(meta["bubble_id"]).astype(np.int64)
    F, K = cn.shape
    print(f"  cn: F={F} K={K:,} nnz={cn.nnz:,}", flush=True)

    counts = np.load(args.counts)
    assert counts.shape[0] == K, f"counts K={counts.shape[0]} != cn K={K}"
    print(f"  counts: nz={(counts>0).sum():,}/{K:,} sum={counts.sum():,} "
          f"max={counts.max()}", flush=True)

    # Per-k-mer bubble size m_b (same as np.bincount(bubble_id)[bubble_id]).
    mb = np.bincount(bubble_id, minlength=int(bubble_id.max()) + 1)[bubble_id]
    mb = mb.astype(np.float64)
    print(f"  m_b per k-mer: median={np.median(mb):.0f} mean={mb.mean():.1f} "
          f"max={mb.max():.0f}", flush=True)

    # Densify cn (float32), mirror sweep_shape_norm_h_only.
    t = time.time()
    cn_dense = np.asarray(cn.todense()).astype(np.float32)
    print(f"  densify cn -> {cn_dense.nbytes/1e9:.1f} GB  [{time.time()-t:.0f}s]",
          flush=True)
    del cn
    gc.collect()

    # nz filter (smaller EM matrix) -- IDENTICAL ordering to the mirror.
    nz = counts > 0
    cn_em = np.ascontiguousarray(cn_dense[:, nz])
    counts_em = counts[nz].astype(np.float32)
    mb_em = mb[nz]
    cov = counts.sum() * F / max(cn_dense.sum(), 1)   # lambda/cov bookkeeping
    print(f"  EM matrix: {cn_em.shape}, counts_sum={counts_em.sum():.0f}, "
          f"cov~{cov:.3f}", flush=True)
    del cn_dense
    gc.collect()

    n_rho = len(rhos)
    h_per_rho = np.zeros((n_rho, F), dtype=np.float64)
    info_per_rho = []
    for ir, rho in enumerate(rhos):
        # omega_k = 1 / (1 + (m_b - 1) * rho).  rho=0 -> omega=1 (plain filt2).
        omega = (1.0 / (1.0 + (mb_em - 1.0) * rho)).astype(np.float32)
        t = time.time()
        h, info = weighted_solve_em(counts_em, cn_em, omega,
                                    max_iter=args.em_max_iter, tol=1e-7)
        h_per_rho[ir] = h
        eff_n = float(1.0 / (h ** 2).sum())
        info_per_rho.append({"rho": float(rho),
                             "iterations": int(info["iterations"]),
                             "converged": bool(info["converged"]),
                             "omega_min": float(omega.min()),
                             "omega_med": float(np.median(omega)),
                             "h_eff_n": eff_n})
        print(f"  rho={rho:>5.2f}: omega in [{omega.min():.4f},1] "
              f"med={np.median(omega):.4f}, EM iter={info['iterations']}, "
              f"eff_n={eff_n:.1f}, h_max={h.max():.4f}  [{time.time()-t:.0f}s]",
              flush=True)

    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    np.savez(args.out,
             founders=founders,
             rhos=np.array(rhos),
             h_per_rho=h_per_rho,
             info=np.array(info_per_rho, dtype=object),
             variant="c2_overdispersion",
             sample=args.sample)
    print(f"\nWrote {args.out}", flush=True)
    print(f"  TOTAL {time.time()-t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
