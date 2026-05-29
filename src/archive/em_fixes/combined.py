"""
Combined per-bubble x inverse-AC LIKELIHOOD weighting EM (h-only).

Discrete MuSiC / W-NNLS analogue. Each k-mer's LIKELIHOOD TERM is weighted by

    omega_k = 1 / (m_b * a_k^gamma)

where
    m_b = #k-mers in k's bubble  (per-bubble pseudo-replication de-weighting)
    a_k = #founders carrying k (>=2 after ac=1 drop)  (inverse-AC down-weight of
          the low-AC k-mers that dominate the Fisher information and are 2.4x
          more abundant in cactus founders).

The design rate mu_k = lambda * (h^T cn)_k is UNCHANGED. We do NOT scale cn rows,
and there is NO post-EM h/w recovery step. We only re-weight the contribution of
each k-mer to the M-step sufficient statistics:

    h_f <- h_f * ( sum_k omega_k * cn_fk * c_k / mu_k ) / ( sum_k omega_k * cn_fk )
    renormalize h to the simplex.

This is the standard weighted-likelihood Poisson EM. For any FIXED omega_k it is
a fixed-point map whose stationary point is the weighted-MLE; it is consistent
(does not delete data, never escapes the simplex interior from a uniform start).

gamma = 0  -> per-bubble weighting only.

Mirrors src/em_solver.py (lambda/coverage handling, nz filter, densify,
simplex normalization) and the sweep skeleton of
src/archive/sweep_shape_norm_h_only.py.

Inputs are PRE-COMPUTED filt2-indexed read counts (no jellyfish): reuse the
existing scratch/g0_sweep_h_test/filt2_*.counts.npy and
scratch/seedmix_h_test/filt2_S*.counts.npy.

Outputs:
  <out_prefix>.h_sweep.npz  keys: founders, gammas,
      h_per_gamma (n_gamma x F), info (list of dicts)
"""
from __future__ import annotations
import argparse, gc, os, sys, time
import numpy as np
from scipy.sparse import load_npz


def solve_em_weighted(
    counts: np.ndarray,      # K-vector observed counts (nz-filtered)
    cn: np.ndarray,          # F x K float32 (binary cn, nz-filtered) -- NOT scaled
    omega: np.ndarray,       # K-vector likelihood weights (nz-filtered)
    max_iter: int = 200,
    tol: float = 1e-7,
    verbose: bool = False,
):
    """Weighted-likelihood multiplicative Poisson EM.

    mu_k = (h^T cn)_k  (lambda cancels in the multiplicative M-step ratio).
    M-step:
        num_f   = sum_k omega_k * cn_fk * (c_k / mu_k)
        den_f   = sum_k omega_k * cn_fk
        h_f    <- h_f * num_f / den_f      (then renormalize to simplex)
    """
    F, K = cn.shape
    cn = cn.astype(np.float32, copy=False)
    counts = counts.astype(np.float32)
    omega = omega.astype(np.float32)

    # weighted counts (omega_k * c_k) reused every iter
    wc = (omega * counts).astype(np.float32)
    # constant per-founder denominator: sum_k omega_k * cn_fk
    den = (cn @ omega).astype(np.float32)           # F-vec
    den = np.maximum(den, np.float32(1e-12))

    h = np.full(F, 1.0 / F, dtype=np.float32)
    history = []
    it = 0
    for it in range(max_iter):
        mu = np.maximum(h @ cn, np.float32(1e-7))   # K-vec
        # num_f = sum_k cn_fk * (omega_k * c_k / mu_k)
        num = h * (cn @ (wc / mu))                  # F-vec  (h_f * sum cn_fk * wc_k/mu_k)
        h_new = num / den
        h_new = h_new / h_new.sum()
        delta = float(np.linalg.norm(h_new - h))
        history.append(delta)
        if verbose and it % 20 == 0:
            print(f"    EM iter {it}: ||dh||={delta:.2e} hmax={h_new.max():.4f}", flush=True)
        h = h_new
        if delta < tol:
            break
    return h.astype(np.float64), {"iterations": it + 1,
                                  "converged": delta < tol,
                                  "delta_history": history}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cn-prefix", required=True,
                    help="prefix; expects {prefix}.cn.npz and {prefix}.meta.npz")
    ap.add_argument("--counts", required=True,
                    help="precomputed filt2-indexed counts .npy (K-vector)")
    ap.add_argument("--sample", required=True)
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--gammas", default="0,0.25,0.5,0.75,1.0",
                    help="comma-separated gamma values for a_k^gamma")
    ap.add_argument("--em-max-iter", type=int, default=200)
    args = ap.parse_args()

    gammas = [float(g) for g in args.gammas.split(",")]
    print(f"=== combined per-bubble x inverse-AC EM: sample={args.sample} "
          f"gammas={gammas} ===", flush=True)
    t0 = time.time()

    cn = load_npz(args.cn_prefix + ".cn.npz").tocsr()
    meta = np.load(args.cn_prefix + ".meta.npz", allow_pickle=True)
    founders = np.asarray(meta["founders"]).astype(str)
    bubble_id = np.asarray(meta["bubble_id"]).astype(np.int64)
    F, K = cn.shape
    print(f"  cn: F={F} K={K:,} nnz={cn.nnz:,}", flush=True)

    counts = np.load(args.counts)
    assert counts.shape[0] == K, f"counts K={counts.shape[0]} != cn K={K}"
    print(f"  counts: nz={(counts > 0).sum():,}/{K:,} sum={counts.sum():,.0f}", flush=True)

    # ---- per-kmer weight components (full K, on the panel) ----
    # m_b = #k-mers in k's bubble
    t = time.time()
    bubble_size = np.bincount(bubble_id)           # per bubble
    m_b = bubble_size[bubble_id].astype(np.float64)  # per-kmer
    # a_k = #founders carrying k  (>=2 since ac=1 already dropped)
    a_k = np.asarray(cn.sum(axis=0)).flatten().astype(np.float64)
    print(f"  m_b: med={np.median(m_b):.1f} max={m_b.max():.0f}; "
          f"a_k: med={np.median(a_k):.1f} min={a_k.min():.0f} max={a_k.max():.0f} "
          f"[{time.time()-t:.0f}s]", flush=True)

    # ---- densify, nz filter (mirror reference) ----
    t = time.time()
    cn_dense = np.asarray(cn.todense()).astype(np.float32)
    print(f"  densify -> {cn_dense.nbytes/1e9:.1f} GB float32 [{time.time()-t:.0f}s]", flush=True)
    del cn; gc.collect()

    nz = counts > 0
    cn_em = np.ascontiguousarray(cn_dense[:, nz])
    counts_em = counts[nz].astype(np.float32)
    m_b_em = m_b[nz]
    a_k_em = a_k[nz]
    del cn_dense; gc.collect()
    print(f"  EM matrix: {cn_em.shape}", flush=True)

    h_per_gamma = np.zeros((len(gammas), F), dtype=np.float64)
    info_per_gamma = []
    for ig, g in enumerate(gammas):
        # omega_k = 1 / (m_b * a_k^gamma)
        omega = 1.0 / (m_b_em * np.power(a_k_em, g))
        t = time.time()
        h, info = solve_em_weighted(counts_em, cn_em, omega,
                                    max_iter=args.em_max_iter, tol=1e-7)
        h_per_gamma[ig] = h
        info_per_gamma.append({
            "gamma": float(g),
            "iterations": int(info["iterations"]),
            "converged": bool(info["converged"]),
            "eff_n": float(1.0 / (h ** 2).sum()),
            "omega_min": float(omega.min()), "omega_max": float(omega.max()),
        })
        print(f"  gamma={g:>4.2f}: EM iter={info['iterations']} "
              f"eff_n={1/(h**2).sum():.1f} hmax={h.max():.4f} [{time.time()-t:.0f}s]",
              flush=True)

    out_path = args.out_prefix + ".h_sweep.npz"
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    np.savez(out_path,
             founders=founders,
             gammas=np.array(gammas),
             h_per_gamma=h_per_gamma,
             info=np.array(info_per_gamma, dtype=object))
    print(f"\nWrote {out_path}  TOTAL {time.time()-t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
