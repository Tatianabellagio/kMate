"""
Per-bubble normalization fix for the Poisson k-mer founder-frequency EM.

PROBLEM
-------
h (231-founder simplex) is estimated by c_k ~ Poisson(lambda * (h^T cn)_k).
The filt2 cn_full over-credits the 80 "cactus" founders because of within-bubble
pseudo-replication: linked k-mers inside one bubble are counted m_b times instead
of as ~1 independent observation, and m_b is larger for cactus founders.

FIX = PER-BUBBLE NORMALIZATION (no tuning parameter)
----------------------------------------------------
Weight each k-mer's likelihood contribution by
    omega_k = 1 / m_b ,   m_b = number of k-mers in k's bubble.
This makes each genomic locus (bubble) contribute ~1 observation rather than m_b.

The rate / mean is UNCHANGED:
    mu_k = lambda * (h^T cn)_k.
Only the M-step accumulation is reweighted (weighted multiplicative EM update):
    h_f <- h_f * ( sum_k omega_k * cn_fk * c_k / mu_k ) / ( sum_k omega_k * cn_fk )
followed by renormalization to the simplex.

This mirrors src/em_solver.py (solve_em) and the densify / nz-filter /
cov-scalar handling in src/archive/sweep_shape_norm_h_only.py, just with
the omega_k weight folded into the E/M accumulation.

Outputs (one .npz per sim) under scratch/h_fixes/perbubble/.
"""
from __future__ import annotations
import argparse, gc, json, os, sys, time
import numpy as np
from scipy.sparse import load_npz


def solve_em_weighted(
    counts: np.ndarray,      # K-vector of observed counts (already nz-filtered)
    cn: np.ndarray,          # F x K binary copy-number (dense float32, nz-filtered)
    omega: np.ndarray,       # K-vector of per-k-mer weights (1/m_b), nz-filtered
    coverage: float,         # lambda scalar (cancels in M-step, kept for parity)
    h_init: np.ndarray | None = None,
    max_iter: int = 200,
    tol: float = 1e-7,
    verbose: bool = False,
) -> tuple[np.ndarray, dict]:
    """Weighted multiplicative Poisson EM.

    mu_k = lambda * (h^T cn)_k  (the cov scalar cancels, so we use the raw h@cn
    just like solve_em; the M-step ratio c_k/mu_k is scale-free in lambda).

    Update:
        cw_k       = omega_k * c_k / (h^T cn)_k             (E-step numerator factor)
        num_f      = h_f * sum_k cn_fk * cw_k                (weighted attributed mass)
        denom_f    = sum_k omega_k * cn_fk                   (weighted carry mass, FIXED)
        h_f        <- num_f / denom_f ; renormalize to simplex
    """
    K = counts.shape[0]
    F = cn.shape[0]
    if cn.dtype != np.float32:
        cn = cn.astype(np.float32)
    counts = counts.astype(np.float32)
    omega = omega.astype(np.float32)

    h = (np.full(F, 1.0 / F, dtype=np.float32) if h_init is None
         else h_init.astype(np.float32))

    # Per-founder weighted carry mass: denom_f = sum_k omega_k * cn_fk  (constant)
    denom_f = np.maximum(cn @ omega, np.float32(1e-12))   # F-vector

    # weighted counts factor: omega_k * c_k  (constant across iterations)
    wcounts = omega * counts                              # K-vector

    history = []
    delta = np.inf
    it = 0
    for it in range(max_iter):
        mu = np.maximum(h @ cn, np.float32(1e-7))         # K-vector  (h^T cn)_k
        cw = wcounts / mu                                 # omega_k * c_k / mu_k
        num_f = h * (cn @ cw)                             # F-vector weighted attributed mass
        h_new = num_f / denom_f                           # weighted multiplicative update
        s = h_new.sum()
        if s <= 0 or not np.isfinite(s):
            # degenerate; bail with current h
            break
        h_new = h_new / s

        delta = float(np.linalg.norm(h_new - h))
        history.append(delta)
        if verbose and it % 20 == 0:
            print(f"    EM iter {it}: ||dh||={delta:.2e} hmin={h_new.min():.4g} "
                  f"hmax={h_new.max():.4g}", flush=True)
        h = h_new
        if delta < tol:
            break

    return h.astype(np.float64), {"iterations": it + 1,
                                  "delta_history": history,
                                  "converged": delta < tol}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cn-prefix", required=True,
                    help="prefix; expects {prefix}.cn.npz and {prefix}.meta.npz")
    ap.add_argument("--counts", required=True,
                    help="precomputed c_k .npy (indexed to filt2 kmer_index)")
    ap.add_argument("--sample", required=True)
    ap.add_argument("--out", required=True, help="output .npz path")
    ap.add_argument("--em-max-iter", type=int, default=300)
    args = ap.parse_args()

    t0 = time.time()
    print(f"=== per-bubble EM: sample={args.sample} ===", flush=True)

    # --- load cn + meta ---
    cn = load_npz(args.cn_prefix + ".cn.npz").tocsr()
    meta = np.load(args.cn_prefix + ".meta.npz", allow_pickle=True)
    founders = np.asarray(meta["founders"]).astype(str)
    bubble_id = np.asarray(meta["bubble_id"]).astype(np.int64)
    F, K = cn.shape
    print(f"  cn: F={F} K={K:,} nnz={cn.nnz:,}  [{time.time()-t0:.0f}s]", flush=True)

    # --- per-bubble multiplicity m_b and weight omega_k = 1/m_b ---
    counts_per_bubble = np.bincount(bubble_id)
    m_b = counts_per_bubble[bubble_id].astype(np.float64)        # K-vector
    omega_full = 1.0 / np.maximum(m_b, 1.0)                       # K-vector
    print(f"  bubbles={counts_per_bubble.size:,}  "
          f"m_b: median={np.median(m_b):.1f} mean={m_b.mean():.2f} "
          f"max={int(m_b.max())}", flush=True)

    # --- a_k allele count (for diagnostics) ---
    a_k = np.asarray(cn.sum(0)).flatten().astype(np.int64)

    # --- load precomputed counts ---
    counts = np.load(args.counts)
    if counts.shape[0] != K:
        raise SystemExit(f"counts length {counts.shape[0]} != K {K}")
    print(f"  counts: nonzero={ (counts>0).sum():,}/{K:,} sum={counts.sum():.0f}",
          flush=True)

    # --- densify cn (float32) ---
    t = time.time()
    cn_dense = np.asarray(cn.todense()).astype(np.float32)
    print(f"  densify cn -> {cn_dense.nbytes/1e9:.1f} GB  [{time.time()-t:.0f}s]",
          flush=True)
    del cn
    gc.collect()

    # --- nz filter (mirror sweep_shape_norm_h_only) ---
    nz = counts > 0
    cn_em = np.ascontiguousarray(cn_dense[:, nz])
    counts_em = counts[nz].astype(np.float32)
    omega_em = omega_full[nz].astype(np.float32)
    cov = counts.sum() * F / max(cn_dense.sum(), 1)
    print(f"  EM matrix: {cn_em.shape} counts_sum={counts_em.sum():.0f} "
          f"cov(lambda)={cov:.4g}", flush=True)
    del cn_dense
    gc.collect()

    # --- run weighted EM ---
    t = time.time()
    h, info = solve_em_weighted(counts_em, cn_em, omega_em, cov,
                                max_iter=args.em_max_iter, tol=1e-7,
                                verbose=True)
    print(f"  EM done: iter={info['iterations']} converged={info['converged']} "
          f"eff_n={1.0/(h**2).sum():.1f}  [{time.time()-t:.0f}s]", flush=True)

    # --- save ---
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    np.savez(args.out,
             sample=args.sample,
             founders=founders,
             h=h,
             cov=cov,
             iterations=info["iterations"],
             converged=info["converged"],
             m_b_median=float(np.median(m_b)),
             m_b_mean=float(m_b.mean()))
    print(f"  wrote {args.out}  [TOTAL {time.time()-t0:.0f}s]", flush=True)


if __name__ == "__main__":
    main()
