"""
IRLS / GLS whitening for k-mer founder-frequency deconvolution.
=================================================================

PROBLEM
-------
Estimate h (231-founder simplex) from k-mer counts under the generative model
    c_k ~ Poisson( lambda * (h^T cn)_k ),     cn in {0,1}^{F x K}
The filt2 cn_full over-credits the 80 "cactus" founders because the design is
imbalanced/correlated: many cactus founders share long, highly-multi-k-mer
bubbles, so plain EM (= unweighted Poisson MLE) lets those founders soak up
mass.  The classic deconvolution remedy (MuSiC = inverse-variance weighted NNLS;
CIBERSORT = nu-SVR on whitened features) is to whiten / down-weight the
correlated, high-variance observations.

METHOD = IRLS / GLS WHITENING
-----------------------------
Treat the deconvolution as weighted least squares (a Gaussian / IRLS
approximation to the Poisson GLM):

    h_hat = argmin_{h>=0, sum h = 1}  (c - mu)^T  Sigma^{-1}  (c - mu),
            mu = lambda * cn^T h

We use a TRACTABLE block-diagonal Sigma^{-1} that unifies the two corrections
the field uses:

  (1) Poisson inverse-variance:  Var(c_k) ~ mu_k  ->  weight 1/mu_k
  (2) within-bubble correlation deflation: a bubble of m_b k-mers carries
      ~one independent genomic observation, not m_b.  So divide each k-mer's
      weight by its bubble size m_b.

  =>  effective per-k-mer GLS weight   omega_k = 1 / ( m_b(k) * mu_k )

This is exactly the GLS that combines per-bubble normalization + Poisson
variance, with a diagonal-times-block-scalar approximation to Sigma^{-1}
(we approximate the within-bubble correlation block by 1/m_b * identity,
i.e. (1/m_b) * I rather than the full equicorrelated inverse -- the cheap,
positive, tractable surrogate).

IRLS LOOP
---------
Iterate to convergence:
  (i)   mu_k   = lambda * (cn^T h)_k                       [current fit]
  (ii)  omega_k = 1 / ( m_b(k) * max(mu_k, eps) )          [GLS weights]
  (iii) re-solve weighted NNLS for h on the simplex, then renormalize:
            h <- argmin_{h>=0} sum_k omega_k ( c_k - lambda*(cn^T h)_k )^2
        solved by a multiplicative (MU) / projected update that keeps h>=0
        and is monotone for weighted-NLS with a nonnegative design.
  (iv)  update lambda = sum_k omega_k c_k (cn^T h)_k / sum_k omega_k (cn^T h)_k^2
        (the closed-form weighted-LS scale given the current shape cn^T h)

The weights multiply the residual / likelihood, NOT the design rows: mu stays
= lambda * cn^T h, so there is no h/w recovery step (unlike the shape-norm
sweep, which scaled cn rows and had to divide h back out).

INNER WEIGHTED-NNLS SOLVER
--------------------------
We avoid materializing the K x F dense weighted design (K ~ 3.7M, F = 231 ->
3.4 GB and a dense lstsq is wasteful).  Instead we use a multiplicative update
for nonnegative weighted least squares (NNLS), which only needs sparse
matvecs cn @ v and cn^T @ v:

    let  y = c  (targets),  A = lambda * cn^T  (K x F design, applied via cn)
    define  a_f = (A^T W y)_f       = lambda * ( cn @ (omega * c) )_f
            b_f = (A^T W A h)_f      = lambda^2 * ( cn @ (omega * (cn^T h)) )_f
    multiplicative step (nonneg-preserving):
            h_f <- h_f * a_f / max(b_f, eps)
    then project to the simplex (renormalize) for the mass constraint.
This is the standard NMF/ NNLS multiplicative update (Lee & Seung) specialized
to a single mixing vector; it monotonically decreases the weighted SSE and
keeps h >= 0.  A few inner MU steps per IRLS outer step suffice.

CONVERGENCE
-----------
Outer IRLS stops when ||h_new - h|| < tol (default 1e-7) or max_iter reached.

INPUTS (reuse, NO jellyfish)
----------------------------
  cn      data/cn_full_231_v3qc_v3_filt2/cn_Chr1.cn.npz  (231 x ~11.15M CSR)
  meta    .../cn_Chr1.meta.npz   (kmer_index, bubble_id, founders)
  counts  scratch/g0_sweep_h_test/filt2_<SIM>.counts.npy
          scratch/seedmix_h_test/filt2_S{1..8}.counts.npy   (filt2-indexed, len == K)

We work on the nz subset (counts > 0) to keep the matrix manageable; the
bubble-size m_b is computed on the FULL panel (np.bincount(bubble_id)[bubble_id])
before subsetting, so deflation reflects true bubble sizes.
"""
from __future__ import annotations
import argparse, gc, os, sys, time
import numpy as np
from scipy.sparse import load_npz, csc_matrix


def whitened_irls(
    counts_nz: np.ndarray,     # K' observed counts (>0 subset), float64
    cnT_nz: csc_matrix,        # K' x F  (counts-major) copy-number, csr/csc
    cn_nz,                     # F x K'  copy-number (csr) for cn @ v matvecs
    mb_nz: np.ndarray,         # K' bubble sizes (deflation denom)
    F: int,
    max_iter: int = 200,
    inner_iter: int = 5,
    tol: float = 1e-7,
    eps: float = 1e-9,
    verbose: bool = True,
):
    """IRLS/GLS whitening solve. Returns (h, info)."""
    y = counts_nz.astype(np.float64)
    total_c = y.sum()
    # init: uniform simplex
    h = np.full(F, 1.0 / F, dtype=np.float64)

    # lambda init from a plain (unweighted) scale: total counts / total expected
    # shape with uniform h. mu0 = (cn^T h); lambda0 = sum(c)/sum(mu0)
    shape0 = cn_nz.T @ h            # K'
    lam = total_c / max(shape0.sum(), eps)

    hist = []
    for it in range(max_iter):
        h_prev = h.copy()
        # (i) current fit shape and mu
        shape = cn_nz.T @ h                      # (cn^T h)_k , K'
        mu = lam * shape
        # (ii) GLS weights: omega_k = 1/(m_b * mu_k)
        omega = 1.0 / (mb_nz * np.maximum(mu, eps))

        # (iii) inner weighted-NNLS multiplicative updates for h (h>=0),
        #       then simplex projection. Design A = lam * cn^T.
        #   a_f = (A^T W y)_f      = lam * cn @ (omega * y)
        #   b_f = (A^T W A h)_f    = lam^2 * cn @ (omega * (cn^T h))
        wy = omega * y
        a = lam * (cn_nz @ wy)                   # F-vec, fixed within inner loop
        for _ in range(inner_iter):
            shp = cn_nz.T @ h                    # K'
            b = (lam * lam) * (cn_nz @ (omega * shp))   # F-vec
            h = h * a / np.maximum(b, eps)
            s = h.sum()
            if s > 0:
                h = h / s                        # simplex projection (mass=1)

        # (iv) update lambda: closed-form weighted-LS scale given shape cn^T h
        shape = cn_nz.T @ h
        num = np.sum(omega * y * shape)
        den = np.sum(omega * shape * shape)
        lam = num / max(den, eps)

        delta = float(np.linalg.norm(h - h_prev))
        hist.append(delta)
        if verbose and (it % 10 == 0 or delta < tol):
            cm_dbg = h.max()
            print(f"  IRLS it={it:3d} ||dh||={delta:.3e} lam={lam:.4g} "
                  f"h_max={cm_dbg:.4f} eff_n={1.0/(h**2).sum():.1f}", flush=True)
        if delta < tol:
            break

    info = {"iterations": it + 1, "converged": delta < tol,
            "lambda": float(lam), "delta_history": hist}
    return h, info


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cn-prefix", required=True,
                    help="expects {prefix}.cn.npz and {prefix}.meta.npz")
    ap.add_argument("--counts", required=True, help="filt2-indexed counts .npy (len==K)")
    ap.add_argument("--sample", required=True)
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--max-iter", type=int, default=200)
    ap.add_argument("--inner-iter", type=int, default=5)
    ap.add_argument("--tol", type=float, default=1e-7)
    args = ap.parse_args()

    t0 = time.time()
    print(f"=== whitening IRLS/GLS: sample={args.sample} ===", flush=True)

    cn = load_npz(args.cn_prefix + ".cn.npz").tocsr()       # F x K
    meta = np.load(args.cn_prefix + ".meta.npz", allow_pickle=True)
    founders = np.asarray(meta["founders"]).astype(str)
    bubble_id = np.asarray(meta["bubble_id"]).astype(np.int64)
    F, K = cn.shape
    print(f"  cn F={F} K={K:,} nnz={cn.nnz:,}  [{time.time()-t0:.0f}s]", flush=True)

    counts = np.load(args.counts)
    assert counts.shape[0] == K, f"counts len {counts.shape[0]} != K {K}"

    # bubble sizes on FULL panel, then subset
    mb = np.bincount(bubble_id, minlength=bubble_id.max() + 1)[bubble_id].astype(np.float64)

    nz = counts > 0
    print(f"  nz counts: {nz.sum():,}/{K:,}  sum={counts.sum():,}", flush=True)
    cn_nz = cn[:, nz].tocsr()                  # F x K'
    cn_nz.sort_indices()
    cnT_nz = cn_nz.T.tocsc()                    # K' x F (not strictly needed; matvecs via cn_nz)
    counts_nz = counts[nz].astype(np.float64)
    mb_nz = mb[nz]
    del cn
    gc.collect()
    print(f"  EM matrix F x K' = {cn_nz.shape}, nnz={cn_nz.nnz:,}  [{time.time()-t0:.0f}s]",
          flush=True)

    h, info = whitened_irls(
        counts_nz, cnT_nz, cn_nz, mb_nz, F,
        max_iter=args.max_iter, inner_iter=args.inner_iter, tol=args.tol,
    )
    print(f"  done: iters={info['iterations']} converged={info['converged']} "
          f"lambda={info['lambda']:.4g}  [{time.time()-t0:.0f}s]", flush=True)

    out = args.out_prefix + ".whiten.npz"
    os.makedirs(os.path.dirname(out), exist_ok=True)
    np.savez(out,
             founders=founders,
             h=h.astype(np.float64),
             sample=args.sample,
             iterations=info["iterations"],
             converged=info["converged"],
             lam=info["lambda"],
             delta_history=np.array(info["delta_history"]))
    print(f"  wrote {out}  [TOTAL {time.time()-t0:.0f}s]", flush=True)


if __name__ == "__main__":
    main()
