"""
INVERSE-AC LIKELIHOOD WEIGHTING for the Poisson k-mer founder-frequency EM.

PROBLEM
-------
h (231-founder simplex) is estimated by  c_k ~ Poisson(mu_k),  mu_k = lambda * (h^T cn)_k.
The filt2 cn_full over-credits the 80 "cactus" founders because, in the Fisher
information  G_fg = sum_k cn_fk cn_gk / mu_k,  the low-AC (rare) k-mers dominate
and cactus founders carry 2.4x more of them. The simplex constraint then rectifies
the resulting variance into an upward class bias.

METHOD
------
Weight each k-mer's LIKELIHOOD TERM by  omega_k = 1 / a_k^gamma,  where
a_k = number of founders carrying k-mer k (a_k >= 2 in filt2, since ac=1 dropped).

The weight multiplies the likelihood term, NOT the cn design rows:
    mu_k = lambda * (h^T cn)_k        <-- UNCHANGED
There is therefore NO h/w "mass-recovery" step.

Weighted EM update (multiplicative, then renormalize to simplex):
    h_f  <-  h_f * ( sum_k omega_k * cn_fk * c_k / mu_k )
                   / ( sum_k omega_k * cn_fk )

This is the multiplicative-EM solution of the weighted estimating equation
    sum_k omega_k * cn_fk * ( c_k / mu_k - lambda ) = 0   (per founder f)
which is mean-zero at truth for any fixed omega_k (weighted score), hence
consistent. gamma = 0 recovers the standard (unweighted) EM.

Mirrors src/em_solver.solve_em and the archived shape-norm sweep for the
exact lambda / coverage scalar, nz=(counts>0) filtering, densify, simplex handling.
Adds only omega_k.
"""
from __future__ import annotations
import argparse, gc, json, os, sys, time
import numpy as np
from scipy.sparse import load_npz


def weighted_em(counts_em, cn_em, omega, cov, max_iter=300, tol=1e-7, verbose=False):
    """Weighted multiplicative Poisson-EM on the nz-filtered design.

    counts_em : (Knz,) float32   observed k-mer counts (counts>0 subset)
    cn_em     : (F, Knz) float32  dense binary copy-number (UNWEIGHTED rows)
    omega     : (Knz,) float32   per-k-mer likelihood weight 1/a_k^gamma
    cov       : float            lambda = counts.sum()*F / cn.sum()  (full-K scalar)

    mu_k = cov * (h @ cn_em)_k.  Returns (h float64, info dict).
    """
    F = cn_em.shape[0]
    h = np.full(F, 1.0 / F, dtype=np.float32)

    cov = np.float32(cov)
    # Denominator of the update is constant across iterations:
    #   D_f = sum_k omega_k * cn_fk
    # (no mu, no counts).  cn_em is binary so this is the omega-weighted
    # number of carried nz-k-mers for founder f.
    D = (cn_em @ omega).astype(np.float32)            # (F,)
    D = np.maximum(D, np.float32(1e-12))

    omega_c = (omega * counts_em).astype(np.float32)  # (Knz,) = omega_k * c_k

    history = []
    it = 0
    for it in range(max_iter):
        mu = np.maximum(cov * (h @ cn_em), np.float32(1e-7))   # (Knz,)
        # numerator N_f = h_f * sum_k omega_k * cn_fk * c_k / mu_k
        N = h * (cn_em @ (omega_c / mu))                       # (F,)
        h_new = N / D
        s = h_new.sum()
        if s <= 0:
            break
        h_new = h_new / s
        delta = float(np.linalg.norm(h_new - h))
        history.append(delta)
        if verbose and it % 20 == 0:
            print(f"    iter {it}: ||dh||={delta:.2e} h_min={h_new.min():.5f} "
                  f"h_max={h_new.max():.5f}", flush=True)
        h = h_new
        if delta < tol:
            break
    return h.astype(np.float64), {"iterations": it + 1,
                                  "converged": history and history[-1] < tol,
                                  "delta_last": history[-1] if history else None}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cn-prefix", required=True,
                    help="prefix; expects {prefix}.cn.npz and {prefix}.meta.npz")
    ap.add_argument("--counts", required=True, nargs="+",
                    help="one or more *.counts.npy (indexed to filt2 kmer_index)")
    ap.add_argument("--labels", required=True, nargs="+",
                    help="label per counts file (same length/order)")
    ap.add_argument("--gammas", default="0.0,0.25,0.5,0.75,1.0",
                    help="comma-separated gamma values to sweep (0=unweighted)")
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--em-max-iter", type=int, default=300)
    args = ap.parse_args()

    assert len(args.counts) == len(args.labels), "counts/labels length mismatch"
    gammas = [float(g) for g in args.gammas.split(",")]
    print(f"=== invac sweep: gammas={gammas} samples={args.labels} ===", flush=True)
    t0 = time.time()

    # --- load cn + meta ---
    cn = load_npz(args.cn_prefix + ".cn.npz").tocsr()
    meta = np.load(args.cn_prefix + ".meta.npz", allow_pickle=True)
    founders = np.asarray(meta["founders"]).astype(str)
    F, K = cn.shape
    print(f"  cn: F={F} K={K:,} nnz={cn.nnz:,}", flush=True)

    # a_k = number of founders carrying k-mer k (full-K)
    a_k = np.asarray(cn.sum(axis=0)).flatten().astype(np.float64)   # (K,)
    print(f"  a_k: min={a_k.min():.0f} median={np.median(a_k):.0f} "
          f"max={a_k.max():.0f}  (frac a_k==2: {(a_k==2).mean():.3f})", flush=True)

    # densify cn once (float32)
    t = time.time()
    cn_dense = np.asarray(cn.todense()).astype(np.float32)
    print(f"  densify cn -> {cn_dense.nbytes/1e9:.1f} GB float32  [{time.time()-t:.0f}s]",
          flush=True)
    cn_full_sum = float(cn.sum())
    del cn
    gc.collect()

    results = {}  # label -> {gamma -> h}
    for cpath, label in zip(args.counts, args.labels):
        counts = np.load(cpath).astype(np.float64)
        assert counts.shape[0] == K, f"{label}: counts K={counts.shape[0]} != cn K={K}"
        nz = counts > 0
        cn_em = np.ascontiguousarray(cn_dense[:, nz])
        counts_em = counts[nz].astype(np.float32)
        a_nz = a_k[nz]
        cov = counts.sum() * F / max(cn_full_sum, 1.0)
        print(f"\n[{label}] nz={nz.sum():,}/{K:,} counts_sum={counts.sum():.0f} "
              f"cov(lambda)={cov:.4f}", flush=True)

        results[label] = {}
        for g in gammas:
            omega = (1.0 / np.power(a_nz, g)).astype(np.float32)   # gamma=0 -> all ones
            t = time.time()
            h, info = weighted_em(counts_em, cn_em, omega, cov,
                                  max_iter=args.em_max_iter, tol=1e-7)
            eff_n = 1.0 / (h ** 2).sum()
            print(f"  gamma={g:>4.2f}: iter={info['iterations']:>3d} "
                  f"conv={info['converged']} eff_n={eff_n:.1f} "
                  f"h_max={h.max():.4f}  [{time.time()-t:.0f}s]", flush=True)
            results[label][f"{g:g}"] = h

        del cn_em, counts_em
        gc.collect()

    # --- save: one npz per label, keys = "h_gamma<g>" ---
    os.makedirs(os.path.dirname(args.out_prefix) or ".", exist_ok=True)
    for label in results:
        out = {"founders": founders, "gammas": np.array(gammas)}
        for gk, h in results[label].items():
            out[f"h_gamma{gk}"] = h
        outp = f"{args.out_prefix}_{label}.invac.npz"
        np.savez(outp, **out)
        print(f"  wrote {outp}", flush=True)

    print(f"\nTOTAL: {time.time()-t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
