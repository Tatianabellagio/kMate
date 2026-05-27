"""
Shape-norm α-sweep — h-only test (no AF projection).

Per-founder weight: w_f = min(1, (target / frac_v_low_ac[f])^α)
where frac_v_low_ac[f] = fraction of founder f's carried k-mers with ac_k <= AC_THRESHOLD.

Apply at EM-load: cn_dense[f, :] *= w_f
Post-EM mass-domain recovery: h_proj[f] = h_em[f] / w_f, then renormalize.

The "intrinsic" property: w_f depends only on cn (the panel), not on side labels
or external coverage data. The clamp at 1 prevents low-evidence founders from
being amplified.

Outputs:
  <out_prefix>.h_sweep.npz    — keys: founders, alphas, ac_thresh,
                                       h_per_alpha (n_alpha × F),
                                       w_per_alpha (n_alpha × F),
                                       frac_v_low_ac (F-vec)
"""
from __future__ import annotations
import argparse, gc, os, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from scipy.sparse import load_npz
from em_solver import solve_em
from kmer_count import count_kmers_in_fasta


def compute_frac_low_ac(cn_csr, ac_threshold: int = 4):
    """Per founder, fraction of carried k-mers with ac_k <= ac_threshold."""
    ac_k = np.asarray(cn_csr.sum(axis=0)).flatten().astype(np.int32)
    rare_mask = (ac_k <= ac_threshold)  # bool K-vec
    F = cn_csr.shape[0]
    frac = np.zeros(F, dtype=np.float64)
    for f in range(F):
        s, e = cn_csr.indptr[f], cn_csr.indptr[f + 1]
        if e <= s: continue
        cols = cn_csr.indices[s:e]
        frac[f] = rare_mask[cols].mean()
    return frac, ac_k


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cn-prefix", required=True,
                    help="prefix; expects {prefix}.cn.npz and {prefix}.meta.npz")
    ap.add_argument("--reads", required=True, nargs="+",
                    help="paired fastqs (1 or 2 files)")
    ap.add_argument("--sample", required=True)
    ap.add_argument("--out-prefix", required=True,
                    help="output prefix; writes {prefix}.h_sweep.npz")
    ap.add_argument("--alphas", default="0,0.3,0.5,1.0,1.5,2.0",
                    help="comma-separated α values to sweep")
    ap.add_argument("--ac-threshold", type=int, default=4,
                    help="k-mers with ac_k <= N are 'very rare' (default 4)")
    ap.add_argument("--target", default="median",
                    choices=["median", "pg_median", "min", "p25"],
                    help="reference target for shape ratio (default panel median)")
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--em-max-iter", type=int, default=200)
    ap.add_argument("--counts-cache", default=None,
                    help="optional: cache c_k to this npy path "
                         "(reuse across α; cheap reruns)")
    args = ap.parse_args()

    alphas = [float(a) for a in args.alphas.split(",")]
    print(f"=== shape-norm α-sweep: sample={args.sample} alphas={alphas} ac<= {args.ac_threshold} ===",
          flush=True)
    t0 = time.time()

    # Load cn + meta
    cn_path = args.cn_prefix + ".cn.npz"
    meta_path = args.cn_prefix + ".meta.npz"
    cn = load_npz(cn_path).tocsr()
    meta = np.load(meta_path, allow_pickle=True)
    founders = np.asarray(meta["founders"]).astype(str)
    kmer_index = list(np.asarray(meta["kmer_index"]).astype(str))
    F, K = cn.shape
    print(f"  cn: F={F} K={K:,} nnz={cn.nnz:,}", flush=True)

    # frac_v_low_ac per founder
    t = time.time()
    frac, ac_k = compute_frac_low_ac(cn, ac_threshold=args.ac_threshold)
    print(f"  frac_v_low_ac (ac<= {args.ac_threshold}): median={np.median(frac):.4f}, "
          f"min={frac.min():.4f}, max={frac.max():.4f}, [{time.time()-t:.0f}s]",
          flush=True)

    # Target
    if args.target == "median":
        target = float(np.median(frac))
    elif args.target == "min":
        target = float(frac.min())
    elif args.target == "p25":
        target = float(np.quantile(frac, 0.25))
    elif args.target == "pg_median":
        # Hack: not intrinsic, exposed for sanity-check parity vs side-aware
        # only useful if --founder-split-json is wired in; left for future
        target = float(np.median(frac))
    print(f"  target = {target:.4f}", flush=True)

    # Count k-mers in reads
    counts = None
    if args.counts_cache and os.path.exists(args.counts_cache):
        counts = np.load(args.counts_cache)
        if counts.shape[0] == K:
            print(f"  counts cache hit: {args.counts_cache} ({K:,} entries)", flush=True)
        else:
            print(f"  counts cache size mismatch ({counts.shape[0]} vs {K}); recounting", flush=True)
            counts = None
    if counts is None:
        t = time.time()
        cd = count_kmers_in_fasta(args.reads, kmer_index, k=31,
                                   threads=args.threads, hash_size="3G")
        counts = np.array([cd[km] for km in kmer_index], dtype=np.int64)
        print(f"  jellyfish: {(counts > 0).sum():,}/{K:,} nonzero  [{time.time()-t:.0f}s]",
              flush=True)
        if args.counts_cache:
            np.save(args.counts_cache, counts)
            print(f"  cached counts → {args.counts_cache}", flush=True)

    # Densify cn once (float32)
    t = time.time()
    cn_dense_base = np.asarray(cn.todense()).astype(np.float32)
    print(f"  densify cn → {cn_dense_base.nbytes/1e9:.1f} GB float32  [{time.time()-t:.0f}s]",
          flush=True)
    del cn
    gc.collect()

    # Filter to nz counts (smaller EM matrix)
    nz = counts > 0
    cn_em_base = np.ascontiguousarray(cn_dense_base[:, nz])
    counts_em = counts[nz].astype(np.float32)
    print(f"  EM matrix: {cn_em_base.shape}, counts_sum={counts_em.sum():.0f}", flush=True)
    cov = counts.sum() * F / max(cn_dense_base.sum(), 1)
    del cn_dense_base
    gc.collect()

    h_per_alpha = np.zeros((len(alphas), F), dtype=np.float64)
    w_per_alpha = np.zeros((len(alphas), F), dtype=np.float64)
    h_proj_per_alpha = np.zeros((len(alphas), F), dtype=np.float64)
    info_per_alpha = []

    for ia, a in enumerate(alphas):
        # Compute w_f with the clamp
        if a == 0.0:
            w = np.ones(F, dtype=np.float32)
        else:
            ratio = target / np.maximum(frac, 1e-9)  # >1 for rare-poor, <1 for rare-rich
            w = np.minimum(1.0, ratio ** a).astype(np.float32)
        w_per_alpha[ia] = w

        # Apply weight to cn rows (scale per-row)
        cn_em = cn_em_base * w[:, None]

        # Run EM
        t = time.time()
        h, info = solve_em(counts_em, cn_em, cov,
                            max_iter=args.em_max_iter, tol=1e-7)
        # Mass-domain recovery: h_proj = h_em / w_f, renormalize
        h_proj = h / np.maximum(w.astype(np.float64), 1e-12)
        h_proj = h_proj / h_proj.sum()
        h_per_alpha[ia] = h
        h_proj_per_alpha[ia] = h_proj
        info_per_alpha.append({
            "alpha": float(a),
            "iterations": int(info["iterations"]),
            "w_min": float(w.min()), "w_max": float(w.max()), "w_med": float(np.median(w)),
            "h_em_eff_n": float(1.0 / (h ** 2).sum()),
            "h_proj_eff_n": float(1.0 / (h_proj ** 2).sum()),
        })
        print(f"  α={a:>4.1f}: w∈[{w.min():.3f},{w.max():.3f}] med={np.median(w):.3f}, "
              f"EM iter={info['iterations']}, "
              f"eff_n (EM)={1/(h**2).sum():.1f}, eff_n (proj)={1/(h_proj**2).sum():.1f} "
              f"[{time.time()-t:.0f}s]", flush=True)

    # Save
    out_path = args.out_prefix + ".h_sweep.npz"
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    np.savez(out_path,
             founders=founders,
             alphas=np.array(alphas),
             ac_threshold=args.ac_threshold,
             target=target,
             frac_v_low_ac=frac,
             h_per_alpha=h_per_alpha,
             h_proj_per_alpha=h_proj_per_alpha,
             w_per_alpha=w_per_alpha,
             info=np.array(info_per_alpha, dtype=object),
             )
    print(f"\nWrote {out_path}", flush=True)
    print(f"  TOTAL: {time.time()-t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
