"""
BIAS-FREE ANCHOR via BALANCED-BUBBLE SUB-DESIGN for the Poisson k-mer
founder-frequency EM.

PROBLEM
-------
h (231-founder simplex) is estimated by  c_k ~ Poisson(mu_k),  mu_k = lambda * (h^T cn)_k.
The filt2 cn_full over-credits the 80 "cactus" founders (real-SEEDMIX cactus mass
0.518 vs true ~0.35) because the DESIGN is IMBALANCED: at some bubbles (local
variant clusters given by meta bubble_id) cactus founders carry systematically
MORE k-mers than PG founders (long-read assemblies expose extra structure), so
those loci over-weight cactus. hapFIRE / HARP avoid this because their per-SNP
design is balanced by construction.

METHOD = ESTIMATE h ON A BALANCED-BUBBLE SUB-DESIGN, then STANDARD UNWEIGHTED EM
--------------------------------------------------------------------------------
The selection criterion is BALANCE, not variant type. For each bubble b:

    n_{f,b} = number of b's k-mers carried by founder f               (231 x B)

Split founders into the cactus / PG classes (data/founder_split_cactus_pg.json).
Per-founder class summaries:

    med_cact_b = median_{f in cactus} n_{f,b}
    med_pg_b   = median_{f in PG}     n_{f,b}

Keep bubble b iff the two classes are balanced:

    ratio_b = med_cact_b / med_pg_b  in  [1/T, T]

SWEEP T in {1.25, 1.5, 2.0}  (tighter T = stricter balance, fewer bubbles).

A SNP bubble naturally passes (every founder carries one allele's k-mers, so the
class medians are equal); but indels/SVs ALSO pass if both classes are uniformly
represented -- that is the point: BALANCE, not SNP-ness.

The kept bubbles select a COLUMN MASK over the full filt2 kmer_index. We then run
the STANDARD UNWEIGHTED EM (mirror src/em_solver.solve_em and the
densify / nz-filter / cov-scalar handling of
src/archive/sweep_shape_norm_h_only.py) on the kept-column subset, using
the EXISTING filt2 counts boolean-masked to the kept columns. NO re-counting, NO
jellyfish; the counts are already aligned to filt2's full kmer_index.

This produces a provably-less-biased REFERENCE h: no per-founder weighting, no
deletion of rich founders' rows -- only loci where the two classes are balanced
contribute, so the design itself can no longer over-weight cactus.

Outputs (one .npz per sim) under scratch/h_fixes/balancedbubble/.
  <out>_<label>.balbub.npz keys:
      founders, T_list, h_T<T> (per-T h-vector), n_kept_bubbles_T<T>,
      n_kept_kmers_T<T>, kept_cact_med_T<T>, kept_pg_med_T<T> (diagnostics)
Plus once, the design-level subset stats are saved to <out>_subset_stats.json.
"""
from __future__ import annotations
import argparse, gc, json, os, time
import numpy as np
from scipy.sparse import load_npz, csc_matrix


# ----------------------------------------------------------------------
# Standard UNWEIGHTED multiplicative Poisson EM (mirror em_solver.solve_em).
# cov (lambda) cancels in the M-step, kept only for parity/diagnostics.
# ----------------------------------------------------------------------
def solve_em(counts, cn, h_init=None, max_iter=300, tol=1e-7, verbose=False):
    F = cn.shape[0]
    if cn.dtype != np.float32:
        cn = cn.astype(np.float32)
    counts = counts.astype(np.float32)
    h = (np.full(F, 1.0 / F, dtype=np.float32) if h_init is None
         else h_init.astype(np.float32))
    total_c = counts.sum()
    history = []
    it = 0
    for it in range(max_iter):
        denom = np.maximum(h @ cn, np.float32(1e-7))     # (h^T cn)_k
        cw = counts / denom
        em_term = h * (cn @ cw)
        h_new = em_term / max(total_c, np.float32(1e-12))
        h_new = h_new / h_new.sum()
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


def per_founder_per_bubble_counts(cn_csr, bubble_id, n_bubbles):
    """n_{f,b} = number of b's k-mers carried by founder f.  Returns (F, B) int64.

    For each founder row f of the CSR, bincount the bubble_id of its carried
    columns. cn is binary (filt2), so a carried column contributes exactly 1.
    """
    F = cn_csr.shape[0]
    N = np.zeros((F, n_bubbles), dtype=np.int64)
    indptr, indices = cn_csr.indptr, cn_csr.indices
    for f in range(F):
        s, e = indptr[f], indptr[f + 1]
        if e <= s:
            continue
        cols = indices[s:e]
        N[f] = np.bincount(bubble_id[cols], minlength=n_bubbles)
    return N


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cn-prefix", required=True,
                    help="prefix; expects {prefix}.cn.npz and {prefix}.meta.npz")
    ap.add_argument("--counts", required=True, nargs="+",
                    help="one or more *.counts.npy (indexed to filt2 kmer_index)")
    ap.add_argument("--labels", required=True, nargs="+",
                    help="label per counts file (same length/order)")
    ap.add_argument("--split-json", required=True,
                    help="founder_split_cactus_pg.json")
    ap.add_argument("--tols", default="1.25,1.5,2.0",
                    help="comma-separated balance tolerances T")
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--em-max-iter", type=int, default=300)
    args = ap.parse_args()

    assert len(args.counts) == len(args.labels), "counts/labels length mismatch"
    T_list = [float(t) for t in args.tols.split(",")]
    print(f"=== balanced-bubble sub-design EM: T={T_list} samples={args.labels} ===",
          flush=True)
    t0 = time.time()

    # --- load cn + meta ---
    cn = load_npz(args.cn_prefix + ".cn.npz").tocsr()
    meta = np.load(args.cn_prefix + ".meta.npz", allow_pickle=True)
    founders = np.asarray(meta["founders"]).astype(str)
    bubble_id = np.asarray(meta["bubble_id"]).astype(np.int64)
    F, K = cn.shape
    n_bubbles = int(bubble_id.max()) + 1
    print(f"  cn: F={F} K={K:,} nnz={cn.nnz:,}  bubbles={n_bubbles:,}  "
          f"[{time.time()-t0:.0f}s]", flush=True)

    # --- class indices ---
    split = json.load(open(args.split_json))
    cact_ids = set(map(str, split["cactus"]))
    pg_ids = set(map(str, split["PG"]))
    cact_idx = np.array([i for i, f in enumerate(founders) if f in cact_ids])
    pg_idx = np.array([i for i, f in enumerate(founders) if f in pg_ids])
    print(f"  classes: cactus={cact_idx.size} PG={pg_idx.size}", flush=True)

    # --- per-founder per-bubble counts n_{f,b} ---
    t = time.time()
    N = per_founder_per_bubble_counts(cn, bubble_id, n_bubbles)   # (F, B) int64
    print(f"  n_{{f,b}}: {N.shape} built  [{time.time()-t:.0f}s]", flush=True)

    # k-mers per bubble (size m_b) for reporting
    m_b = np.bincount(bubble_id, minlength=n_bubbles)            # (B,)

    # --- per-class per-bubble median over founders (zeros included) ---
    # med over ALL founders in the class: a founder carrying none of b's k-mers
    # contributes n=0, so the median reflects "is this bubble uniformly carried
    # across the class".
    med_cact = np.median(N[cact_idx, :], axis=0)                 # (B,)
    med_pg = np.median(N[pg_idx, :], axis=0)                     # (B,)
    # also mean for diagnostics
    mean_cact = N[cact_idx, :].mean(axis=0)
    mean_pg = N[pg_idx, :].mean(axis=0)

    # balance ratio med_cact / med_pg ; guard zero denominators.
    # A bubble where one class has median 0 is by definition NOT balanced
    # (one class is not uniformly represented) -> ratio set to inf/0 -> rejected
    # unless BOTH are 0 (empty/degenerate bubble) -> exclude.
    both_zero = (med_cact == 0) & (med_pg == 0)
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = med_cact / med_pg                                # (B,)
    # bubbles with med_pg==0 (but med_cact>0) -> ratio=inf -> rejected
    # bubbles with med_cact==0 (but med_pg>0) -> ratio=0   -> rejected

    print(f"  per-bubble class medians: med_cact(med over B)="
          f"{np.median(med_cact):.2f} med_pg={np.median(med_pg):.2f}", flush=True)
    print(f"  both-zero bubbles (excluded): {int(both_zero.sum()):,}", flush=True)

    # --- densify cn once (float32) for the EM ---
    t = time.time()
    cn_dense = np.asarray(cn.todense()).astype(np.float32)       # (F, K)
    print(f"  densify cn -> {cn_dense.nbytes/1e9:.1f} GB float32  "
          f"[{time.time()-t:.0f}s]", flush=True)
    del cn
    gc.collect()

    # ------------------------------------------------------------------
    # Build per-T column masks + run unweighted EM per (T, sample)
    # ------------------------------------------------------------------
    subset_stats = {}
    # results[label][T] = dict(h=..., n_kept_bubbles, n_kept_kmers, ...)
    results = {lab: {} for lab in args.labels}

    # pre-load counts
    counts_by_label = {}
    for cpath, label in zip(args.counts, args.labels):
        c = np.load(cpath).astype(np.float64)
        assert c.shape[0] == K, f"{label}: counts K={c.shape[0]} != cn K={K}"
        counts_by_label[label] = c

    for T in T_list:
        Tk = f"{T:g}"
        keep_bub = (~both_zero) & (ratio >= 1.0 / T) & (ratio <= T)   # (B,) bool
        n_kept_bub = int(keep_bub.sum())
        # column (k-mer) mask: a k-mer is kept iff its bubble is kept
        col_keep = keep_bub[bubble_id]                               # (K,) bool
        n_kept_kmers = int(col_keep.sum())

        # design-level balance of the KEPT sub-design:
        # per-founder total carried k-mers within kept columns, by class
        kept_per_founder = cn_dense[:, col_keep].sum(axis=1)        # (F,)
        cact_pf = kept_per_founder[cact_idx]
        pg_pf = kept_per_founder[pg_idx]
        bal_med = float(np.median(cact_pf) / max(np.median(pg_pf), 1e-9))
        bal_mean = float(cact_pf.mean() / max(pg_pf.mean(), 1e-9))

        subset_stats[Tk] = dict(
            T=T,
            n_kept_bubbles=n_kept_bub,
            frac_kept_bubbles=n_kept_bub / float((~both_zero).sum()),
            n_kept_kmers=n_kept_kmers,
            frac_kept_kmers=n_kept_kmers / float(K),
            # per-founder k-mer balance cactus/PG in kept sub-design:
            kept_perfounder_cact_med=float(np.median(cact_pf)),
            kept_perfounder_pg_med=float(np.median(pg_pf)),
            kept_perfounder_balance_median=bal_med,
            kept_perfounder_balance_mean=bal_mean,
            # full-design balance for reference (T=inf):
        )
        print(f"\n[T={Tk}] kept bubbles={n_kept_bub:,}/"
              f"{int((~both_zero).sum()):,} "
              f"({100*subset_stats[Tk]['frac_kept_bubbles']:.1f}%)  "
              f"kept k-mers={n_kept_kmers:,}/{K:,} "
              f"({100*subset_stats[Tk]['frac_kept_kmers']:.1f}%)", flush=True)
        print(f"       per-founder k-mer balance cactus/PG (kept): "
              f"median={bal_med:.3f} mean={bal_mean:.3f}", flush=True)

        # densify the kept sub-design ONCE per T (shared across samples)
        cn_keep = np.ascontiguousarray(cn_dense[:, col_keep])      # (F, n_kept_kmers)
        cn_keep_sum = float(cn_keep.sum())

        for label in args.labels:
            counts = counts_by_label[label]
            c_keep = counts[col_keep]                              # masked to kept cols
            # nz filter within the kept sub-design (mirror sweep_shape_norm)
            nz = c_keep > 0
            cn_em = np.ascontiguousarray(cn_keep[:, nz])
            counts_em = c_keep[nz].astype(np.float32)
            cov = c_keep.sum() * F / max(cn_keep_sum, 1.0)
            t = time.time()
            h, info = solve_em(counts_em, cn_em, max_iter=args.em_max_iter,
                               tol=1e-7)
            cm = float(h[cact_idx].sum())
            eff_n = 1.0 / (h ** 2).sum()
            print(f"    [{label}] nz={int(nz.sum()):,} counts_sum="
                  f"{counts_em.sum():.0f} cov={cov:.3g} -> iter={info['iterations']} "
                  f"cact_mass={cm:.4f} eff_n={eff_n:.1f}  [{time.time()-t:.0f}s]",
                  flush=True)
            results[label][Tk] = dict(
                h=h, n_kept_bubbles=n_kept_bub, n_kept_kmers=n_kept_kmers,
                kept_perfounder_balance_median=bal_med,
                kept_perfounder_balance_mean=bal_mean, cact_mass=cm,
                eff_n=float(eff_n), iterations=int(info["iterations"]))
            del cn_em
            gc.collect()
        del cn_keep
        gc.collect()

    # --- save per-label npz ---
    os.makedirs(os.path.dirname(args.out_prefix) or ".", exist_ok=True)
    for label in args.labels:
        out = {"founders": founders, "T_list": np.array(T_list)}
        for Tk, r in results[label].items():
            out[f"h_T{Tk}"] = r["h"]
            out[f"n_kept_bubbles_T{Tk}"] = r["n_kept_bubbles"]
            out[f"n_kept_kmers_T{Tk}"] = r["n_kept_kmers"]
            out[f"balance_median_T{Tk}"] = r["kept_perfounder_balance_median"]
            out[f"cact_mass_T{Tk}"] = r["cact_mass"]
        outp = f"{args.out_prefix}_{label}.balbub.npz"
        np.savez(outp, **out)
        print(f"  wrote {outp}", flush=True)

    # --- save subset-level stats once ---
    with open(f"{args.out_prefix}_subset_stats.json", "w") as f:
        json.dump(subset_stats, f, indent=2)
    print(f"  wrote {args.out_prefix}_subset_stats.json", flush=True)
    print(f"\nTOTAL: {time.time()-t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
