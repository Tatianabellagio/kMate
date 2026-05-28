"""
BIAS-FREE ANCHOR via a BALANCED SUB-DESIGN for the Poisson k-mer founder-freq EM.

PROBLEM
-------
We estimate h (231-founder simplex) under c_k ~ Poisson(lambda * (h^T cn)_k).
The filt2 cn_full OVER-CREDITS the 80 "cactus" founders (real-SEEDMIX cactus
mass ~0.518 vs true ~0.35) because the DESIGN is IMBALANCED: cactus founders
carry richer / rarer k-mer sets, so they accrue more attributable count-mass
per unit true frequency. hapFIRE / HARP avoid this because their SNP design is
BALANCED BY CONSTRUCTION: at each biallelic SNP every founder carries exactly one
allele, hence exactly the k-mer set of that allele -> equal per-founder support.

METHOD = ESTIMATE h ON A BALANCED SUB-DESIGN (the hapFIRE analogue), then run
the STANDARD UNWEIGHTED EM on it. No weighting, no deletion of rich founders.

BALANCED SUBSET CONSTRUCTION
----------------------------
Restrict to k-mers whose bubble is a SINGLE-REFERENCE-BASE bubble, i.e. a
biallelic-SNP locus. The coordinate convention in the meta is:
    span_b = bubble_end[b] - bubble_start[b]
and span_b == 1 means the bubble spans exactly ONE reference base (a SNP).
At such a locus each founder carries exactly the k-mers of its single allele
(ref or alt), so per-founder k-mer support is symmetric across the cactus/PG
classes by construction -- this is the direct analogue of hapFIRE's one-allele-
per-founder-per-SNP design.

We keep every k-mer k with span[bubble_id[k]] == SPAN (default SPAN=1).
The kept-column mask indexes directly into filt2's full kmer_index, so the
precomputed counts (.npy, length K) are subset by the same mask -- NO recount.

EM
--
The STANDARD unweighted multiplicative Poisson EM (mirrors em_solver.solve_em),
on the kept-column subset, with the same densify / nz-filter / cov-scalar /
simplex handling as archive/sweep_shape_norm_h_only.py.

Outputs (under scratch/h_fixes/balanced/):
  balanced_diagnostic.npz   -- subset definition + per-founder / per-class balance
  balanced_<SIM>.npz        -- one per sim: h (231-vec), founders, cov, info
  balanced_scores.tsv       -- score table (g0 5 sims + SEEDMIX S1-8 + summary)
"""
from __future__ import annotations
import argparse, gc, json, os, sys, time
import numpy as np
from scipy.sparse import load_npz

# import the EXISTING standard unweighted EM (do not reimplement)
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from em_solver import solve_em  # noqa: E402


# --------------------------------------------------------------------------- #
# EM driver on a precomputed kept-column subset (mirrors sweep_shape_norm)     #
# --------------------------------------------------------------------------- #
def run_em_on_subset(cn_dense_sub: np.ndarray, counts_sub: np.ndarray,
                     em_max_iter: int = 300):
    """Standard unweighted EM on the (already column-subset) dense cn + counts.

    cn_dense_sub : F x Ksub dense float32   (kept-columns of filt2 cn)
    counts_sub   : Ksub vector              (kept-columns of filt2 counts)
    """
    F = cn_dense_sub.shape[0]
    # nz filter (mirror sweep_shape_norm_h_only): only columns with c_k > 0
    nz = counts_sub > 0
    cn_em = np.ascontiguousarray(cn_dense_sub[:, nz])
    counts_em = counts_sub[nz].astype(np.float32)
    # cov scalar (lambda) as in the reference driver; cancels in M-step
    cov = counts_sub.sum() * F / max(cn_dense_sub.sum(), 1)
    h, info = solve_em(counts_em, cn_em, cov, max_iter=em_max_iter, tol=1e-7)
    info["cov"] = float(cov)
    info["nz_cols"] = int(nz.sum())
    info["counts_sum"] = float(counts_em.sum())
    return h, info


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cn-prefix",
                    default="poolfreq/data/cn_full_231_v3qc_v3_filt2/cn_Chr1",
                    help="prefix; expects {prefix}.cn.npz and {prefix}.meta.npz")
    ap.add_argument("--g0-dir", default="scratch/g0_sweep_h_test")
    ap.add_argument("--seedmix-dir", default="scratch/seedmix_h_test")
    ap.add_argument("--out-dir", default="scratch/h_fixes/balanced")
    ap.add_argument("--split-json", default="data/founder_split_cactus_pg.json")
    ap.add_argument("--span", type=int, default=1,
                    help="keep k-mers whose bubble span (end-start) == this "
                         "(1 = single-base biallelic SNP). Use 0 to disable "
                         "and keep span <= --span-max instead.")
    ap.add_argument("--span-max", type=int, default=0,
                    help="if >0, keep k-mers with 1 <= span <= span-max "
                         "(relaxation; overrides --span exact match)")
    ap.add_argument("--em-max-iter", type=int, default=300)
    args = ap.parse_args()

    t0 = time.time()
    os.makedirs(args.out_dir, exist_ok=True)
    print(f"=== balanced sub-design EM  span={args.span} span_max={args.span_max} ===",
          flush=True)

    # ---- load cn + meta ----
    cn = load_npz(args.cn_prefix + ".cn.npz").tocsr()
    meta = np.load(args.cn_prefix + ".meta.npz", allow_pickle=True)
    founders = np.asarray(meta["founders"]).astype(str)
    bubble_id = np.asarray(meta["bubble_id"]).astype(np.int64)
    bstart = np.asarray(meta["bubble_start"]).astype(np.int64)
    bend = np.asarray(meta["bubble_end"]).astype(np.int64)
    F, K = cn.shape
    print(f"  cn: F={F} K={K:,} nnz={cn.nnz:,}  [{time.time()-t0:.0f}s]", flush=True)

    # ---- balanced subset mask (per-kmer via bubble_id) ----
    span_b = bend - bstart                       # per-bubble span (B-vector)
    span_k = span_b[bubble_id]                   # per-kmer span    (K-vector)
    if args.span_max > 0:
        keep = (span_k >= 1) & (span_k <= args.span_max)
        subset_def = f"1 <= bubble_span <= {args.span_max}"
    else:
        keep = (span_k == args.span)
        subset_def = f"bubble_span == {args.span}  (single-ref-base SNP locus)"
    n_keep = int(keep.sum())
    n_bub_keep = int(np.unique(bubble_id[keep]).size)
    print(f"  BALANCED SUBSET: {subset_def}", flush=True)
    print(f"    kept k-mers: {n_keep:,} / {K:,} ({100*n_keep/K:.2f}%)  "
          f"in {n_bub_keep:,} bubbles", flush=True)

    # ---- classes ----
    cls = json.load(open(args.split_json))
    cact_set = set(cls["cactus"]); pg_set = set(cls["PG"])
    is_cact = np.array([f in cact_set for f in founders])
    is_pg = np.array([f in pg_set for f in founders])
    print(f"    classes: cactus={is_cact.sum()} PG={is_pg.sum()}", flush=True)

    # ---- BALANCE DIAGNOSTIC: per-founder support on subset vs full ----
    # per-founder allele count = number of k-mers founder carries (cn row nnz),
    # computed on full design and on the kept subset.
    keep_idx = np.flatnonzero(keep)
    cn_keep = cn[:, keep_idx].tocsr()            # F x Ksub sparse (cheap)
    supp_full = np.asarray((cn > 0).sum(axis=1)).flatten().astype(np.float64)
    supp_sub = np.asarray((cn_keep > 0).sum(axis=1)).flatten().astype(np.float64)

    def cls_stats(v):
        c = v[is_cact]; p = v[is_pg]
        return (float(c.mean()), float(p.mean()),
                float(c.mean()/max(p.mean(), 1e-9)))
    fc, fp, fr = cls_stats(supp_full)
    sc, sp, sr = cls_stats(supp_sub)
    print(f"  --- per-founder k-mer support (carried-kmer count) ---", flush=True)
    print(f"    FULL design: cactus_mean={fc:,.0f}  PG_mean={fp:,.0f}  "
          f"cactus/PG ratio={fr:.3f}", flush=True)
    print(f"    SNP subset : cactus_mean={sc:,.0f}  PG_mean={sp:,.0f}  "
          f"cactus/PG ratio={sr:.3f}", flush=True)
    np.savez(os.path.join(args.out_dir, "balanced_diagnostic.npz"),
             founders=founders, subset_def=subset_def,
             span=args.span, span_max=args.span_max,
             n_keep=n_keep, n_bubbles_keep=n_bub_keep, K=K,
             keep_idx=keep_idx.astype(np.int64),
             is_cact=is_cact, is_pg=is_pg,
             supp_full=supp_full, supp_sub=supp_sub,
             full_cact_mean=fc, full_pg_mean=fp, full_ratio=fr,
             sub_cact_mean=sc, sub_pg_mean=sp, sub_ratio=sr)
    del cn_keep
    gc.collect()

    # ---- densify FULL cn once, then column-subset for EM ----
    t = time.time()
    cn_dense = np.asarray(cn.todense()).astype(np.float32)
    print(f"  densify FULL cn -> {cn_dense.nbytes/1e9:.1f} GB  [{time.time()-t:.0f}s]",
          flush=True)
    del cn
    gc.collect()
    cn_sub = np.ascontiguousarray(cn_dense[:, keep_idx])   # F x Ksub
    del cn_dense
    gc.collect()
    print(f"  cn_sub: {cn_sub.shape}  {cn_sub.nbytes/1e9:.1f} GB", flush=True)

    # ---- sims to run ----
    g0_sims = ["g0_n231_rep0_rand", "g0_n200_rep0_rand",
               "g0_n50_rep0_cact", "g0_n50_rep1_bal", "g0_n50_rep2_pg"]
    sm_sims = [f"S{i}" for i in range(1, 9)]

    jobs = []  # (name, counts_path)
    for s in g0_sims:
        jobs.append((s, os.path.join(args.g0_dir, f"filt2_{s}.counts.npy")))
    for s in sm_sims:
        jobs.append((s, os.path.join(args.seedmix_dir, f"filt2_{s}.counts.npy")))

    results = {}
    for name, cpath in jobs:
        if not os.path.exists(cpath):
            print(f"  !! MISSING counts: {cpath} -- skipping {name}", flush=True)
            continue
        counts_full = np.load(cpath)
        if counts_full.shape[0] != K:
            raise SystemExit(f"counts {name} len {counts_full.shape[0]} != K {K}")
        counts_sub = counts_full[keep_idx].astype(np.float32)
        t = time.time()
        h, info = run_em_on_subset(cn_sub, counts_sub, em_max_iter=args.em_max_iter)
        cm = float(h[is_cact].sum())
        results[name] = h
        np.savez(os.path.join(args.out_dir, f"balanced_{name}.npz"),
                 sample=name, founders=founders, h=h,
                 cov=info["cov"], iterations=info["iterations"],
                 converged=info["converged"],
                 nz_cols=info["nz_cols"], counts_sum=info["counts_sum"],
                 cactus_mass=cm, subset_def=subset_def, n_keep=n_keep)
        print(f"  [{name}] iter={info['iterations']} conv={info['converged']} "
              f"nz={info['nz_cols']:,} cactus_mass={cm:.4f} eff_n={1/(h**2).sum():.1f} "
              f"[{time.time()-t:.0f}s]", flush=True)

    np.savez(os.path.join(args.out_dir, "balanced_all_h.npz"),
             founders=founders,
             names=np.array(list(results.keys())),
             H=np.array([results[n] for n in results.keys()]))
    print(f"\nDONE balanced EM. {len(results)} sims. [TOTAL {time.time()-t0:.0f}s]",
          flush=True)


if __name__ == "__main__":
    main()
