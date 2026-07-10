#!/usr/bin/env python3
"""Window-mode counterpart to diag_convergence.py: does coverage matter per window?

GLOBAL mode pools k-mer evidence across the whole chromosome, so ĥ error was
coverage-INDEPENDENT (pure identifiability bias ~0.0485; see RESULTS_global_*).
WINDOW mode fits each ~10 kb window from ONLY its local k-mers (thousands, not
8.1M), so two things should change and are measured here:

  (A) IDENTIFIABILITY is worse locally — even noiseless (coverage=∞), per-window
      ĥ should be much further from h_true than the global 0.0485, because a
      handful of local k-mers cannot pin 80 founders.
  (B) COVERAGE now bites — at low λ many of a window's few k-mers get zero counts,
      so support collapses and per-window error grows as λ falls (unlike global).
      This is the mechanism behind low-info windows being NaN'd.

Clean generative test on p80, g0 pool (no recombination ⇒ local ancestry = global,
so the per-window truth is the SAME h_true at every window). For each window:
  noiseless:  c_k = μ_k(h_true)        over the window's k-mers
  coverage λ: c_k ~ Poisson(λ μ_k)     (R replicates, averaged)
  ĥ_w = per-window weighted-Poisson EM (ω=1/m_b); error = ||ĥ_w - h_true||.

Reports per-window errors + support, and the aggregate error-vs-coverage and
error-vs-support relationships. Run via sbatch.
"""
import argparse, sys, time
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import sparse

ROOT = Path("/global/scratch/users/tbellg/kmate")
sys.path.insert(0, str(ROOT / "src"))
from kmate.em_solver import solve_em                                   # noqa: E402
from kmate.block_em import define_windows, assign_kmers_to_blocks      # noqa: E402

PREFIX = ROOT / "benchmarks/p80/data/kmer_pa_p80_filt2"


def load_all():
    kp = sparse.load_npz(PREFIX / "kmer_pa_Chr1.kmer_pa.npz")
    kmer_pa = np.asarray(kp.todense(), dtype=np.float32) if sparse.issparse(kp) \
        else np.asarray(kp, dtype=np.float32)
    m = np.load(PREFIX / "kmer_pa_Chr1.meta.npz", allow_pickle=True)
    founders = np.asarray(m["founders"]).astype(str)
    bid = np.asarray(m["bubble_id"])
    m_b = np.bincount(bid)[bid].astype(np.float64)
    omega = (1.0 / m_b).astype(np.float32)
    meta = dict(bubble_id=bid, bubble_chrom=np.asarray(m["bubble_chrom"]).astype(str),
                bubble_start=np.asarray(m["bubble_start"]),
                bubble_end=np.asarray(m["bubble_end"]))
    return kmer_pa, founders, omega, meta


def load_truth_h(pool, founders):
    w = pd.read_csv(ROOT / f"benchmarks/p80/sims/{pool}/pool_weights.tsv", sep="\t")
    wmap = {str(f): float(x) for f, x in zip(w["founder"], w["weight"])}
    h = np.array([wmap.get(str(f), 0.0) for f in founders], dtype=np.float64)
    s = h.sum()
    return h / s if s > 0 else h


def win_em(c, K, omega, max_iter=500, tol=1e-9):
    """Per-window EM on c>0 k-mers only. Returns ĥ (F,) or None if no evidence."""
    nz = c > 0
    if nz.sum() < 1:
        return None
    om = None if omega is None else omega[nz]
    h, _ = solve_em(c[nz].astype(np.float32), K[:, nz], 1.0,
                    max_iter=max_iter, tol=tol, omega=om)
    return h


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pool", default="cov10_n50_g0_s42_hotspots_p80_chr1")
    ap.add_argument("--window-bp", type=int, default=10000)
    ap.add_argument("--covs", default="3,10,30")
    ap.add_argument("--reps", type=int, default=3)
    ap.add_argument("--min-kmers", type=int, default=200,
                    help="window support floor (windows below this would be NaN'd)")
    ap.add_argument("--out", default=str(ROOT / "benchmarks/h_uncertainty/results"))
    args = ap.parse_args()
    covs = [float(x) for x in args.covs.split(",")]
    outdir = Path(args.out); outdir.mkdir(parents=True, exist_ok=True)

    kmer_pa, founders, omega, meta = load_all()
    F, K = kmer_pa.shape
    h_true = load_truth_h(args.pool, founders)
    supp_true = h_true > 0
    mu_true = (h_true @ kmer_pa.astype(np.float64))     # per-k-mer expected rate
    print(f"p80 F={F} K={K:,}; pool={args.pool} true_support={int(supp_true.sum())}",
          flush=True)

    # windows + k-mer→window assignment
    blocks = define_windows(meta["bubble_chrom"], meta["bubble_start"],
                            meta["bubble_end"], window_bp=args.window_bp)
    kblock = assign_kmers_to_blocks(meta["bubble_id"], meta["bubble_chrom"],
                                    meta["bubble_start"], meta["bubble_end"], blocks)
    print(f"{len(blocks)} windows of {args.window_bp} bp; "
          f"{int((kblock>=0).sum()):,}/{K:,} k-mers assigned", flush=True)

    # group k-mers by window (sort once, take contiguous slices)
    order = np.argsort(kblock, kind="stable")
    kb_sorted = kblock[order]
    starts = np.searchsorted(kb_sorted, np.arange(len(blocks)), side="left")
    ends = np.searchsorted(kb_sorted, np.arange(len(blocks)), side="right")

    rng = np.random.default_rng(20260629)
    rows = []
    l2 = lambda h: float(np.linalg.norm(h - h_true))
    t0 = time.time()
    for wi in range(len(blocks)):
        idx = order[starts[wi]:ends[wi]]
        nk = idx.size
        if nk == 0:
            continue
        Kw = kmer_pa[:, idx]
        muw = mu_true[idx]
        omw = omega[idx]
        # how many TRUE founders even have a local k-mer here (local resolvability)
        n_loc_founders = int((Kw[supp_true].sum(axis=1) > 0).sum())

        rec = dict(window=wi, start=blocks[wi].start, n_kmers=nk,
                   n_local_truefounders=n_loc_founders,
                   below_floor=int(nk < args.min_kmers))
        # noiseless: identifiability floor for this window
        h_nl = win_em(muw.astype(np.float32), Kw, omw)
        rec["err_noiseless"] = l2(h_nl) if h_nl is not None else np.nan
        rec["nsupp_noiseless"] = int((h_nl > 1e-3).sum()) if h_nl is not None else 0
        # coverage sweep (R replicate Poisson draws, averaged)
        for lam in covs:
            errs, nzk = [], []
            for _ in range(args.reps):
                c = rng.poisson(lam * muw)
                nzk.append(int((c > 0).sum()))
                h = win_em(c.astype(np.float32), Kw, omw)
                errs.append(l2(h) if h is not None else np.nan)
            rec[f"err_cov{lam:g}"] = float(np.nanmean(errs))
            rec[f"nzk_cov{lam:g}"] = float(np.mean(nzk))
        rows.append(rec)
        if wi % 500 == 0:
            print(f"  window {wi}/{len(blocks)}  {time.time()-t0:.0f}s", flush=True)

    df = pd.DataFrame(rows)
    df.to_csv(outdir / "window_diag.tsv", sep="\t", index=False)

    # ---- aggregate summary ----
    print("\n=== WINDOW-MODE ĥ ERROR vs COVERAGE (per-window ||ĥ-h_true||) ===",
          flush=True)
    res = df[df.n_kmers >= args.min_kmers]      # resolvable windows
    def desc(col):
        x = df[col].dropna()
        return (f"median={x.median():.3f}  mean={x.mean():.3f}  "
                f"p90={x.quantile(.9):.3f}  frac>0.1={np.mean(x>0.1):.2f}")
    print(f"all {len(df)} windows ({(df.n_kmers<args.min_kmers).mean()*100:.0f}% "
          f"below {args.min_kmers}-kmer floor):", flush=True)
    print(f"  noiseless : {desc('err_noiseless')}", flush=True)
    for lam in covs:
        print(f"  cov {lam:g}    : {desc(f'err_cov{lam:g}')}", flush=True)
    print(f"\nGLOBAL-mode reference (whole chrom): ||ĥ-h_true||₂ ≈ 0.0485 (flat in cov)",
          flush=True)
    # coverage effect on resolvable windows: does low cov inflate error vs noiseless?
    print(f"\nOn the {len(res)} windows ≥{args.min_kmers} k-mers "
          f"(coverage-effect, identifiability subtracted):", flush=True)
    for lam in covs:
        d = (res[f"err_cov{lam:g}"] - res["err_noiseless"]).dropna()
        print(f"  cov {lam:g}: median(err_cov - err_noiseless) = {d.median():+.3f}",
              flush=True)
    # correlation error vs support
    sub = df.dropna(subset=["err_cov3"]) if "err_cov3" in df else df
    if len(sub) > 10:
        c = np.corrcoef(np.log10(sub.n_kmers + 1), sub["err_cov" + f"{covs[0]:g}"])[0, 1]
        print(f"\ncorr(log10 n_kmers, err_cov{covs[0]:g}) = {c:.3f} "
              f"(negative ⇒ more local k-mers → lower error)", flush=True)
    print(f"\nwrote {outdir}/window_diag.tsv", flush=True)


if __name__ == "__main__":
    main()
