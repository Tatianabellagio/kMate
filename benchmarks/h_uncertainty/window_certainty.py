#!/usr/bin/env python3
"""WINDOW-mode certainty: decompose per-window ĥ error and test what predicts it.

Global mode is identifiability-bias dominated and coverage-independent, so neither
Fisher nor bootstrap SE is a usable certainty there (both blind to the bias; see
RESULTS_global_convergence + Fisher-vs-bootstrap summary). WINDOW mode is the
opposite regime: each ~10 kb window is fit from only its local k-mers, so BOTH
(a) local identifiability (how many founders the window's k-mers can separate) and
(b) coverage/support (how many of those k-mers actually got counts) drive the error.

This script measures, per window, whether the two error components are PREDICTABLE
from local quantities — i.e. whether we can REPORT a certainty:

  err_noiseless(w)              = ||ĥ_w(c=μ) - h_true||      identifiability floor
  err_cov(w,λ)                 = ||ĥ_w(Poisson λ) - h_true|| total error at cov λ
  coverage_excess(w,λ)         = err_cov - err_noiseless     the variance part

  identifiability predictors:  n_kmers, n_local_truefounders, Fisher eff_rank, cond(J_w)
  variance predictor:          Fisher total SE  sqrt(Σ diag Σ_w)   [+ window bootstrap on a sample]

Aggregate question 1 (identifiability): does eff_rank / cond(J_w) predict err_noiseless?
Aggregate question 2 (variance):       does Fisher/boot SE predict coverage_excess, and
                                       does a per-window 95% interval COVER h_true better
                                       than global's 6%?

g0 pool (no recomb ⇒ per-window truth = global h_true). p80, Chr1. Run via sbatch.
"""
import argparse, sys, time
from pathlib import Path
import numpy as np, pandas as pd
from scipy import sparse

ROOT = Path("/global/scratch/users/tbellg/kmate")
sys.path.insert(0, str(ROOT / "src"))
from kmate.em_solver import solve_em                                    # noqa
from kmate.block_em import define_windows, assign_kmers_to_blocks       # noqa
from kmate.h_uncertainty import (fisher_information_h,                  # noqa
                                 _tangent_pinv_on_support, bootstrap_cov_h)
PRE = ROOT / "benchmarks/p80/data/kmer_pa_p80_filt2"


def load_all():
    kp = sparse.load_npz(PRE / "kmer_pa_Chr1.kmer_pa.npz")
    K = np.asarray(kp.todense(), dtype=np.float32) if sparse.issparse(kp) else np.asarray(kp, np.float32)
    m = np.load(PRE / "kmer_pa_Chr1.meta.npz", allow_pickle=True)
    fnd = np.asarray(m["founders"]).astype(str); bid = np.asarray(m["bubble_id"])
    omega = (1.0 / np.bincount(bid)[bid]).astype(np.float32)
    meta = dict(bid=bid, chrom=np.asarray(m["bubble_chrom"]).astype(str),
                start=np.asarray(m["bubble_start"]), end=np.asarray(m["bubble_end"]))
    return K, fnd, omega, meta


def truth(pool, fnd):
    w = pd.read_csv(ROOT / f"benchmarks/p80/sims/{pool}/pool_weights.tsv", sep="\t")
    wm = {str(f): float(x) for f, x in zip(w["founder"], w["weight"])}
    h = np.array([wm.get(str(f), 0.0) for f in fnd], dtype=np.float64); return h / h.sum()


def eff_rank_cond(J, supp):
    s = np.asarray(supp)
    if s.size <= 1:
        return float(s.size), np.inf
    ev = np.clip(np.linalg.eigvalsh(J[np.ix_(s, s)]), 0, None)
    pos = ev[ev > ev.max() * 1e-12] if ev.max() > 0 else ev[:0]
    if pos.size == 0:
        return 0.0, np.inf
    p = pos / pos.sum()
    return float(np.exp(-(p * np.log(p)).sum())), float(pos.max() / pos.min())


def win_em(c, K, om, max_iter=400, tol=1e-9):
    nz = c > 0
    if nz.sum() < 1:
        return None
    h, _ = solve_em(c[nz].astype(np.float32), K[:, nz], 1.0, max_iter=max_iter, tol=tol,
                    omega=None if om is None else om[nz])
    return h


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pool", default="cov10_n50_g0_s42_hotspots_p80_chr1")
    ap.add_argument("--window-bp", type=int, default=10000)
    ap.add_argument("--covs", default="3,10,30")
    ap.add_argument("--min-kmers", type=int, default=200)
    ap.add_argument("--boot-sample", type=int, default=150,
                    help="# windows (>=min-kmers) to also get a window-bootstrap SE on")
    ap.add_argument("--boot-B", type=int, default=40)
    ap.add_argument("--out", default=str(ROOT / "benchmarks/h_uncertainty/results"))
    args = ap.parse_args()
    covs = [float(x) for x in args.covs.split(",")]
    outdir = Path(args.out); outdir.mkdir(parents=True, exist_ok=True)

    K, fnd, omega, meta = load_all()
    F, Kn = K.shape
    h_true = truth(args.pool, fnd); supp_true = np.flatnonzero(h_true > 0)
    mu = (h_true @ K).astype(np.float32)
    print(f"p80 F={F} K={Kn:,} support={supp_true.size}", flush=True)

    blocks = define_windows(meta["chrom"], meta["start"], meta["end"], window_bp=args.window_bp)
    kb = assign_kmers_to_blocks(meta["bid"], meta["chrom"], meta["start"], meta["end"], blocks)
    order = np.argsort(kb, kind="stable"); kbs = kb[order]
    st = np.searchsorted(kbs, np.arange(len(blocks)), "left")
    en = np.searchsorted(kbs, np.arange(len(blocks)), "right")
    print(f"{len(blocks)} windows of {args.window_bp}bp", flush=True)

    rng = np.random.default_rng(20260629)
    boot_idx = set()  # windows to also bootstrap
    l2 = lambda h: float(np.linalg.norm(h - h_true))
    rows = []; t0 = time.time()
    # choose bootstrap sample windows up front (evenly spaced, >=min-kmers)
    elig = [w for w in range(len(blocks)) if (en[w]-st[w]) >= args.min_kmers]
    if elig:
        step = max(1, len(elig)//max(1, args.boot_sample))
        boot_idx = set(elig[::step][:args.boot_sample])

    for wi in range(len(blocks)):
        idx = order[st[wi]:en[wi]]; nk = idx.size
        if nk == 0:
            continue
        Kw = K[:, idx]; muw = mu[idx]; omw = omega[idx]
        nloc = int((Kw[supp_true].sum(axis=1) > 0).sum())
        rec = dict(window=wi, n_kmers=nk, n_local_truefounders=nloc,
                   below_floor=int(nk < args.min_kmers))
        # identifiability: noiseless fit + Fisher eigenstructure (at noiseless counts)
        h_nl = win_em(muw.astype(np.float32), Kw, omw)
        rec["err_noiseless"] = l2(h_nl) if h_nl is not None else np.nan
        if h_nl is not None:
            supp_w = np.flatnonzero(h_nl > 1e-3)
            J = fisher_information_h(h_nl, Kw, muw, omega=omw)
            er, cd = eff_rank_cond(J, supp_w)
            rec["eff_rank"] = er; rec["cond"] = cd; rec["nsupp"] = int(supp_w.size)
        # coverage sweep: error + Fisher SE prediction
        for lam in covs:
            c = rng.poisson(lam * muw)
            h = win_em(c.astype(np.float32), Kw, omw)
            rec[f"err_cov{lam:g}"] = l2(h) if h is not None else np.nan
            rec[f"nzk_cov{lam:g}"] = int((c > 0).sum())
            if h is not None:
                sw = np.flatnonzero(h > 1e-3)
                J = fisher_information_h(h, Kw, c.astype(np.float64), omega=omw)
                Sig = _tangent_pinv_on_support(J, sw, rcond=1e-2)
                se = np.sqrt(np.clip(np.diag(Sig), 0, None))
                rec[f"fisherSEtot_cov{lam:g}"] = float(np.sqrt((se**2).sum()))
                # coverage of truth by per-founder 95% interval (on support)
                if sw.size:
                    cov_ok = np.abs(h[sw] - h_true[sw]) <= 1.96*np.maximum(se[sw], 1e-12)
                    rec[f"cover95_cov{lam:g}"] = float(cov_ok.mean())
        # window bootstrap SE (gold) on the sampled windows, at the middle coverage
        if wi in boot_idx:
            lam = covs[len(covs)//2]
            c = rng.poisson(lam * muw)
            h = win_em(c.astype(np.float32), Kw, omw)
            if h is not None and (c > 0).sum() >= 2:
                Sb, _ = bootstrap_cov_h(h, Kw, c.astype(np.float64), omega=omw,
                                        B=args.boot_B, coverage=lam, seed=7,
                                        max_iter=200, tol=1e-5)
                seb = np.sqrt(np.clip(np.diag(Sb), 0, None))
                rec[f"bootSEtot_cov{lam:g}"] = float(np.sqrt((seb**2).sum()))
        rows.append(rec)
        if wi % 500 == 0:
            print(f"  window {wi}/{len(blocks)}  {time.time()-t0:.0f}s", flush=True)

    df = pd.DataFrame(rows); df.to_csv(outdir / "window_certainty.tsv", sep="\t", index=False)
    res = df[df.n_kmers >= args.min_kmers].copy()
    def corr(a, b):
        d = res[[a, b]].replace([np.inf, -np.inf], np.nan).dropna()
        return np.corrcoef(d[a], d[b])[0, 1] if len(d) > 5 else np.nan
    midlam = covs[len(covs)//2]
    print("\n=== WINDOW CERTAINTY (windows ≥%d k-mers, n=%d) ===" % (args.min_kmers, len(res)), flush=True)
    print("[identifiability] does local structure predict the noiseless floor err_noiseless?", flush=True)
    print(f"  corr(log10 n_kmers,        err_noiseless) = {corr_log('n_kmers','err_noiseless',res):.3f}", flush=True)
    print(f"  corr(eff_rank,             err_noiseless) = {corr('eff_rank','err_noiseless'):.3f}", flush=True)
    print(f"  corr(log10 cond,           err_noiseless) = {corr_log('cond','err_noiseless',res):.3f}", flush=True)
    print(f"  corr(n_local_truefounders, err_noiseless) = {corr('n_local_truefounders','err_noiseless'):.3f}", flush=True)
    print("\n[variance] does Fisher/boot SE predict the coverage excess (err_cov - err_noiseless)?", flush=True)
    for lam in covs:
        res[f"excess_cov{lam:g}"] = res[f"err_cov{lam:g}"] - res["err_noiseless"]
        print(f"  cov {lam:g}: corr(fisherSEtot, excess) = {corr(f'fisherSEtot_cov{lam:g}', f'excess_cov{lam:g}'):.3f}"
              f" ; mean cover95 = {res[f'cover95_cov{lam:g}'].mean():.2f}", flush=True)
    bc = f"bootSEtot_cov{midlam:g}"
    if bc in res:
        d = res[[f"fisherSEtot_cov{midlam:g}", bc]].dropna()
        if len(d) > 5:
            print(f"\n  Fisher vs window-bootstrap SE (cov {midlam:g}, n={len(d)}): "
                  f"corr={np.corrcoef(d.iloc[:,0], d.iloc[:,1])[0,1]:.3f} "
                  f"med ratio={np.median(d.iloc[:,0]/np.maximum(d.iloc[:,1],1e-9)):.2f}", flush=True)
    print(f"\nGLOBAL ref: err 0.0485 (coverage-flat); global bootstrap 95% cover ≈ 0.06.", flush=True)
    print(f"wrote {outdir}/window_certainty.tsv", flush=True)


def corr_log(a, b, res):
    import numpy as np
    d = res[[a, b]].replace([np.inf, -np.inf], np.nan).dropna()
    return np.corrcoef(np.log10(d[a] + 1), d[b])[0, 1] if len(d) > 5 else np.nan


if __name__ == "__main__":
    main()
