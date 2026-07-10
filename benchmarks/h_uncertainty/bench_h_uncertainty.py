#!/usr/bin/env python3
"""Benchmark: analytic Fisher SE on ĥ vs the parametric bootstrap (gold).

Clean generative test on the p80 panel (F=80, ~8.1M k-mers). For a known true
mixture h_true and a chosen coverage λ:

  counts  c_k ~ Poisson(λ · μ_k(h_true))          (the data)
  ĥ      = per_founder-normalized Poisson EM(counts, ω=uniform)  (production
           GLOBAL-mode estimator, 2026-07-06: Kf_w fix, drops ω=1/m_b)
  Fisher  SE = sqrt diag of fisher_cov_h(ĥ, ...)  (analytic, cheap)
  Boot    SE = sqrt diag cov of B EM-refits on c* ~ Poisson(λ̂ μ_k(ĥ))   (gold)
  TruthMC SE = same but resampling around h_true  (the true sampling dist)

Bootstrap and TruthMC both assume the same Poisson model, so they bound what
ANY honest SE should report; Fisher is the analytic approximation we want to
validate (and find the rcond that makes it track). Sweeps coverage to confirm
low-cov → wider SE and that Fisher follows. Writes per-founder SEs + a summary.

Heavy (B EM refits on 8.1M k-mers per coverage) → run via sbatch.
"""
import argparse, sys, time
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import sparse

ROOT = Path("/global/scratch/users/tbellg/kmate")
sys.path.insert(0, str(ROOT / "src"))
from kmate.em_solver import solve_em                                   # noqa: E402
from kmate.h_uncertainty import (fisher_information_h,                 # noqa: E402
                                 _tangent_pinv_on_support,
                                 bootstrap_cov_h, identifiability)

PREFIX = ROOT / "benchmarks/p80/data/kmer_pa_p80_filt2"


def load_kmer_pa():
    kp = sparse.load_npz(PREFIX / "kmer_pa_Chr1.kmer_pa.npz")
    kmer_pa = np.asarray(kp.todense(), dtype=np.float32) if sparse.issparse(kp) \
        else np.asarray(kp, dtype=np.float32)
    meta = np.load(PREFIX / "kmer_pa_Chr1.meta.npz", allow_pickle=True)
    founders = np.asarray(meta["founders"]).astype(str)
    # GLOBAL mode: uniform kmer-weight (omega=None) per PIPELINE_STATE.md Sec.0
    # (2026-07-06) -- supersedes the old omega=1/m_b production weighting.
    omega = None
    return kmer_pa, founders, omega


def load_truth_h(pool, founders):
    w = pd.read_csv(ROOT / f"benchmarks/p80/sims/{pool}/pool_weights.tsv", sep="\t")
    wmap = {str(f): float(x) for f, x in zip(w["founder"], w["weight"])}
    h = np.array([wmap.get(str(f), 0.0) for f in founders], dtype=np.float64)
    s = h.sum()
    return h / s if s > 0 else h


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pool", default="cov10_n50_g0_s42_hotspots_p80_chr1")
    ap.add_argument("--covs", default="3,10,30")
    ap.add_argument("--B", type=int, default=120)
    ap.add_argument("--rconds", default="1e-10,1e-4,1e-2,5e-2")
    ap.add_argument("--support-eps", type=float, default=1e-3)
    ap.add_argument("--out", default=str(ROOT / "benchmarks/h_uncertainty/results"))
    args = ap.parse_args()
    covs = [float(x) for x in args.covs.split(",")]
    rconds = [float(x) for x in args.rconds.split(",")]
    outdir = Path(args.out); outdir.mkdir(parents=True, exist_ok=True)

    kmer_pa, founders, omega = load_kmer_pa()
    F, K = kmer_pa.shape
    h_true = load_truth_h(args.pool, founders)
    print(f"p80 kmer_pa F={F} K={K:,}; pool={args.pool} "
          f"true support={int((h_true>0).sum())}", flush=True)

    rng = np.random.default_rng(20260629)
    rows, summ = [], []
    for lam in covs:
        t0 = time.time()
        mu_t = h_true @ kmer_pa
        counts = rng.poisson(lam * mu_t).astype(np.float32)
        nnz = int((counts > 0).sum())
        h_hat, info = solve_em(counts, kmer_pa, lam, max_iter=300, tol=1e-7,
                               omega=omega)
        print(f"[cov {lam:g}] nz={nnz:,} EM {info['iterations']} it; "
              f"||ĥ-h_true||={np.linalg.norm(h_hat-h_true):.4f}", flush=True)

        # gold: bootstrap (around ĥ) and truth-MC (around h_true)
        Sig_b, _ = bootstrap_cov_h(h_hat, kmer_pa, counts, omega=omega,
                                   B=args.B, coverage=lam, seed=101,
                                   max_iter=120, tol=1e-5)
        se_b = np.sqrt(np.clip(np.diag(Sig_b), 0, None))
        Sig_t, _ = bootstrap_cov_h(h_hat, kmer_pa, counts, omega=omega,
                                   B=args.B, coverage=lam, seed=202, rate_h=h_true,
                                   max_iter=120, tol=1e-5)
        se_t = np.sqrt(np.clip(np.diag(Sig_t), 0, None))

        # analytic Fisher at a grid of rcond
        J = fisher_information_h(h_hat, kmer_pa, counts, omega=omega)
        supp = np.flatnonzero(h_hat > args.support_eps)
        se_f = {}
        for rc in rconds:
            Sig_f = _tangent_pinv_on_support(J, supp, rcond=rc)
            se_f[rc] = np.sqrt(np.clip(np.diag(Sig_f), 0, None))

        ident = identifiability(h_hat, kmer_pa, counts, omega=omega,
                                support_eps=args.support_eps)
        for f in range(F):
            row = dict(cov=lam, founder=founders[f], h_true=h_true[f],
                       h_hat=h_hat[f], se_boot=se_b[f], se_truthmc=se_t[f],
                       in_support=int(f in supp))
            for rc in rconds:
                row[f"se_fisher_rc{rc:g}"] = se_f[rc][f]
            rows.append(row)

        # summary: corr & ratio Fisher-vs-bootstrap on support; calibration
        m = (se_b > 1e-7) & np.isin(np.arange(F), supp)
        for rc in rconds:
            corr = np.corrcoef(se_f[rc][m], se_b[m])[0, 1] if m.sum() > 2 else np.nan
            ratio = float(np.median(se_f[rc][m] / np.maximum(se_b[m], 1e-12)))
            summ.append(dict(cov=lam, rcond=rc, n_support=int(m.sum()),
                             corr_fisher_boot=corr, med_ratio_fisher_boot=ratio,
                             eff_rank=ident["eff_rank"], cond=ident["cond"]))
        # bootstrap calibration: does ĥ±1.96·se_boot cover h_true?
        cov_in = np.abs(h_hat - h_true) <= 1.96 * np.maximum(se_b, 1e-12)
        boot_vs_truth = np.corrcoef(se_b[m], se_t[m])[0, 1] if m.sum() > 2 else np.nan
        print(f"[cov {lam:g}] boot~truthMC corr={boot_vs_truth:.3f}; "
              f"95% cover(h_true) on support={cov_in[supp].mean():.2f}; "
              f"{time.time()-t0:.0f}s", flush=True)

    pd.DataFrame(rows).to_csv(outdir / "per_founder_se.tsv", sep="\t", index=False)
    sdf = pd.DataFrame(summ)
    sdf.to_csv(outdir / "summary.tsv", sep="\t", index=False)
    print("\n=== SUMMARY (Fisher vs bootstrap, on support) ===", flush=True)
    print(sdf.to_string(index=False), flush=True)
    print(f"\nwrote {outdir}/per_founder_se.tsv + summary.tsv", flush=True)


if __name__ == "__main__":
    main()
