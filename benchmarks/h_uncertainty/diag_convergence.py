#!/usr/bin/env python3
"""Diagnose the EM convergence / bias floor seen in the h-uncertainty benchmark.

The benchmark showed ||ĥ - h_true|| ≈ 0.048 at BOTH cov 3 and cov 10 (flat),
with the float32 EM running to max_iter (300) without tripping tol=1e-7. Three
candidate causes, separated here:

  (1) UNDER-ITERATION   — EM not run long enough. Test: many more iters, watch
      ||ĥ - h_true|| vs iteration. If it keeps dropping, it was under-iterated.
  (2) FLOAT32 PRECISION — tol=1e-7 on an 80-vector may sit below the float32
      round-off floor, so `delta < tol` never fires and the "non-convergence"
      is cosmetic. Test: same EM in float64; compare the ||Δh|| plateau and the
      final error.
  (3) IDENTIFIABILITY BIAS — the inverse problem is ill-posed (collinear
      founders); even with NOISELESS data (c_k = μ_k(h_true), no Poisson) and
      full convergence, ĥ ≠ h_true. This is the floor that no SE method around ĥ
      can ever cover. Test: the noiseless run. Its residual = pure bias.

Error decomposition per coverage:
  total error      = ||ĥ_cov - h_true||           (what the benchmark reported)
  identifiability  = ||ĥ_noiseless - h_true||      (coverage = ∞, no Poisson)
  coverage/variance≈ total - identifiability (in quadrature, roughly)

Heavy (p80: F=80, K=8.1M; many iters, float64 copy) → run via sbatch.
"""
import argparse, sys, time
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import sparse

ROOT = Path("/global/scratch/users/tbellg/kmate")
sys.path.insert(0, str(ROOT / "src"))
from kmate.em_solver import solve_em                                   # noqa: E402

PREFIX = ROOT / "benchmarks/p80/data/kmer_pa_p80_filt2"


def load_kmer_pa():
    kp = sparse.load_npz(PREFIX / "kmer_pa_Chr1.kmer_pa.npz")
    kmer_pa = np.asarray(kp.todense(), dtype=np.float32) if sparse.issparse(kp) \
        else np.asarray(kp, dtype=np.float32)
    meta = np.load(PREFIX / "kmer_pa_Chr1.meta.npz", allow_pickle=True)
    founders = np.asarray(meta["founders"]).astype(str)
    # GLOBAL mode: uniform kmer-weight (omega=None) per PIPELINE_STATE.md Sec.0
    # (2026-07-06) -- supersedes the old omega=1/m_b production weighting. Note
    # the "0.048 bias floor" this script investigates was later root-caused to
    # the Kf_w global-normalization bug (per_founder is now solve_em's default),
    # not pure identifiability -- re-run to see the floor mostly vanish.
    omega = None
    return kmer_pa, founders, omega


def load_truth_h(pool, founders):
    w = pd.read_csv(ROOT / f"benchmarks/p80/sims/{pool}/pool_weights.tsv", sep="\t")
    wmap = {str(f): float(x) for f, x in zip(w["founder"], w["weight"])}
    h = np.array([wmap.get(str(f), 0.0) for f in founders], dtype=np.float64)
    s = h.sum()
    return h / s if s > 0 else h


def em_f64(counts, kmer_pa64, omega64, max_iter, tol, h_true=None, trace_every=0):
    """Plain weighted-Poisson EM in float64 (same update as solve_em, no priors).
    Optionally records ||ĥ-h_true|| every `trace_every` iters."""
    F, K = kmer_pa64.shape
    h = np.full(F, 1.0 / F, dtype=np.float64)
    wc = counts.astype(np.float64) if omega64 is None else (omega64 * counts)
    total_c = wc.sum()
    dh_hist, err_hist = [], []
    it = 0
    for it in range(max_iter):
        denom = np.maximum(h @ kmer_pa64, 1e-12)
        cw = wc / denom
        em_term = h * (kmer_pa64 @ cw)
        h_new = em_term / max(total_c, 1e-300)
        h_new = h_new / h_new.sum()
        delta = float(np.linalg.norm(h_new - h))
        dh_hist.append(delta)
        if trace_every and (it % trace_every == 0) and h_true is not None:
            err_hist.append((it, float(np.linalg.norm(h_new - h_true))))
        h = h_new
        if delta < tol:
            break
    return h, {"iterations": it + 1, "delta_history": dh_hist,
               "err_history": err_hist, "converged": delta < tol}


def nz_filter(counts, kmer_pa, kmer_pa64, omega, omega64):
    """Restrict to k-mers with count>0 (EM M-step is exact on these; production
    does the same). Big speedup at low coverage; identical solution."""
    nz = counts > 0
    om = None if omega is None else omega[nz]
    om64 = None if omega64 is None else omega64[nz]
    return counts[nz], kmer_pa[:, nz], kmer_pa64[:, nz], om, om64, int(nz.sum())


def errs(h, h_true, supp_true):
    d = h - h_true
    return dict(
        l1=float(np.abs(d).sum()), l2=float(np.linalg.norm(d)),
        max_abs=float(np.abs(d).max()),
        l2_on_truesupp=float(np.linalg.norm(d[supp_true])),
        mass_off_truesupp=float(h[~supp_true].sum()),
        n_hat_supp=int((h > 1e-3).sum()),
    )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pool", default="cov10_n50_g0_s42_hotspots_p80_chr1")
    ap.add_argument("--covs", default="3,10,30,100")
    ap.add_argument("--max-iter", type=int, default=2000)
    ap.add_argument("--f64-covs", default="10",
                    help="coverages (besides noiseless) to also run in float64")
    ap.add_argument("--out", default=str(ROOT / "benchmarks/h_uncertainty/results"))
    args = ap.parse_args()
    covs = [float(x) for x in args.covs.split(",")]
    f64_covs = {float(x) for x in args.f64_covs.split(",") if x.strip()}
    outdir = Path(args.out); outdir.mkdir(parents=True, exist_ok=True)

    kmer_pa, founders, omega = load_kmer_pa()
    F, K = kmer_pa.shape
    h_true = load_truth_h(args.pool, founders)
    supp_true = h_true > 0
    mu_true = (h_true @ kmer_pa.astype(np.float64))   # expected per-kmer rate (λ=1)
    print(f"p80 F={F} K={K:,}; pool={args.pool} true_support={int(supp_true.sum())}",
          flush=True)

    kmer_pa64 = kmer_pa.astype(np.float64)            # for the f64 EM
    omega64 = None if omega is None else omega.astype(np.float64)
    rng = np.random.default_rng(20260629)
    rows, traces = [], {}

    # ---- (3) NOISELESS: pure identifiability floor (coverage = ∞) ----
    t0 = time.time()
    counts_nl = mu_true.astype(np.float32)            # c_k = μ_k(h_true), no Poisson
    c_n, kp_n, kp64_n, om_n, om64_n, nnz_n = nz_filter(
        counts_nl, kmer_pa, kmer_pa64, omega, omega64)
    h_nl_32, i32 = solve_em(c_n, kp_n, 1.0, max_iter=args.max_iter,
                            tol=1e-12, omega=om_n)
    h_nl_64, i64 = em_f64(c_n, kp64_n, om64_n, args.max_iter, 1e-14,
                          h_true=h_true, trace_every=25)
    e32, e64 = errs(h_nl_32, h_true, supp_true), errs(h_nl_64, h_true, supp_true)
    rows.append(dict(cov="noiseless", prec="f32", iters=i32["iterations"],
                     conv=i32["converged"], final_dh=i32["delta_history"][-1], **e32))
    rows.append(dict(cov="noiseless", prec="f64", iters=i64["iterations"],
                     conv=i64["converged"], final_dh=i64["delta_history"][-1], **e64))
    traces["noiseless_f64_dh"] = np.array(i64["delta_history"])
    traces["noiseless_f64_err"] = np.array(i64["err_history"])
    print(f"[noiseless] f32 err_l2={e32['l2']:.4f} (it {i32['iterations']}, "
          f"conv={i32['converged']}, dh={i32['delta_history'][-1]:.1e}) | "
          f"f64 err_l2={e64['l2']:.4f} (it {i64['iterations']}, "
          f"conv={i64['converged']}, dh={i64['delta_history'][-1]:.1e}) "
          f"mass_off_truesupp f64={e64['mass_off_truesupp']:.4f}  {time.time()-t0:.0f}s",
          flush=True)

    # ---- (1)+(2) COVERAGE SWEEP with Poisson noise (f32 always; f64 on request) ----
    for lam in covs:
        t0 = time.time()
        counts = rng.poisson(lam * mu_true).astype(np.float32)
        c_f, kp_f, kp64_f, om_f, om64_f, nnz = nz_filter(
            counts, kmer_pa, kmer_pa64, omega, omega64)
        h32, j32 = solve_em(c_f, kp_f, lam, max_iter=args.max_iter,
                            tol=1e-12, omega=om_f)
        e32 = errs(h32, h_true, supp_true)
        rows.append(dict(cov=lam, prec="f32", iters=j32["iterations"],
                         conv=j32["converged"], final_dh=j32["delta_history"][-1],
                         nnz=nnz, **e32))
        msg = (f"[cov {lam:g}] nz={nnz:,} | f32 err_l2={e32['l2']:.4f} "
               f"(it {j32['iterations']}, conv={j32['converged']}, "
               f"dh={j32['delta_history'][-1]:.1e})")
        if lam in f64_covs:
            h64, j64 = em_f64(c_f, kp64_f, om64_f, args.max_iter, 1e-14,
                              h_true=h_true, trace_every=25)
            e64 = errs(h64, h_true, supp_true)
            rows.append(dict(cov=lam, prec="f64", iters=j64["iterations"],
                             conv=j64["converged"], final_dh=j64["delta_history"][-1],
                             nnz=nnz, **e64))
            traces[f"cov{lam:g}_f64_dh"] = np.array(j64["delta_history"])
            traces[f"cov{lam:g}_f64_err"] = np.array(j64["err_history"])
            msg += (f" | f64 err_l2={e64['l2']:.4f} (it {j64['iterations']}, "
                    f"conv={j64['converged']}, dh={j64['delta_history'][-1]:.1e})")
        print(msg + f"  {time.time()-t0:.0f}s", flush=True)

    df = pd.DataFrame(rows)
    df.to_csv(outdir / "convergence_diag.tsv", sep="\t", index=False)
    np.savez(outdir / "convergence_traces.npz", **traces)
    print("\n=== CONVERGENCE / BIAS DIAGNOSTIC ===", flush=True)
    print(df.to_string(index=False), flush=True)
    nl_rows = df[(df["cov"].astype(str) == "noiseless") & (df.prec == "f64")]["l2"]
    if len(nl_rows):
        print(f"\nIdentifiability floor (noiseless f64 ||ĥ-h_true||_2) = "
              f"{nl_rows.iloc[0]:.4f}", flush=True)
    print("If coverage-sweep f64 errors ≈ this floor → bias-dominated, not "
          "coverage. If they fall toward 0 as cov↑ → coverage-dominated.", flush=True)
    print(f"\nwrote {outdir}/convergence_diag.tsv + convergence_traces.npz", flush=True)


if __name__ == "__main__":
    main()
