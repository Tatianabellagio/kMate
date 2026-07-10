#!/usr/bin/env python3
"""AF-functional certainty: does the EM bias that kills h-certainty CANCEL in the
identifiable AF functional, and is the Fisher delta-method AF SE calibrated?

Premise (from the metric/CRLB result): per-FOUNDER ĥ error is dominated by a
coverage-independent EM bias on the flat likelihood manifold (global 0.0485), which
NO variance method (Fisher/bootstrap/CRLB) captures (h bootstrap 95% covers truth 6%).
But AF_r = hᵀv_r / hᵀu_r sums over carriers, so swaps within the flat manifold cancel
=> AF should be (a) almost UNbiased at infinite coverage and (b) variance-dominated,
so the Fisher delta-method SE (af_se_from_cov) should be CALIBRATED.

This validates that claim, decomposing per-record AF error into
  bias/identifiability = AF(ĥ_noiseless) - AF(h_true)        target ~0 (vs h's 0.0485)
  coverage/variance    = AF(ĥ_cov λ)     - AF(ĥ_noiseless)
and measuring the Fisher AF-SE 95% coverage of the TRUE AF (vs h's 6%), overall and
for well-called records. v=var_pa, u=var_called. p80 Chr1 global. Run via sbatch.
"""
import sys, time
from pathlib import Path
import numpy as np, pandas as pd
from scipy import sparse

ROOT = Path("/global/scratch/users/tbellg/kmate")
sys.path.insert(0, str(ROOT / "src"))
from kmate.em_solver import solve_em                                    # noqa
from kmate.h_uncertainty import fisher_cov_h, af_se_from_cov           # noqa
KPRE = ROOT / "benchmarks/p80/data/kmer_pa_p80_filt2"
VPRE = ROOT / "benchmarks/p80/data"
COVS = [3.0, 10.0, 30.0]
REPS = 3
SEED = 20260630


def load_kmers():
    kp = sparse.load_npz(KPRE / "kmer_pa_Chr1.kmer_pa.npz")
    K = np.asarray(kp.todense(), np.float32) if sparse.issparse(kp) else np.asarray(kp, np.float32)
    m = np.load(KPRE / "kmer_pa_Chr1.meta.npz", allow_pickle=True)
    fnd = np.asarray(m["founders"]).astype(str)
    # GLOBAL mode: uniform kmer-weight (omega=None) per PIPELINE_STATE.md Sec.0
    # (2026-07-06) -- supersedes the old omega=1/m_b production weighting.
    omega = None
    return K, fnd, omega


def truth_h(fnd, pool="cov10_n50_g0_s42_hotspots_p80_chr1"):
    w = pd.read_csv(ROOT / f"benchmarks/p80/sims/{pool}/pool_weights.tsv", sep="\t")
    wm = {str(f): float(x) for f, x in zip(w["founder"], w["weight"])}
    h = np.array([wm.get(str(f), 0.0) for f in fnd], np.float64); return h / h.sum()


def load_var(fnd_k):
    V0 = np.asarray(sparse.load_npz(VPRE / "var_pa_p80.var_pa.npz").todense(), np.float64)
    U0 = np.asarray(sparse.load_npz(VPRE / "var_pa_p80.var_called.npz").todense(), np.float64)
    fnd_v = np.asarray(np.load(VPRE / "var_pa_p80.meta.npz", allow_pickle=True)["founders"]).astype(str)
    pos = {f: i for i, f in enumerate(fnd_v)}
    ridx = np.array([pos[f] for f in fnd_k])          # reorder var rows -> kmer founder order
    return V0[ridx], U0[ridx]


def project(h, V, U, eps=1e-12):
    den = np.maximum(h @ U, eps); return (h @ V) / den, (h @ U)


def fit(c, K, omega, cov):
    nz = c > 0
    om_nz = None if omega is None else omega[nz]
    h, _ = solve_em(c[nz].astype(np.float32), K[:, nz], cov,
                    max_iter=1500, tol=1e-9, omega=om_nz)
    return h


def cover_stats(af_hat, af_true, se, mask):
    ok = np.abs(af_hat - af_true) <= 1.96 * np.maximum(se, 1e-12)
    d = np.abs(af_hat - af_true); m = mask & np.isfinite(se) & np.isfinite(d)
    disc = np.corrcoef(se[m], d[m])[0, 1] if m.sum() > 5 else np.nan
    return float(ok[mask].mean()), float(disc)


def main():
    t0 = time.time()
    K, fnd, omega = load_kmers()
    h_true = truth_h(fnd); supp = np.flatnonzero(h_true > 0)
    mu = (h_true @ K).astype(np.float64)
    V, U = load_var(fnd)
    print(f"p80 F={K.shape[0]} K={K.shape[1]:,} records={V.shape[1]:,} "
          f"support={supp.size}  load {time.time()-t0:.0f}s", flush=True)

    # truth AF + record masks
    af_true, den_true = project(h_true, V, U)
    defined = den_true > 0                              # called by >=1 present founder
    callrate = (U[supp] > 0).sum(0) / supp.size         # fraction of present founders calling r
    wellcalled = defined & (callrate >= 0.99)
    print(f"records defined={defined.sum():,}  well-called(>=0.99)={wellcalled.sum():,}", flush=True)

    rng = np.random.default_rng(SEED)
    rows = []

    # (1) noiseless: the IDENTIFIABILITY/BIAS floor of the AF functional
    h_nl = fit(mu, K, omega, 1.0)
    af_nl, _ = project(h_nl, V, U)
    h_bias = float(np.linalg.norm(h_nl - h_true))
    af_bias = np.abs(af_nl - af_true)
    print(f"\n=== (1) IDENTIFIABILITY (noiseless, cov=inf) ===", flush=True)
    print(f"  h  bias  ||ĥ-h_true|| = {h_bias:.4f}   (the un-fixable per-founder floor)", flush=True)
    print(f"  AF bias  mean|ΔAF| defined   = {af_bias[defined].mean():.5f}   "
          f"median = {np.median(af_bias[defined]):.5f}   max = {af_bias[defined].max():.4f}", flush=True)
    print(f"  AF bias  mean|ΔAF| wellcalled= {af_bias[wellcalled].mean():.5f}", flush=True)
    print(f"  -> does the h bias CANCEL in AF?  ratio AFmean/hbias = {af_bias[defined].mean()/h_bias:.3f}", flush=True)

    # (2) coverage sweep: variance + Fisher AF-SE calibration
    print(f"\n=== (2) COVERAGE + Fisher AF-SE calibration ===", flush=True)
    for lam in COVS:
        cov_all, well_all, disc_all, mae_all, exc_all = [], [], [], [], []
        for r in range(REPS):
            c = rng.poisson(lam * mu)
            h = fit(c.astype(np.float32), K, omega, lam)
            af_hat, _ = project(h, V, U)
            Sig, sp = fisher_cov_h(h, K, c.astype(np.float64), omega=omega)
            _, se = af_se_from_cov(h, Sig, V, U, support=sp)
            c95, disc = cover_stats(af_hat, af_true, se, defined)
            w95, _ = cover_stats(af_hat, af_true, se, wellcalled)
            cov_all.append(c95); well_all.append(w95); disc_all.append(disc)
            mae_all.append(float(np.abs(af_hat - af_true)[defined].mean()))
            exc_all.append(float((np.abs(af_hat-af_true)-af_bias)[defined].mean()))
        rows.append(dict(cov=lam, af_mae=np.mean(mae_all), af_cov_excess=np.mean(exc_all),
                         fisher95_defined=np.mean(cov_all), fisher95_wellcalled=np.mean(well_all),
                         se_discrim_corr=np.nanmean(disc_all)))
        print(f"  cov {lam:>4g}:  AF MAE={np.mean(mae_all):.5f}  cov-excess={np.mean(exc_all):+.5f}  "
              f"Fisher95(defined)={np.mean(cov_all):.3f}  Fisher95(wellcalled)={np.mean(well_all):.3f}  "
              f"corr(SE,|err|)={np.nanmean(disc_all):+.3f}", flush=True)

    df = pd.DataFrame(rows)
    df.to_csv(ROOT / "benchmarks/h_uncertainty/results/af_certainty.tsv", sep="\t", index=False)
    print(f"\n=== VERDICT ===", flush=True)
    print(f"  h: bias {h_bias:.3f}, bootstrap 95% covers truth ~6%  (variance != error)", flush=True)
    print(f"  AF: noiseless bias {af_bias[defined].mean():.5f} (well-called {af_bias[wellcalled].mean():.5f}); "
          f"Fisher 95% covers ~{df.fisher95_wellcalled.mean():.2f}", flush=True)
    print(f"  => if AF bias<<h bias AND Fisher95~0.95: certify AF with delta-method SE. ({time.time()-t0:.0f}s)", flush=True)


if __name__ == "__main__":
    main()
