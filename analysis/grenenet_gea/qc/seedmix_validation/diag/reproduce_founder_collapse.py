#!/usr/bin/env python3
"""Minimal, read-free reproducer for kMate global-mode founder collapse.

Ideal noiseless experiment on Chr1: with h_true = uniform (1/231), the ideal
k-mer counts are c_k = mu_k = sum_f h_f K[f,k] = ac_k / 231.  Feeding these
EXACT expected counts back into the same EM the production runner uses isolates
the collapse to its cause: if uniform is not recovered from its own noiseless
expectation, the collapse is *pure identifiability* of the panel geometry, not
read noise or coverage.

Ladder:
  1. filt2inv panel (2<=ac<=230), omega=None       -> production panel, MLE
  2. filt2inv panel, omega=1/m_b                    -> production config (bubble wt)
  3. raw p231 panel (ac>=1, keeps private), None    -> does keeping ac==1 rescue?
  4. raw p231 panel, omega=1/m_b
  (+) filt2inv + Poisson read noise at realistic coverage -> does noise add collapse?

Also quantifies the cactus(long-read,80) vs PG(short-read,151) k-mer-budget
imbalance and the discriminative (low-ac) k-mers the absorbed founders lose to
the ac==1 filter.

Usage: python reproduce_founder_collapse.py
Writes results to founder_collapse_diag.npz in this directory.
DO NOT edit src/kmate/* — this imports solve_em read-only.
"""
import os, sys, json, time
import numpy as np
from scipy.sparse import load_npz

ROOT = "/global/scratch/users/tbellg/kmate"
sys.path.insert(0, f"{ROOT}/src")
from kmate.em_solver import solve_em

OUTDIR = f"{ROOT}/analysis/grenenet_gea/qc/seedmix_validation/diag"
FILT = f"{ROOT}/data/kmer_pa_231_arch3_filt2inv/kmer_pa_Chr1"
RAW  = f"{ROOT}/benchmarks/p231/data/kmer_pa_p231/kmer_pa_Chr1"
SPLIT = json.load(open(f"{ROOT}/data/founder_split_cactus_pg.json"))
CACTUS = set(map(str, SPLIT["cactus"]))   # 80 long-read
PG     = set(map(str, SPLIT["PG"]))       # 151 short-read
WATCH = ["9977", "9985", "10013", "9941", "9507", "9761"]

MAXIT, TOL = 400, 1e-8


def load_panel(prefix):
    t = time.time()
    K = load_npz(f"{prefix}.kmer_pa.npz").tocsr().astype(np.float32)
    m = np.load(f"{prefix}.meta.npz", allow_pickle=True)
    fo = m["founders"].astype(str)
    bid = np.asarray(m["bubble_id"]).astype(np.int64)
    print(f"  loaded {prefix.split('/')[-1]}  K={K.shape} nnz={K.nnz:,} "
          f"({time.time()-t:.0f}s)", flush=True)
    return K, fo, bid


def run_em_ideal(K, omega, tag):
    F = K.shape[0]
    h_true = np.full(F, 1.0 / F, dtype=np.float64)
    c = np.asarray(K.T.dot(h_true)).ravel().astype(np.float32)   # ideal counts = ac/F
    t = time.time()
    h, info = solve_em(c, K, coverage=1.0,
                       h_init=np.full(F, 1.0 / F, dtype=np.float32),
                       max_iter=MAXIT, tol=TOL, omega=omega)
    print(f"  [{tag}] EM {info['iterations']} it, conv={info['converged']} "
          f"eff_n={1/np.sum(h**2):.1f}  ({time.time()-t:.0f}s)", flush=True)
    return h


def run_em_noisy(K, omega, tag, cov=30.0, seed=0):
    """Poisson read noise: expected counts scaled to per-founder coverage cov,
    then Poisson-sampled. Total genome-copy 'reads' ~ cov per founder."""
    F = K.shape[0]
    h_true = np.full(F, 1.0 / F, dtype=np.float64)
    mu = np.asarray(K.T.dot(h_true)).ravel() * (cov * F)   # scale: each founder ~cov
    rng = np.random.default_rng(seed)
    c = rng.poisson(mu).astype(np.float32)
    nz = c > 0
    Kn = K[:, nz]
    om = None if omega is None else omega[nz]
    h, info = solve_em(c[nz], Kn, coverage=cov,
                       h_init=np.full(F, 1.0 / F, dtype=np.float32),
                       max_iter=MAXIT, tol=TOL, omega=om)
    print(f"  [{tag}] noisy EM {info['iterations']} it eff_n={1/np.sum(h**2):.1f}",
          flush=True)
    return h


def summarize(h, fo, tag):
    F = len(fo)
    u = 1.0 / F
    order = np.argsort(h)
    n_lo = int((h < 1e-3).sum())
    is_cac = np.array([f in CACTUS for f in fo])
    mass_cac = float(h[is_cac].sum())
    mass_pg = float(h[~is_cac].sum())
    # expected mass if h were uniform: 80/231 vs 151/231
    exp_cac = is_cac.sum() / F
    print(f"\n=== {tag} ===")
    print(f"  h range [{h.min():.2e}, {h.max():.2e}]  uniform={u:.2e}  "
          f"eff_n={1/np.sum(h**2):.1f}")
    print(f"  founders < 1e-3 : {n_lo}  |  < u/10 : {int((h<u/10).sum())}  "
          f"|  < u/2 : {int((h<u/2).sum())}")
    print(f"  cactus mass={mass_cac:.3f} (exp {exp_cac:.3f}, ratio "
          f"{mass_cac/exp_cac:.2f})  PG mass={mass_pg:.3f} (exp {1-exp_cac:.3f}, "
          f"ratio {mass_pg/(1-exp_cac):.2f})")
    print("  10 lowest-h founders:")
    for j in order[:10]:
        print(f"    {fo[j]:>7} {'LR' if is_cac[j] else 'SR'}  h={h[j]:.2e}  "
              f"(h/u={h[j]/u:.3f})")
    print("  WATCH founders:")
    for w in WATCH:
        if w in fo:
            j = list(fo).index(w)
            print(f"    {w:>7} {'LR' if is_cac[j] else 'SR'}  h={h[j]:.2e}  "
                  f"(h/u={h[j]/u:.3f})")
    return dict(n_lo=n_lo, mass_cac=mass_cac, mass_pg=mass_pg,
                exp_cac=exp_cac, hmin=float(h.min()), hmax=float(h.max()),
                eff_n=float(1/np.sum(h**2)))


def budget(K, fo, bid, tag):
    """Per-founder k-mer budget split by cactus/PG + low-ac discriminative budget."""
    ac = np.asarray(K.sum(axis=0)).ravel()          # carriers/kmer
    Kf = np.asarray(K.sum(axis=1)).ravel()          # kmers/founder
    Kcsr = K.tocsr()
    # per founder counts of ac==1 and ac<=5 kmers
    priv = np.zeros(len(fo)); lowac = np.zeros(len(fo)); meanac = np.zeros(len(fo))
    ac1_mask = (ac == 1); ac5_mask = (ac <= 5)
    for f in range(len(fo)):
        cols = Kcsr[f].indices
        priv[f] = ac1_mask[cols].sum()
        lowac[f] = ac5_mask[cols].sum()
        meanac[f] = ac[cols].mean() if cols.size else 0
    is_cac = np.array([f in CACTUS for f in fo])
    print(f"\n=== BUDGET {tag} ===")
    print(f"  total kmers/founder: cactus median={np.median(Kf[is_cac]):.0f} "
          f"PG median={np.median(Kf[~is_cac]):.0f}  "
          f"(ratio {np.median(Kf[is_cac])/np.median(Kf[~is_cac]):.2f})")
    print(f"  ac==1 (private)/founder: cactus median={np.median(priv[is_cac]):.0f} "
          f"PG median={np.median(priv[~is_cac]):.0f}")
    print(f"  ac<=5 (discriminative)/founder: cactus median={np.median(lowac[is_cac]):.0f} "
          f"PG median={np.median(lowac[~is_cac]):.0f}")
    print(f"  mean_ac/founder: cactus median={np.median(meanac[is_cac]):.1f} "
          f"PG median={np.median(meanac[~is_cac]):.1f}")
    print("  WATCH founder budgets (Kf | ac==1 | ac<=5 | mean_ac):")
    for w in WATCH:
        if w in fo:
            j = list(fo).index(w)
            print(f"    {w:>7} {'LR' if is_cac[j] else 'SR'}  "
                  f"{Kf[j]:>8.0f} | {priv[j]:>7.0f} | {lowac[j]:>7.0f} | {meanac[j]:.1f}")
    return dict(Kf=Kf, priv=priv, lowac=lowac, meanac=meanac, ac_hist_le5=lowac)


def main():
    out = {}
    # ---------- FILT2INV ----------
    Kf_, fo, bidf = load_panel(FILT)
    m_bf = np.bincount(bidf)[bidf].astype(np.float32)
    omf = (1.0 / m_bf).astype(np.float32)
    out["founders"] = fo
    bf = budget(Kf_, fo, bidf, "filt2inv")
    for k in ["Kf", "priv", "lowac", "meanac"]:
        out[f"filt_{k}"] = bf[k]

    h1 = run_em_ideal(Kf_, None, "filt2inv/None")
    out["h_filt_none"] = h1; out["sum_filt_none"] = summarize(h1, fo, "1. filt2inv  omega=None (MLE)")

    h2 = run_em_ideal(Kf_, omf, "filt2inv/omega")
    out["h_filt_omega"] = h2; out["sum_filt_omega"] = summarize(h2, fo, "2. filt2inv  omega=1/m_b (PRODUCTION)")

    hns = run_em_noisy(Kf_, omf, "filt2inv/omega+noise", cov=30.0)
    out["h_filt_omega_noisy"] = hns
    summarize(hns, fo, "(+) filt2inv omega=1/m_b + Poisson noise (cov=30)")

    del Kf_

    # ---------- RAW P231 ----------
    Kr_, fo2, bidr = load_panel(RAW)
    assert (fo == fo2).all(), "founder order mismatch"
    m_br = np.bincount(bidr)[bidr].astype(np.float32)
    omr = (1.0 / m_br).astype(np.float32)
    br = budget(Kr_, fo2, bidr, "raw_p231")
    for k in ["Kf", "priv", "lowac", "meanac"]:
        out[f"raw_{k}"] = br[k]

    h3 = run_em_ideal(Kr_, None, "raw/None")
    out["h_raw_none"] = h3; out["sum_raw_none"] = summarize(h3, fo2, "3. raw p231  omega=None (keeps ac==1)")

    h4 = run_em_ideal(Kr_, omr, "raw/omega")
    out["h_raw_omega"] = h4; out["sum_raw_omega"] = summarize(h4, fo2, "4. raw p231  omega=1/m_b")

    np.savez_compressed(f"{OUTDIR}/founder_collapse_diag.npz",
                        **{k: v for k, v in out.items() if not isinstance(v, dict)},
                        summaries=json.dumps({k: v for k, v in out.items()
                                              if isinstance(v, dict)}))
    print(f"\nsaved -> {OUTDIR}/founder_collapse_diag.npz")


if __name__ == "__main__":
    main()
