"""Controlled synthetic before/after test of commit 9669be7 (full-panel Kf_w fix),
run through the REAL production solver (src/kmate/em_solver.py), not a prototype.

Both "OLD" and "NEW" call the SAME current solve_em(); the only difference is how
Kf_w is supplied, replicating exactly what changed in per_sample_per_chrom.py:
  OLD (pre-fix): kmer_pa pre-sliced to nz columns, kfw=None
                 -> solve_em falls back to Kf_w = kmer_pa_nz @ w, i.e. Kf_w summed
                    over the OBSERVED (c_k>0) set only (the survivorship bias).
  NEW (fix):     kfw = kfw_full, computed on the FULL panel (incl. c_k=0 columns)
                 BEFORE the nz slice -> exactly what per_sample_per_chrom.py now does.

Truth = uniform h=1/231 (equimolar seed-mix design). Two depths (T, 0.3T) x
{noiseless, noisy} on the production filt2inv Chr1 panel, current production
config: normalize="per_founder", omega=None (--kmer-weight uniform, the current
global-mode default).
"""
import json, sys, time
import numpy as np
from scipy.sparse import load_npz

ROOT = "/global/scratch/users/tbellg/kmate"
OUT = f"{ROOT}/analysis/grenenet_gea/qc/seedmix_validation/fix_kfw_fullpanel"
PANEL = f"{ROOT}/data/kmer_pa_231_arch3_filt2inv/kmer_pa_Chr1"
sys.path.insert(0, f"{ROOT}/src")
from kmate.em_solver import solve_em

T_FULL = 77.0  # "realistic seedmix-like pool depth" (matches fix_norm/run_experiment.py)


def solve_old(counts, K_full):
    nz = counts > 0
    h, info = solve_em(counts[nz], K_full[:, nz], coverage=1.0,
                        max_iter=400, tol=1e-7, normalize="per_founder", kfw=None)
    return h, info


def solve_new(counts, K_full):
    kfw_full = K_full.sum(axis=1).astype(np.float32)
    nz = counts > 0
    h, info = solve_em(counts[nz], K_full[:, nz], coverage=1.0,
                        max_iter=400, tol=1e-7, normalize="per_founder", kfw=kfw_full)
    return h, info


def main():
    t0 = time.time()
    print("loading Chr1 filt2inv panel ...", flush=True)
    Ksp = load_npz(f"{PANEL}.kmer_pa.npz").astype(np.float32).tocsr()
    K = np.asarray(Ksp.todense(), dtype=np.float32)
    F, Kn = K.shape
    print(f"F={F} K={Kn:,} dense {K.nbytes/1e9:.1f}GB [{time.time()-t0:.0f}s]", flush=True)

    h_true = np.full(F, 1.0 / F, dtype=np.float64)
    mu = h_true.astype(np.float32) @ K   # per-kmer expected rate at uniform truth

    scenarios = {
        "1x_noiseless": (T_FULL, None),
        "1x_noisy": (T_FULL, 0),
        "0.3x_noiseless": (0.3 * T_FULL, None),
        "0.3x_noisy": (0.3 * T_FULL, 1),
    }

    results = {}
    hsave = {}
    for name, (T, seed) in scenarios.items():
        rate = T * mu
        if seed is None:
            counts = rate.astype(np.float32)
        else:
            rng = np.random.default_rng(seed)
            counts = rng.poisson(rate).astype(np.float32)

        for tag, fn in [("old_observed_only", solve_old), ("new_full_panel", solve_new)]:
            t = time.time()
            h, info = fn(counts, K)
            rmse = float(np.sqrt(np.mean((h - h_true) ** 2)))
            nabs = int((h < 1e-3).sum())
            key = f"{name}__{tag}"
            results[key] = dict(rmse=rmse, n_abs=nabs, iters=info["iterations"],
                                 conv=info["converged"], secs=round(time.time() - t, 1))
            hsave[f"h__{key}"] = h
            print(f"[{key:32s}] rmse={rmse:.2e} n_abs={nabs:3d} "
                  f"{info['iterations']}it {results[key]['secs']}s", flush=True)

    np.savez_compressed(f"{OUT}/synthetic_h.npz", founders=np.arange(F), h_true=h_true, **hsave)
    json.dump(results, open(f"{OUT}/synthetic_summary.json", "w"), indent=2)
    print(f"\nsaved -> synthetic_h.npz + synthetic_summary.json ({time.time()-t0:.0f}s total)", flush=True)


if __name__ == "__main__":
    main()
