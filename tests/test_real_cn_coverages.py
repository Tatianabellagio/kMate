"""
Test EM and WLS on the REAL 200-bubble kmer_pa matrix, with simulated counts at
multiple coverages (5×, 10×, 20×, 30×).

Uses the existing uniform82 kmer_pa and skewed5 truths but generates new Poisson
counts at the requested coverage so we can compare across coverages.

Reports R² and RMSE.
"""
import sys, os, time
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import numpy as np
from scipy.sparse import load_npz
from em_solver import solve_em
from block_solver import solve_block_wls


DATA = os.path.join(os.path.dirname(__file__), "..", "data")


def metrics(h_hat, h_true):
    h_hat = h_hat / h_hat.sum() if h_hat.sum() > 0 else h_hat
    rmse = float(np.sqrt(np.mean((h_hat - h_true)**2)))
    if h_true.std() > 1e-9:
        ss_res = ((h_true - h_hat)**2).sum()
        ss_tot = ((h_true - h_true.mean())**2).sum()
        r2 = float(1 - ss_res / ss_tot)
    else:
        r2 = float("nan")
    return r2, rmse


def main():
    # Load real 200-bubble kmer_pa
    cn_sparse = load_npz(os.path.join(DATA, "test_chr1_first200.kmer_pa.npz"))
    meta = np.load(os.path.join(DATA, "test_chr1_first200.meta.npz"), allow_pickle=True)
    founders = meta["founders"]
    F, K = cn_sparse.shape
    kmer_pa = np.asarray(cn_sparse.todense()).astype(np.int8)
    ac = kmer_pa.sum(axis=0)
    print(f"Real kmer_pa: F={F}, K={K:,}")
    print(f"  AC dist: AC=1: {(ac==1).sum():,}  AC 2-4: {((ac>=2)&(ac<=4)).sum():,}  "
          f"AC 5-10: {((ac>=5)&(ac<=10)).sum():,}  AC>10: {(ac>10).sum():,}")

    # Define test pools
    rng = np.random.default_rng(42)

    pools = []

    # Pool 1: UNIFORM 1/82
    h_true = np.full(F, 1.0/F)
    pools.append(("UNIFORM 82", h_true))

    # Pool 2: SKEWED 5 carriers
    h_true = np.zeros(F)
    h_true[[0, 1, 2, 3, 4]] = [0.40, 0.25, 0.15, 0.10, 0.10]
    pools.append(("SKEWED 5", h_true))

    # Pool 3: NEAR-UNIFORM (Dirichlet α=10, mild variation)
    h_true = rng.dirichlet(np.full(F, 10.0))
    pools.append(("NEAR-UNIFORM α=10", h_true))

    # Pool 4: SEEDMIX-like (Dirichlet α=2, max/min ~3-5x)
    h_true = rng.dirichlet(np.full(F, 2.0))
    pools.append(("SEEDMIX-LIKE α=2", h_true))

    # Pool 5: SPARSE (Dirichlet α=0.3)
    h_true = rng.dirichlet(np.full(F, 0.3))
    pools.append(("SPARSE α=0.3", h_true))

    # Coverages to test
    coverages = [5, 10, 20, 30]

    # Header
    print(f"\n{'Pool':<22} {'cov':>5} {'method':<6} {'R²':>10} {'RMSE':>10} {'iter':>6}  pool stats")
    print("-" * 92)

    for label, h_true in pools:
        # Pool stats
        h_max = h_true.max()
        h_min = h_true[h_true > 1e-10].min() if (h_true > 1e-10).any() else 0
        ratio = h_max / h_min if h_min > 0 else float("inf")
        eff_n = 1 / np.sum(h_true**2)

        # For each coverage, generate counts and solve
        for cov in coverages:
            mu = cov * (h_true @ kmer_pa)
            counts = rng.poisson(np.maximum(mu, 1e-6))

            # EM
            t = time.time()
            h_em, info = solve_em(counts, kmer_pa, cov, max_iter=200, tol=1e-7)
            t_em = time.time() - t
            r2_em, rmse_em = metrics(h_em, h_true)

            # WLS
            t = time.time()
            h_wls, _ = solve_block_wls(counts, kmer_pa, cov)
            t_wls = time.time() - t
            r2_wls, rmse_wls = metrics(h_wls, h_true)

            stats = f"  ratio={ratio:.1f} eff_n={eff_n:.1f}"
            r2_em_s = f"{r2_em:.3f}" if not np.isnan(r2_em) else "  n/a"
            r2_wls_s = f"{r2_wls:.3f}" if not np.isnan(r2_wls) else "  n/a"
            print(f"{label:<22} {cov:>5}  EM     {r2_em_s:>10} {rmse_em:>10.4f} {info['iterations']:>6}{stats if cov==coverages[0] else ''}")
            print(f"{label:<22} {cov:>5}  WLS    {r2_wls_s:>10} {rmse_wls:>10.4f} {'':>6}")
        print()


if __name__ == "__main__":
    main()
