"""
Quick unit test of the EM solver on synthetic data.
"""
import sys, os, time
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import numpy as np
from em_solver import solve_em
from block_solver import solve_block_wls


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


def gen_synth(F, K, h_true, coverage, seed=0):
    rng = np.random.default_rng(seed)
    af = rng.uniform(0.1, 0.5, K)
    cn = (rng.random((F, K)) < af[None, :]).astype(np.int8)
    mu = coverage * (h_true @ cn)
    counts = rng.poisson(np.maximum(mu, 1e-6))
    return cn, counts


def main():
    print("="*72)
    print("EM solver quick tests")
    print("="*72)

    rng = np.random.default_rng(0)

    print("\n--- T1: F=20, K=500, cov=50, sparse h ---")
    F, K, cov = 20, 500, 50
    h_true = rng.dirichlet(np.full(F, 0.3))
    cn, counts = gen_synth(F, K, h_true, cov, seed=1)
    t = time.time()
    h_em, info = solve_em(counts, cn, cov, max_iter=200)
    r2, rmse = metrics(h_em, h_true)
    print(f"  EM:  r²={r2:.4f}  rmse={rmse:.4f}  iters={info['iterations']}  [{time.time()-t:.2f}s]")
    h_wls, _ = solve_block_wls(counts, cn, cov)
    r2, rmse = metrics(h_wls, h_true)
    print(f"  WLS: r²={r2:.4f}  rmse={rmse:.4f}")

    print("\n--- T2: F=82, K=18523, cov=18.6, UNIFORM truth ---")
    F = 82
    h_true = np.full(F, 1.0/F)
    cn, counts = gen_synth(F, 18523, h_true, 18.6, seed=2)
    t = time.time()
    h_em, info = solve_em(counts, cn, 18.6, max_iter=200)
    r2, rmse = metrics(h_em, h_true)
    print(f"  EM:  r²={r2:.4f}  rmse={rmse:.4f}  iters={info['iterations']}  [{time.time()-t:.2f}s]")
    h_wls, _ = solve_block_wls(counts, cn, 18.6)
    r2, rmse = metrics(h_wls, h_true)
    print(f"  WLS: r²={r2:.4f}  rmse={rmse:.4f}")

    print("\n--- T3: F=82, K=18523, cov=10, UNIFORM truth (LOW COV) ---")
    cn, counts = gen_synth(F, 18523, h_true, 10, seed=3)
    t = time.time()
    h_em, info = solve_em(counts, cn, 10, max_iter=200)
    r2, rmse = metrics(h_em, h_true)
    print(f"  EM:  r²={r2:.4f}  rmse={rmse:.4f}  iters={info['iterations']}  [{time.time()-t:.2f}s]")
    h_wls, _ = solve_block_wls(counts, cn, 10)
    r2, rmse = metrics(h_wls, h_true)
    print(f"  WLS: r²={r2:.4f}  rmse={rmse:.4f}")

    print("\n--- T4: F=82, K=18523, cov=10, SKEWED truth (5 carriers) ---")
    h_true = np.zeros(F)
    h_true[[2, 13, 27, 41, 55]] = [0.40, 0.25, 0.15, 0.10, 0.10]
    cn, counts = gen_synth(F, 18523, h_true, 10, seed=4)
    t = time.time()
    h_em, info = solve_em(counts, cn, 10, max_iter=200)
    r2, rmse = metrics(h_em, h_true)
    print(f"  EM:  r²={r2:.4f}  rmse={rmse:.4f}  iters={info['iterations']}  [{time.time()-t:.2f}s]")
    h_wls, _ = solve_block_wls(counts, cn, 10)
    r2, rmse = metrics(h_wls, h_true)
    print(f"  WLS: r²={r2:.4f}  rmse={rmse:.4f}")

    print("\n--- T5: F=231, K=50000, cov=10, near-uniform but with structure ---")
    F = 231
    rng = np.random.default_rng(5)
    h_true = rng.dirichlet(np.full(F, 5.0))  # near-uniform with mild variation
    h_true = h_true / h_true.sum()
    cn, counts = gen_synth(F, 50000, h_true, 10, seed=5)
    t = time.time()
    h_em, info = solve_em(counts, cn, 10, max_iter=200)
    r2, rmse = metrics(h_em, h_true)
    print(f"  EM:  r²={r2:.4f}  rmse={rmse:.4f}  iters={info['iterations']}  [{time.time()-t:.2f}s]")
    print(f"       h_true range: [{h_true.min():.4f}, {h_true.max():.4f}]")
    print(f"       h_em   range: [{h_em.min():.4f}, {h_em.max():.4f}]")


if __name__ == "__main__":
    main()
