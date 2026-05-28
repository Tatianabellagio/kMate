"""
Test the JOINT hapFIRE-style solver against the others on the same battery.

Adds a new column: JOINT — per-bubble haplotype reduction + joint CVXPY.
Hypothesis: should beat both flat-EM and per-bubble-average for non-uniform
truth, because it pools per-bubble evidence into ONE global founder vector.
"""
import sys, os, time
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import numpy as np
from scipy.sparse import load_npz
from block_solver import solve_block_wls
from em_solver import solve_em
from joint_solver import solve_joint


DATA = os.path.join(os.path.dirname(__file__), "..", "data")


def metrics(h_hat, h_true):
    h_hat = np.asarray(h_hat)
    if h_hat.sum() > 0:
        h_hat = h_hat / h_hat.sum()
    rmse = float(np.sqrt(np.mean((h_hat - h_true)**2)))
    if h_true.std() > 1e-9:
        ss_res = ((h_true - h_hat)**2).sum()
        ss_tot = ((h_true - h_true.mean())**2).sum()
        r2 = float(1 - ss_res / ss_tot)
    else:
        r2 = float("nan")
    return r2, rmse


def main():
    cn_sparse = load_npz(os.path.join(DATA, "test_chr1_first200.cn.npz"))
    meta = np.load(os.path.join(DATA, "test_chr1_first200.meta.npz"), allow_pickle=True)
    bubble_id = meta["bubble_id"]
    F, K = cn_sparse.shape
    cn = np.asarray(cn_sparse.todense()).astype(np.int8)
    print(f"cn: F={F}, K={K:,}, n_bubbles={len(set(bubble_id))}")

    rng = np.random.default_rng(42)
    pools = [
        ("UNIFORM", np.full(F, 1.0/F)),
        ("NEAR-UNI α=10", rng.dirichlet(np.full(F, 10.0))),
        ("SEEDMIX-LIKE α=2", rng.dirichlet(np.full(F, 2.0))),
        ("SPARSE α=0.3", rng.dirichlet(np.full(F, 0.3))),
    ]
    h = np.zeros(F); h[[0, 1, 2, 3, 4]] = [0.40, 0.25, 0.15, 0.10, 0.10]
    pools.append(("SKEWED 5", h))

    print(f"\n{'Pool':<18} {'cov':>4}  {'flat-EM':>14}  {'flat-WLS':>14}  {'JOINT':>14}")
    print(f"{'':<18} {'':>4}  {'R²    RMSE':>14}  {'R²    RMSE':>14}  {'R²    RMSE':>14}")
    print("-" * 76)

    coverages = [5, 10, 20, 30]
    for label, h_true in pools:
        for cov in coverages:
            mu = cov * (h_true @ cn)
            counts = rng.poisson(np.maximum(mu, 1e-6))

            h_em, _ = solve_em(counts, cn, cov, max_iter=200)
            try:
                h_wls, _ = solve_block_wls(counts, cn, cov)
            except Exception:
                h_wls = np.full(F, 1.0/F)
            t = time.time()
            h_jt, info = solve_joint(counts, cn, bubble_id)
            t_jt = time.time() - t
            if h_jt is None:
                h_jt = np.full(F, 1.0/F)

            r2_em, rmse_em = metrics(h_em, h_true)
            r2_wls, rmse_wls = metrics(h_wls, h_true)
            r2_jt, rmse_jt = metrics(h_jt, h_true)

            def fmt(r, e):
                rs = f"{r:6.3f}" if not np.isnan(r) else "  n/a "
                return f"{rs} {e:7.4f}"

            print(f"{label:<18} {cov:>4}  {fmt(r2_em, rmse_em):>14}  "
                  f"{fmt(r2_wls, rmse_wls):>14}  {fmt(r2_jt, rmse_jt):>14}")
        print()

    # Show top-5 for SKEWED 5 with JOINT solver
    print(f"{'-'*76}\nDiagnostic: SKEWED 5 at cov=10 — top 8 by JOINT")
    h_true = np.zeros(F); h_true[[0, 1, 2, 3, 4]] = [0.40, 0.25, 0.15, 0.10, 0.10]
    counts = rng.poisson(np.maximum(10 * (h_true @ cn), 1e-6))
    h_jt, info = solve_joint(counts, cn, bubble_id, verbose=True)
    h_jt_norm = h_jt / h_jt.sum()
    for i in np.argsort(-h_jt_norm)[:8]:
        marker = "✓" if h_true[i] > 0 else " "
        print(f"  {marker} f{i}: h_jt={h_jt_norm[i]:.4f}  truth={h_true[i]:.4f}")


if __name__ == "__main__":
    main()
