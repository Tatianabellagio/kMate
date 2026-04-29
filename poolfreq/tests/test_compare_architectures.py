"""
Compare all three architectures on real cn matrix at multiple coverages and pool types.

Architectures:
  WLS         baseline weighted least-squares (block_solver.solve_block_wls)
  EM_flat     flat EM directly over F founders (em_solver.solve_em)
  HF_2stage   two-stage hapFIRE-style: reduce to unique haplotypes → EM → project
              (hapfire_solver.solve_block_two_stage)
  HF_2stage+CVXPY same as HF_2stage but project via CVXPY simplex constraint

Metrics: R² and RMSE (no Pearson r per user request).

Tests:
  Pool types: UNIFORM, NEAR-UNIFORM α=10, SEEDMIX-LIKE α=2, SPARSE α=0.3, SKEWED 5
  Coverages: 5, 10, 20, 30 (10× is the GrENE-Net realistic case)
"""
import sys, os, time
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import numpy as np
from scipy.sparse import load_npz

from block_solver import solve_block_wls
from em_solver import solve_em
from hapfire_solver import solve_block_two_stage, reduce_to_unique_haplotypes


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
    F, K = cn_sparse.shape
    cn = np.asarray(cn_sparse.todense()).astype(np.int8)
    ac = cn.sum(axis=0)

    # How many unique haplotype patterns in this cn?
    cn_uniq, f2h, h2f = reduce_to_unique_haplotypes(cn)
    print(f"="*86)
    print(f"Compare architectures on real 200-bubble cn (F={F}, K={K:,}, cov tested 5-30×)")
    print(f"  unique haplotypes in this block: H={cn_uniq.shape[0]} (F={F}, ratio={F/cn_uniq.shape[0]:.1f})")
    print(f"  hap_to_founders sizes: {[len(x) for x in h2f[:10]]}{'...' if len(h2f) > 10 else ''}")
    print(f"="*86)

    rng = np.random.default_rng(42)

    pools = []
    pools.append(("UNIFORM", np.full(F, 1.0/F)))

    h = rng.dirichlet(np.full(F, 10.0))
    pools.append(("NEAR-UNI α=10", h))

    h = rng.dirichlet(np.full(F, 2.0))
    pools.append(("SEEDMIX-LIKE α=2", h))

    h = rng.dirichlet(np.full(F, 0.3))
    pools.append(("SPARSE α=0.3", h))

    h = np.zeros(F)
    h[[0, 1, 2, 3, 4]] = [0.40, 0.25, 0.15, 0.10, 0.10]
    pools.append(("SKEWED 5", h))

    coverages = [5, 10, 20, 30]

    print(f"\n{'Pool':<18}{'cov':>5}  {'WLS':>20}  {'EM_flat':>20}  {'HF_2stage':>20}")
    print(f"{'':<18}{'':>5}  {'R²    RMSE':>20}  {'R²    RMSE':>20}  {'R²    RMSE':>20}")
    print("-" * 92)

    for label, h_true in pools:
        for cov in coverages:
            mu = cov * (h_true @ cn)
            counts = rng.poisson(np.maximum(mu, 1e-6))

            # WLS
            try:
                h_wls, _ = solve_block_wls(counts, cn, cov)
                r2_w, rmse_w = metrics(h_wls, h_true)
            except Exception:
                r2_w, rmse_w = float("nan"), float("nan")

            # EM flat
            h_em, _ = solve_em(counts, cn, cov, max_iter=200, tol=1e-7)
            r2_e, rmse_e = metrics(h_em, h_true)

            # HF two-stage (equal-share projection)
            h_hf = solve_block_two_stage(counts, cn, cov,
                                         project_method="equal_share",
                                         em_max_iter=200)
            r2_h, rmse_h = metrics(h_hf, h_true)

            def fmt(r, e):
                rs = f"{r:6.3f}" if not np.isnan(r) else "  n/a "
                return f"{rs} {e:7.4f}"

            print(f"{label:<18}{cov:>5}  {fmt(r2_w, rmse_w):>20}  {fmt(r2_e, rmse_e):>20}  {fmt(r2_h, rmse_h):>20}")
        print()

    # Quick sanity: the two-stage on a single block can't distinguish founders
    # sharing the same haplotype. Show this:
    print(f"\n{'-'*86}\nDiagnostic: two-stage SKEWED 5 (h_hf vs truth on top founders)")
    h_true = np.zeros(F)
    h_true[[0, 1, 2, 3, 4]] = [0.40, 0.25, 0.15, 0.10, 0.10]
    counts = rng.poisson(np.maximum(10 * (h_true @ cn), 1e-6))
    h_hf, diag = solve_block_two_stage(counts, cn, 10, return_diagnostics=True)
    print(f"  unique haplotypes: H={diag['n_unique_haplotypes']} for F={F}")
    print(f"  truth founders: 0-4")
    print(f"  founder 0 lives in haplotype: {f2h[0]}, sharing with: {h2f[f2h[0]][:10]}")
    # Expected from per-block-equal-share if 5 truth founders are alone in their haplotypes:
    # h_hat[f] = h_true[f]
    print(f"  inferred top 8 founders (h_hf):")
    h_hf_norm = h_hf / h_hf.sum()
    for i in np.argsort(-h_hf_norm)[:8]:
        marker = "✓" if h_true[i] > 0 else " "
        print(f"    {marker} f{i}: h_hf={h_hf_norm[i]:.4f}  truth={h_true[i]:.4f}")


if __name__ == "__main__":
    main()
