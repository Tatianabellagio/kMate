"""
Test the two-stage architecture with different block sizes.

Key insight: at small block sizes (1-5 bubbles), the unique-haplotype
reduction is huge (5-15 unique haplotypes vs 82 founders), making each
per-block EM well-conditioned.

For each pool type and coverage:
  1. Partition into blocks of N bubbles
  2. Run two-stage solve per block (reduce → EM → equal-share project)
  3. Average h across blocks
  4. Compare to flat EM and flat WLS.
"""
import sys, os, time
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import numpy as np
from scipy.sparse import load_npz
from block_solver import solve_block_wls
from em_solver import solve_em
from hapfire_solver import (solve_block_two_stage, reduce_to_unique_haplotypes,
                              solve_em_simple, project_haplotypes_to_founders)


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


def two_stage_aggregated(counts, cn, bubble_id, bubbles_per_block, F):
    """Run two-stage per sub-block; return averaged founder freqs."""
    sub_block = bubble_id // bubbles_per_block
    h_per_block = []
    for b in range(int(sub_block.max()) + 1):
        mask = sub_block == b
        if mask.sum() < 3:
            continue
        cn_sub = cn[:, mask]
        c_sub = counts[mask]
        cn_uniq, founder_to_hap, hap_to_founders = reduce_to_unique_haplotypes(cn_sub)
        if cn_uniq.shape[0] < 2:
            continue  # no haplotype variation, skip
        p_hap, _ = solve_em_simple(c_sub, cn_uniq, max_iter=100, tol=1e-7)
        h_b = project_haplotypes_to_founders(p_hap, hap_to_founders, F, "equal_share")
        if h_b.sum() > 0:
            h_per_block.append(h_b / h_b.sum())
    if not h_per_block:
        return np.full(F, 1.0/F)
    return np.mean(h_per_block, axis=0)


def main():
    cn_sparse = load_npz(os.path.join(DATA, "test_chr1_first200.cn.npz"))
    meta = np.load(os.path.join(DATA, "test_chr1_first200.meta.npz"), allow_pickle=True)
    bubble_id = meta["bubble_id"]
    F, K = cn_sparse.shape
    cn = np.asarray(cn_sparse.todense()).astype(np.int8)

    print(f"="*100)
    print(f"Block-size sweep: two-stage architecture with varying bubbles/block")
    print(f"  F={F}, K={K:,}, n_bubbles={len(set(bubble_id))}")
    print(f"="*100)

    rng = np.random.default_rng(42)
    pools = [
        ("UNIFORM", np.full(F, 1.0/F)),
        ("NEAR-UNI α=10", rng.dirichlet(np.full(F, 10.0))),
        ("SEEDMIX-LIKE α=2", rng.dirichlet(np.full(F, 2.0))),
        ("SPARSE α=0.3", rng.dirichlet(np.full(F, 0.3))),
    ]
    h = np.zeros(F)
    h[[0, 1, 2, 3, 4]] = [0.40, 0.25, 0.15, 0.10, 0.10]
    pools.append(("SKEWED 5", h))

    print(f"\n{'Pool':<18} {'cov':>4}  {'flat-EM':>15}  {'2S 1bub':>15}  {'2S 5bub':>15}  {'2S 20bub':>15}  {'flat-WLS':>15}")
    print(f"{'':<18} {'':>4}  {'R²    RMSE':>15}  {'R²    RMSE':>15}  {'R²    RMSE':>15}  {'R²    RMSE':>15}  {'R²    RMSE':>15}")
    print("-" * 110)

    coverages = [5, 10, 20, 30]

    for label, h_true in pools:
        for cov in coverages:
            mu = cov * (h_true @ cn)
            counts = rng.poisson(np.maximum(mu, 1e-6))

            # flat EM
            h_em, _ = solve_em(counts, cn, cov, max_iter=200, tol=1e-7)
            # 2-stage with various block sizes
            h_2s_1 = two_stage_aggregated(counts, cn, bubble_id, 1, F)
            h_2s_5 = two_stage_aggregated(counts, cn, bubble_id, 5, F)
            h_2s_20 = two_stage_aggregated(counts, cn, bubble_id, 20, F)
            # flat WLS
            try:
                h_wls, _ = solve_block_wls(counts, cn, cov)
            except Exception:
                h_wls = np.full(F, 1.0/F)

            r2_em, rmse_em = metrics(h_em, h_true)
            r2_1, rmse_1 = metrics(h_2s_1, h_true)
            r2_5, rmse_5 = metrics(h_2s_5, h_true)
            r2_20, rmse_20 = metrics(h_2s_20, h_true)
            r2_w, rmse_w = metrics(h_wls, h_true)

            def fmt(r, e):
                rs = f"{r:6.3f}" if not np.isnan(r) else "  n/a "
                return f"{rs} {e:7.4f}"

            print(f"{label:<18} {cov:>4}  {fmt(r2_em, rmse_em):>15}  "
                  f"{fmt(r2_1, rmse_1):>15}  {fmt(r2_5, rmse_5):>15}  "
                  f"{fmt(r2_20, rmse_20):>15}  {fmt(r2_w, rmse_w):>15}")
        print()


if __name__ == "__main__":
    main()
