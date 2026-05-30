"""
End-to-end test: real kmer_pa matrix from build_kmer_pa + synthetic Poisson counts
+ block solver. Verify the solver recovers a known h on real-shaped data.
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import numpy as np
from scipy.sparse import load_npz
from block_solver import solve_block_wls, solve_block_irls


def main():
    DATA = os.path.join(os.path.dirname(__file__), "..", "data")
    cn_sparse = load_npz(os.path.join(DATA, "test_chr1_first200.kmer_pa.npz"))
    meta = np.load(os.path.join(DATA, "test_chr1_first200.meta.npz"), allow_pickle=True)
    founders = meta["founders"]
    bubble_id = meta["bubble_id"]
    F, K = cn_sparse.shape
    print(f"Loaded kmer_pa: {F} founders × {K:,} k-mers")
    print(f"  kmer_pa density: {cn_sparse.nnz / (F*K):.2%}")

    # AC distribution per kmer (how many founders carry each)
    ac_per_kmer = np.asarray(cn_sparse.sum(axis=0)).ravel()
    print(f"  AC distribution: median={np.median(ac_per_kmer):.0f}  "
          f"min={ac_per_kmer.min()}  max={ac_per_kmer.max()}")

    # Take all 10 bubbles together as a single test "block"
    kmer_pa = np.asarray(cn_sparse.todense()).astype(np.int8)

    # Set a sparse h with a few founders dominating
    rng = np.random.default_rng(42)
    h_true = np.zeros(F)
    carriers = rng.choice(F, 5, replace=False)
    h_true[carriers] = rng.dirichlet(np.ones(5))
    print(f"\nTrue h: {h_true.round(3)}")
    print(f"  carriers: {carriers}, weights: {h_true[carriers].round(3)}")

    # Simulate Poisson counts at coverage 50
    lam = 50.0
    mu_true = lam * (kmer_pa.T @ h_true)
    counts = rng.poisson(np.maximum(mu_true, 1e-6))
    print(f"\nSimulated counts: total={counts.sum():,}  nonzero={(counts>0).sum()}/{K}")
    print(f"  median count: {np.median(counts):.0f}  max: {counts.max()}")

    # Solve via WLS
    print("\n--- WLS solve ---")
    h_wls, obj_wls = solve_block_wls(counts, kmer_pa, coverage=lam)
    print(f"  ||h_wls - h_true||: {np.linalg.norm(h_wls - h_true):.4f}")
    print(f"  max |Δh|:          {np.abs(h_wls - h_true).max():.4f}")

    # Top 5 inferred founders
    top_wls = np.argsort(-h_wls)[:5]
    print(f"\nTop 5 inferred (WLS):")
    for f in top_wls:
        print(f"  {founders[f]}: h_hat={h_wls[f]:.3f}, h_true={h_true[f]:.3f}")

    # Solve via IRLS
    print("\n--- IRLS solve ---")
    h_irls, obj_irls = solve_block_irls(counts, kmer_pa, coverage=lam, max_iter=8, verbose=False)
    print(f"  ||h_irls - h_true||: {np.linalg.norm(h_irls - h_true):.4f}")

    # Spot check: AC per kmer used for recovery
    # Ground truth: predicted counts should match observed
    pred = lam * (kmer_pa.T @ h_wls)
    rmse = np.sqrt(np.mean((counts - pred)**2))
    print(f"\nResidual RMSE (counts - λ·kmer_pa.T·h_wls): {rmse:.2f}  (vs sqrt(λ)={np.sqrt(lam):.2f} expected from Poisson)")


if __name__ == "__main__":
    main()
