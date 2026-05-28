"""
Run pipeline on the skewed-5 simulated pool, with multiple solver variants.
Truth: 5 founders at [0.40, 0.25, 0.15, 0.10, 0.10].
"""
from __future__ import annotations
import sys, os, time
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import numpy as np
import pandas as pd
import cvxpy as cp
from scipy.sparse import load_npz
from kmer_count import count_kmers_in_fasta


DATA = os.path.join(os.path.dirname(__file__), "..", "data")
SIM_PREFIX = os.path.join(DATA, "sim_chr1_skewed", "skewed5")


def solve_wls(counts, cn, coverage, weights=None):
    K, F = counts.shape[0], cn.shape[0]
    if weights is None: weights = 1.0 / np.sqrt(counts + 1.0)
    h = cp.Variable(F, nonneg=True)
    mu = coverage * (cn.T @ h)
    residual = cp.multiply(weights, counts - mu)
    prob = cp.Problem(cp.Minimize(cp.sum_squares(residual)), [cp.sum(h) == 1])
    prob.solve(solver="SCS", verbose=False)
    return np.asarray(h.value)


def solve_kl(counts, cn, coverage):
    K, F = counts.shape[0], cn.shape[0]
    h = cp.Variable(F, nonneg=True)
    mu = coverage * (cn.T @ h) + 1e-3
    obj = cp.Minimize(cp.sum(cp.kl_div(counts, mu)))
    prob = cp.Problem(obj, [cp.sum(h) == 1])
    prob.solve(solver="SCS", verbose=False)
    return np.asarray(h.value)


def per_block_avg(counts, cn, coverage, bubble_id, n_blocks=20, solver=solve_wls):
    """Solve each sub-block independently, return averaged h."""
    F = cn.shape[0]
    bubbles_per_block = max(1, len(set(bubble_id)) // n_blocks)
    sub_block = bubble_id // bubbles_per_block
    h_acc = []
    for b in range(sub_block.max() + 1):
        mask = sub_block == b
        if mask.sum() < F: continue
        h = solver(counts[mask], cn[:, mask], coverage)
        h_acc.append(h / h.sum())
    return np.mean(h_acc, axis=0)


def main():
    cn_sparse = load_npz(os.path.join(DATA, "test_chr1_first200.cn.npz"))
    meta = np.load(os.path.join(DATA, "test_chr1_first200.meta.npz"), allow_pickle=True)
    kmer_index = meta["kmer_index"]
    bubble_id = meta["bubble_id"]
    founders = meta["founders"]
    F, K = cn_sparse.shape
    cn = np.asarray(cn_sparse.todense()).astype(np.int8)
    ac = cn.sum(axis=0)
    print(f"cn: F={F} × K={K:,}")

    # Truth
    truth = pd.read_csv(f"{SIM_PREFIX}_truth.tsv", sep="\t")
    truth_dict = dict(zip(truth.founder.astype(str), truth.weight))
    h_true = np.array([truth_dict.get(str(f), 0.0) for f in founders])
    h_true = h_true / h_true.sum()
    n_truth = (h_true > 0).sum()
    print(f"truth: {n_truth} non-zero founders")
    for i in np.argsort(-h_true)[:5]:
        print(f"  {founders[i]}: w={h_true[i]:.3f}")

    # Counts
    counts_path = os.path.join(DATA, "sim_chr1_skewed", "skewed5_counts.npz")
    if os.path.exists(counts_path):
        counts = np.load(counts_path)["counts"]
    else:
        print(f"\nCounting k-mers...")
        t = time.time()
        counts_dict = count_kmers_in_fasta(
            [f"{SIM_PREFIX}_pool_1.fq.gz", f"{SIM_PREFIX}_pool_2.fq.gz"],
            list(kmer_index), k=31, threads=8, hash_size="2G")
        counts = np.array([counts_dict[km] for km in kmer_index], dtype=np.int64)
        np.savez(counts_path, counts=counts)
        print(f"  [{time.time()-t:.0f}s]")
    cov = counts.sum() * F / ac.sum()
    print(f"\ncoverage est: {cov:.1f}×, nonzero kmers: {(counts>0).sum():,}/{K:,}")

    def report(label, h_hat):
        h_hat = np.asarray(h_hat) / np.asarray(h_hat).sum()
        err = np.linalg.norm(h_hat - h_true)
        r = np.corrcoef(h_hat, h_true)[0, 1]
        # which top-5 inferred match true?
        top5 = set(np.argsort(-h_hat)[:5])
        truth5 = set(np.argsort(-h_true)[:5])
        recovered = len(top5 & truth5)
        print(f"  {label:50s}  ||Δh||={err:.4f}  r={r:.3f}  top-5 recovered: {recovered}/5")
        return err

    print(f"\n{'-'*72}\nSolver variants on skewed truth\n{'-'*72}")

    # S1 baseline WLS
    h = solve_wls(counts, cn, cov)
    report("S1  baseline WLS", h)

    # S2 Poisson NLL
    h = solve_kl(counts, cn, cov)
    report("S2  Poisson NLL (KL div)", h)

    # S3 per-block average (10 sub-blocks)
    h = per_block_avg(counts, cn, cov, bubble_id, n_blocks=10, solver=solve_wls)
    report("S3  per-block-avg WLS, 10 sub-blocks", h)

    # S4 per-block average (20 sub-blocks)
    h = per_block_avg(counts, cn, cov, bubble_id, n_blocks=20, solver=solve_wls)
    report("S4  per-block-avg WLS, 20 sub-blocks", h)

    # S5 per-block average (50 sub-blocks)
    h = per_block_avg(counts, cn, cov, bubble_id, n_blocks=50, solver=solve_wls)
    report("S5  per-block-avg WLS, 50 sub-blocks", h)

    # S6 per-block KL
    h = per_block_avg(counts, cn, cov, bubble_id, n_blocks=20, solver=solve_kl)
    report("S6  per-block-avg KL, 20 sub-blocks", h)

    # Show top recovered for best
    print(f"\nTop 8 by S5 (per-block-avg WLS, 50 sub-blocks):")
    h_best = per_block_avg(counts, cn, cov, bubble_id, n_blocks=50, solver=solve_wls)
    h_best = h_best / h_best.sum()
    for i in np.argsort(-h_best)[:8]:
        marker = "✓" if h_true[i] > 0 else " "
        print(f"  {marker} {founders[i]:>10s}: h={h_best[i]:.4f}  truth={h_true[i]:.4f}")


if __name__ == "__main__":
    main()
