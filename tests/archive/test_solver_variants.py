"""
Test multiple solver variants on the existing uniform-82 Chr1 simulation.

Try to find a solver setup that recovers uniform h to ||Δh|| < 0.02.

Variants:
  V0  baseline WLS (current)
  V1  drop AC=1 k-mers
  V2  drop AC<5 k-mers (only common variants)
  V3  L1 loss instead of L2
  V4  Huber loss
  V5  L2 reg toward uniform (Tikhonov toward 1/F)
  V6  weighted L2 with weights = 1/(c + 0.5)  (heavier underweighting of zeros)
  V7  Drop k-mers where observed count = 0 (only "informative" rows)
  V8  Combine: AC>=2 + L2 reg toward uniform
"""
from __future__ import annotations
import sys, os, time
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import numpy as np
import pandas as pd
import cvxpy as cp
from scipy.sparse import load_npz
from kmer_count import count_kmers_in_fasta


def solve_block_custom(counts, kmer_pa, coverage, alpha=1.0, omega=None,
                       weights=None, loss="l2", reg_uniform=0.0,
                       reg_sparsity=0.0):
    """Generic CVXPY block solver supporting several loss + regularization options."""
    K = counts.shape[0]
    F = kmer_pa.shape[0]
    if omega is None: omega = np.zeros(K)
    if np.isscalar(alpha): alpha = np.full(K, alpha)
    if weights is None: weights = 1.0 / np.sqrt(counts + 1.0)

    h = cp.Variable(F, nonneg=True)
    mu = coverage * cp.multiply(alpha, kmer_pa.T @ h) + coverage * omega
    residual = cp.multiply(weights, counts - mu)

    if loss == "l2":
        data_term = cp.sum_squares(residual)
    elif loss == "l1":
        data_term = cp.norm1(residual)
    elif loss == "huber":
        data_term = cp.sum(cp.huber(residual, M=2.0))
    else:
        raise ValueError(loss)

    reg_terms = 0
    if reg_uniform > 0:
        reg_terms = reg_terms + reg_uniform * cp.sum_squares(h - 1.0/F)

    obj = cp.Minimize(data_term + reg_terms)
    prob = cp.Problem(obj, [cp.sum(h) == 1])
    prob.solve(solver="SCS", verbose=False)
    if h.value is None:
        raise RuntimeError(f"solver failed: status={prob.status}")
    return np.asarray(h.value), prob.value


DATA = os.path.join(os.path.dirname(__file__), "..", "data")
SIM_PREFIX = os.path.join(DATA, "sim_chr1", "uniform82")


def main():
    print("="*72)
    print("Solver variants on uniform-82 Chr1 simulation")
    print("="*72)

    # Load kmer_pa
    cn_sparse = load_npz(os.path.join(DATA, "test_chr1_first200.kmer_pa.npz"))
    meta = np.load(os.path.join(DATA, "test_chr1_first200.meta.npz"), allow_pickle=True)
    kmer_index = meta["kmer_index"]
    bubble_id = meta["bubble_id"]
    founders = meta["founders"]
    F, K = cn_sparse.shape
    kmer_pa = np.asarray(cn_sparse.todense()).astype(np.int8)
    ac = kmer_pa.sum(axis=0)
    print(f"\ncn: {F} × {K:,}")
    print(f"AC dist: AC=1: {(ac==1).sum():,}   AC 2-4: {((ac>=2)&(ac<=4)).sum():,}   "
          f"AC 5-10: {((ac>=5)&(ac<=10)).sum():,}   AC>10: {(ac>10).sum():,}")

    # Truth: uniform 1/82
    truth = pd.read_csv(f"{SIM_PREFIX}_truth.tsv", sep="\t")
    truth_dict = dict(zip(truth.founder.astype(str), truth.weight))
    h_true = np.array([truth_dict.get(str(f), 0.0) for f in founders])
    h_true = h_true / h_true.sum()
    print(f"truth: uniform 1/{F} = {1/F:.4f}")

    # Load or compute counts
    counts_path = os.path.join(DATA, "sim_chr1", "uniform82_counts.npz")
    if os.path.exists(counts_path):
        counts = np.load(counts_path)["counts"]
        print(f"\nLoaded cached counts from {counts_path}")
    else:
        print("\nCounting k-mers (one-time)...")
        t0 = time.time()
        fq1 = f"{SIM_PREFIX}_pool_1.fq.gz"
        fq2 = f"{SIM_PREFIX}_pool_2.fq.gz"
        counts_dict = count_kmers_in_fasta([fq1, fq2], list(kmer_index),
                                            k=31, threads=8, hash_size="2G")
        counts = np.array([counts_dict[km] for km in kmer_index], dtype=np.int64)
        np.savez(counts_path, counts=counts)
        print(f"  [took {time.time()-t0:.0f}s]  saved to {counts_path}")

    print(f"counts: nonzero={int((counts>0).sum()):,}/{K:,} ({(counts>0).mean():.1%})  "
          f"median nonzero={float(np.median(counts[counts>0])):.1f}")

    # Coverage estimate from data: total_counts ≈ λ × Σ_k AC[k]/F
    # → λ = total_counts × F / Σ AC[k]
    cov_est = counts.sum() * F / ac.sum()
    print(f"coverage estimate: {cov_est:.1f}×")

    print(f"\n{'-'*72}\nVARIANTS\n{'-'*72}")

    def report(label, h_hat):
        h_hat = np.asarray(h_hat)
        h_hat = h_hat / h_hat.sum()
        err = np.linalg.norm(h_hat - h_true)
        # CV across founders (should be low for uniform recovery)
        cv = h_hat.std() / h_hat.mean()
        n_zero = (h_hat < 1e-4).sum()
        print(f"  {label:55s}  ||Δh||={err:.4f}  CV={cv:.3f}  n_zero={n_zero}/{F}")
        return err

    # V0 baseline
    h0, _ = solve_block_custom(counts, kmer_pa, cov_est)
    report("V0  baseline WLS", h0)

    # V1 drop AC=1
    keep = ac >= 2
    h1, _ = solve_block_custom(counts[keep], kmer_pa[:, keep], cov_est)
    report(f"V1  drop AC=1 (keep {keep.sum():,})", h1)

    # V2 drop AC<5
    keep5 = ac >= 5
    h2, _ = solve_block_custom(counts[keep5], kmer_pa[:, keep5], cov_est)
    report(f"V2  drop AC<5 (keep {keep5.sum():,})", h2)

    # V3 L1 (slow)
    h3, _ = solve_block_custom(counts, kmer_pa, cov_est, loss="l1")
    report("V3  L1 loss", h3)

    # V4 Huber
    h4, _ = solve_block_custom(counts, kmer_pa, cov_est, loss="huber")
    report("V4  Huber loss", h4)

    # V5 L2 + regularization toward uniform (try several lambdas)
    for lam in [0.01, 0.1, 1.0, 10.0]:
        hr, _ = solve_block_custom(counts, kmer_pa, cov_est, reg_uniform=lam)
        report(f"V5  L2 + uniform-reg λ={lam}", hr)

    # V6 different weight scheme: w = 1/(count + 0.5), lighter on zeros
    weights6 = 1.0 / (counts + 0.5)
    h6, _ = solve_block_custom(counts, kmer_pa, cov_est, weights=weights6)
    report("V6  weights = 1/(count+0.5)", h6)

    # V7 drop k-mers with count=0
    nz = counts > 0
    h7, _ = solve_block_custom(counts[nz], kmer_pa[:, nz], cov_est)
    report(f"V7  drop count=0 rows (keep {nz.sum():,})", h7)

    # V8 AC≥2 + uniform reg
    for lam in [0.01, 0.1, 1.0]:
        keep = ac >= 2
        h8, _ = solve_block_custom(counts[keep], kmer_pa[:, keep], cov_est, reg_uniform=lam)
        report(f"V8  AC≥2 + uniform-reg λ={lam}", h8)

    # V9 AC≥5 + Huber
    keep5 = ac >= 5
    h9, _ = solve_block_custom(counts[keep5], kmer_pa[:, keep5], cov_est, loss="huber")
    report(f"V9  AC≥5 + Huber", h9)

    # V10 AC≥2 + drop count=0 + uniform reg
    for lam in [0.01, 0.1, 1.0]:
        keep = (ac >= 2) & (counts > 0)
        if keep.sum() < F: continue
        h10, _ = solve_block_custom(counts[keep], kmer_pa[:, keep], cov_est, reg_uniform=lam)
        report(f"V10 AC≥2 ∧ count>0 + uniform-reg λ={lam} (n={keep.sum():,})", h10)


if __name__ == "__main__":
    main()
