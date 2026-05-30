"""
Run pipeline against a known-truth Chr1 simulated pool.

Truth: uniform 1/82 over the 82 panel founders.
Tests:
- Recovery accuracy (Pearson r between h_hat and truth)
- Whether ω alternate-fit improves accuracy
- Per-block solve vs all-bubbles-as-one-block
- Coverage fit accuracy
"""
from __future__ import annotations
import sys, os, time
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import numpy as np
import pandas as pd
from scipy.sparse import load_npz
from block_solver import solve_block_wls, solve_block_irls, update_omega
from kmer_count import count_kmers_in_fasta


DATA = os.path.join(os.path.dirname(__file__), "..", "data")
SIM_PREFIX = os.path.join(DATA, "sim_chr1", "uniform82")


def main():
    print("="*70)
    print("Chr1 simulation end-to-end test")
    print("="*70)

    # Load kmer_pa
    cn_sparse = load_npz(os.path.join(DATA, "test_chr1_first200.kmer_pa.npz"))
    meta = np.load(os.path.join(DATA, "test_chr1_first200.meta.npz"), allow_pickle=True)
    kmer_index = meta["kmer_index"]
    bubble_id = meta["bubble_id"]
    founders = meta["founders"]
    F, K = cn_sparse.shape
    kmer_pa = np.asarray(cn_sparse.todense()).astype(np.int8)
    print(f"\n[1] kmer_pa matrix: F={F} × K={K:,}")
    print(f"    AC distribution: median={int(np.median(kmer_pa.sum(axis=0)))}  AC=1: {(kmer_pa.sum(axis=0)==1).sum():,}")

    # Load truth
    truth = pd.read_csv(f"{SIM_PREFIX}_truth.tsv", sep="\t")
    truth_dict = dict(zip(truth.founder.astype(str), truth.weight))
    h_true = np.array([truth_dict.get(str(f), 0.0) for f in founders])
    h_true = h_true / h_true.sum()
    print(f"\n[2] truth: {(h_true > 0).sum()} founders with weight, "
          f"min={h_true[h_true>0].min():.4f}, max={h_true.max():.4f}")

    # Count k-mers in simulated pool FASTQs
    fq1 = f"{SIM_PREFIX}_pool_1.fq.gz"
    fq2 = f"{SIM_PREFIX}_pool_2.fq.gz"
    print(f"\n[3] Counting k-mers in simulated pool...")
    t0 = time.time()
    counts_dict = count_kmers_in_fasta([fq1, fq2], list(kmer_index),
                                        k=31, threads=8, hash_size="2G")
    counts = np.array([counts_dict[km] for km in kmer_index], dtype=np.int64)
    print(f"    [took {time.time()-t0:.0f}s]")
    print(f"    nonzero k-mers: {(counts>0).sum():,}/{K:,} ({(counts>0).mean():.1%})")
    print(f"    median nonzero count: {np.median(counts[counts>0]):.1f}, "
          f"max: {counts.max()}")

    # Estimate coverage from total counts: total_count = K * lambda_kmer * mean(kmer_pa @ h)
    # Use the AC distribution as the prior: mean kmer_pa^T @ uniform = mean AC / F
    mean_ac = kmer_pa.sum(axis=0).mean()
    expected_per_kmer = mean_ac / F
    cov_kmer_est = counts.mean() / expected_per_kmer
    print(f"    coverage estimate from data: {cov_kmer_est:.1f}×")

    # ============== Variant A: WLS, no ω ==============
    print(f"\n[4A] WLS, ω=0 (baseline)")
    t0 = time.time()
    h_a, obj = solve_block_wls(counts, kmer_pa, coverage=cov_kmer_est)
    print(f"    [took {time.time()-t0:.0f}s, obj={obj:.0f}]")
    r_a = np.corrcoef(h_a, h_true)[0, 1]
    err_a = np.linalg.norm(h_a - h_true)
    print(f"    Pearson r vs truth: {r_a:.3f}")
    print(f"    ||h - h_true|| = {err_a:.4f}")
    print(f"    h sum: {h_a.sum():.4f}, support>0.001: {(h_a>0.001).sum()}/{F}")

    # ============== Variant B: WLS + ω alternate-fit ==============
    print(f"\n[4B] WLS + ω alternate-fit (one sample, multiple iters)")
    # For a single sample we can't really fit ω across-samples, but we can
    # estimate ω as the residual after fitting h. Iterate.
    omega = np.zeros(K)
    h_b = np.full(F, 1.0/F)
    for outer in range(3):
        h_b, _ = solve_block_wls(counts, kmer_pa, coverage=cov_kmer_est, omega=omega)
        # Estimate omega: counts - cov*kmer_pa.T@h, clipped at 0
        residual = counts - cov_kmer_est * (kmer_pa.T @ h_b)
        omega = np.maximum(0, residual / cov_kmer_est)
        # only count nonzero contamination
        n_contam = (omega > 0.5).sum()
        r_b = np.corrcoef(h_b, h_true)[0, 1]
        print(f"    iter {outer}: r={r_b:.3f}, n_contaminated_kmers (ω>0.5)={n_contam:,}")
    print(f"    final r vs truth: {r_b:.3f}, ||h - h_true|| = {np.linalg.norm(h_b - h_true):.4f}")

    # ============== Variant C: scaled coverage fit (joint) ==============
    print(f"\n[4C] try several coverage scales (no ω) — pick best")
    cov_scales = [0.5, 0.75, 1.0, 1.25, 1.5, 2.0]
    best_r = 0
    for scale in cov_scales:
        cov_try = cov_kmer_est * scale
        h_c, obj_c = solve_block_wls(counts, kmer_pa, coverage=cov_try)
        r_c = np.corrcoef(h_c, h_true)[0, 1]
        if r_c > best_r:
            best_r = r_c
            best_cov = cov_try
            h_c_best = h_c.copy()
        print(f"    cov={cov_try:.1f}× (×{scale}): r={r_c:.3f}, obj={obj_c:.0f}")
    print(f"    best: cov={best_cov:.1f}×, r={best_r:.3f}")

    # ============== Variant D: per-block solve and aggregate ==============
    # Group bubbles into 10 sub-blocks (~20 bubbles each).
    print(f"\n[4D] per-block solve (10 sub-blocks of ~20 bubbles each), then average")
    n_blocks = 10
    # Map each kmer to a sub-block based on its bubble_id
    bubbles_per_block = max(1, len(set(bubble_id)) // n_blocks)
    sub_block = bubble_id // bubbles_per_block
    h_d_per_block = []
    for b in range(sub_block.max() + 1):
        mask = sub_block == b
        if mask.sum() < F: continue  # skip blocks with fewer kmers than founders
        h_b_sub, _ = solve_block_wls(counts[mask], kmer_pa[:, mask], coverage=cov_kmer_est)
        h_d_per_block.append(h_b_sub)
    h_d = np.mean(h_d_per_block, axis=0)
    h_d = h_d / h_d.sum()  # renormalize
    r_d = np.corrcoef(h_d, h_true)[0, 1]
    print(f"    n sub-blocks solved: {len(h_d_per_block)}")
    print(f"    averaged r vs truth: {r_d:.3f}, ||h - h_true|| = {np.linalg.norm(h_d - h_true):.4f}")

    # Summary
    print(f"\n{'='*70}\nSUMMARY")
    print(f"{'='*70}")
    print(f"  A  (WLS, ω=0):                r={r_a:.3f}")
    print(f"  B  (WLS, ω alt-fit):          r={r_b:.3f}")
    print(f"  C  (best cov scale):           r={best_r:.3f}  at cov={best_cov:.1f}×")
    print(f"  D  (per-block average):        r={r_d:.3f}")
    print(f"  TRUTH: uniform 1/82 = {1/82:.4f}")


if __name__ == "__main__":
    main()
