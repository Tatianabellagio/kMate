"""
Deep diagnostic of why the solver doesn't recover uniform from uniform sim.

Diagnostics:
  D1  kmer_pa matrix sanity: per-founder total k-mer count (Σ_k kmer_pa[f,k])
       — if very heterogeneous, "data poverty" hypothesis
  D2  expected vs observed counts under uniform h: are observed counts
       systematically biased?
  D3  ONE-FOUNDER ground-truth check: simulate reads from one founder only,
       verify h_hat concentrates on that founder
  D4  Try unweighted L2 (w=1)
  D5  Try IRLS (proper Poisson weights)
  D6  Try weighted by 1/sqrt(AC) (down-weight high-AC k-mers)
  D7  Try simple linear regression on the predicted-vs-observed plot
"""
from __future__ import annotations
import sys, os, time
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import numpy as np
import pandas as pd
import cvxpy as cp
from scipy.sparse import load_npz
from kmer_count import count_kmers_in_fasta, count_kmers_in_bam


DATA = os.path.join(os.path.dirname(__file__), "..", "data")
SIM_PREFIX = os.path.join(DATA, "sim_chr1", "uniform82")


def solve_custom(counts, kmer_pa, coverage, weights=None, omega=None):
    K = counts.shape[0]
    F = kmer_pa.shape[0]
    if omega is None: omega = np.zeros(K)
    if weights is None: weights = np.ones(K)
    h = cp.Variable(F, nonneg=True)
    mu = coverage * (kmer_pa.T @ h) + coverage * omega
    residual = cp.multiply(weights, counts - mu)
    obj = cp.Minimize(cp.sum_squares(residual))
    prob = cp.Problem(obj, [cp.sum(h) == 1])
    prob.solve(solver="SCS", verbose=False)
    return np.asarray(h.value)


def main():
    print("="*72)
    print("Solver diagnostics")
    print("="*72)

    # Load
    cn_sparse = load_npz(os.path.join(DATA, "test_chr1_first200.kmer_pa.npz"))
    meta = np.load(os.path.join(DATA, "test_chr1_first200.meta.npz"), allow_pickle=True)
    kmer_index = meta["kmer_index"]
    founders = meta["founders"]
    F, K = cn_sparse.shape
    kmer_pa = np.asarray(cn_sparse.todense()).astype(np.int8)
    ac = kmer_pa.sum(axis=0)
    counts = np.load(os.path.join(DATA, "sim_chr1", "uniform82_counts.npz"))["counts"]
    cov = counts.sum() * F / ac.sum()
    print(f"\ncn: {F} × {K:,}, cov={cov:.1f}×")

    h_true = np.full(F, 1.0/F)

    # ============ D1: per-founder total k-mer count ============
    print(f"\n--- D1: per-founder k-mer totals (Σ_k kmer_pa[f,k]) ---")
    per_founder = kmer_pa.sum(axis=1)
    print(f"  min: {per_founder.min()}  max: {per_founder.max()}  median: {int(np.median(per_founder))}")
    print(f"  CV (std/mean): {per_founder.std()/per_founder.mean():.3f}")
    print(f"  founders with fewest k-mers (n=5): {founders[np.argsort(per_founder)[:5]]} → {sorted(per_founder)[:5]}")
    print(f"  founders with most k-mers (n=5):   {founders[np.argsort(-per_founder)[:5]]} → {sorted(per_founder, reverse=True)[:5]}")

    # ============ D2: expected vs observed under uniform h ============
    print(f"\n--- D2: expected vs observed under uniform h ---")
    expected_uniform = cov * (kmer_pa.T @ h_true)  # K-vector
    # Group counts by AC and compare
    for ac_lo, ac_hi, name in [(1, 1, "AC=1"), (2, 4, "AC 2-4"),
                                (5, 10, "AC 5-10"), (11, 1000, "AC>10")]:
        mask = (ac >= ac_lo) & (ac <= ac_hi)
        if mask.sum() == 0: continue
        obs = counts[mask].sum()
        exp_ = expected_uniform[mask].sum()
        print(f"  {name:8s} (n={mask.sum():>6,})  "
              f"sum_obs={obs:>9,}  sum_exp={exp_:>10.0f}  "
              f"obs/exp={obs/max(1,exp_):.3f}")

    # ============ D3: ONE-FOUNDER simulation read check ============
    # Use the per-founder fastqs we already wrote (or re-simulate one)
    # Per-founder FASTQs were deleted after combining. So skip this for now.
    # (Could add later: re-sim ONE founder, count, check h.)
    print(f"\n--- D3: SKIP one-founder check (per-founder fastqs deleted) ---")

    # ============ D4-D6: solver variants ============
    print(f"\n--- D4-D7: solver variants ---")

    def report(label, h_hat):
        h_hat = np.asarray(h_hat) / np.asarray(h_hat).sum()
        err = np.linalg.norm(h_hat - h_true)
        cv_h = h_hat.std() / h_hat.mean()
        nz = (h_hat < 1e-4).sum()
        print(f"  {label:50s}  ||Δh||={err:.4f}  CV={cv_h:.3f}  n_zero={nz}")
        return err

    # D4 unweighted L2
    h4 = solve_custom(counts, kmer_pa, cov, weights=np.ones(K))
    report("D4  unweighted L2 (w=1)", h4)

    # D5 IRLS (proper Poisson weights from current iter's μ)
    h_curr = h_true.copy()  # warm-start with true uniform
    for it in range(8):
        mu = cov * (kmer_pa.T @ h_curr) + 1e-3
        weights = 1.0 / np.sqrt(mu)
        h_new = solve_custom(counts, kmer_pa, cov, weights=weights)
        delta = np.linalg.norm(h_new - h_curr)
        h_curr = h_new
        if delta < 1e-4: break
    report(f"D5  IRLS (warm-start truth, {it+1} iters)", h_curr)

    # D5b IRLS from uniform start
    h_curr = np.full(F, 1.0/F)
    for it in range(8):
        mu = cov * (kmer_pa.T @ h_curr) + 1e-3
        weights = 1.0 / np.sqrt(mu)
        h_new = solve_custom(counts, kmer_pa, cov, weights=weights)
        delta = np.linalg.norm(h_new - h_curr)
        h_curr = h_new
        if delta < 1e-4: break
    report(f"D5b IRLS (uniform start, {it+1} iters)", h_curr)

    # D6 weights = 1/sqrt(AC)
    weights6 = 1.0 / np.sqrt(np.maximum(ac, 1))
    h6 = solve_custom(counts, kmer_pa, cov, weights=weights6)
    report("D6  weights = 1/sqrt(AC)", h6)

    # D7 weights = 1/(AC + 1)
    weights7 = 1.0 / (ac + 1.0)
    h7 = solve_custom(counts, kmer_pa, cov, weights=weights7)
    report("D7  weights = 1/(AC+1)", h7)

    # D8 weights = sqrt(AC)
    weights8 = np.sqrt(ac.astype(float))
    h8 = solve_custom(counts, kmer_pa, cov, weights=weights8)
    report("D8  weights = sqrt(AC)", h8)

    # ============ D9: relative-frequency formulation ============
    # Drop the ω term and the absolute coverage; fit p_k = c_k / sum(c_k) vs kmer_pa^T h / sum(kmer_pa^T h)
    # i.e., regress relative frequencies, not absolute counts
    print(f"\n--- D9: relative-frequency regression ---")
    p_obs = counts / counts.sum()
    pred_normalizer = kmer_pa.sum(axis=0).sum()  # for unit h, sum(kmer_pa^T h) = sum(kmer_pa[f,:]) = scalar
    # Actually for arbitrary h, sum(kmer_pa^T h) = Σ_k Σ_f kmer_pa[f,k] h[f] = Σ_f h[f] Σ_k kmer_pa[f,k] = h^T per_founder
    # Use a non-normalized formulation:
    # p_pred = (kmer_pa^T h) / (h^T per_founder)
    # This is non-linear. Approximate by replacing denominator with predicted total.
    # Simpler: use absolute counts but with PROPER coverage
    h9 = solve_custom(counts, kmer_pa, cov, weights=np.ones(K))
    report("D9  (= D4, unweighted L2)", h9)

    # ============ D10: solve at TRUE coverage ============
    # If we knew exactly what coverage was, would it work?
    # Assuming uniform h, true coverage λ such that mean count = λ × mean(AC)/F
    # mean count obs = counts.mean()  → λ_true = counts.mean() × F / ac.mean()
    cov_alt = counts.mean() * F / ac.mean()
    print(f"\n--- D10: alternative coverage estimates ---")
    for c_try in [cov, cov_alt, cov * 0.5, cov * 2.0, cov * 5.0]:
        hh = solve_custom(counts, kmer_pa, c_try, weights=np.ones(K))
        report(f"D10 cov={c_try:.1f}× (unweighted L2)", hh)


if __name__ == "__main__":
    main()
