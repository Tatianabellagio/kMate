"""
Sanity-check the block solver on synthetic data with known truth.

Tests:
  T1 — clean k-mers, no contamination: solver recovers h to within noise
  T2 — over-determined system: SNP-rich block (K >> F) recovers exact h
  T3 — under-determined: K < F, solver still finds a feasible h on simplex
  T4 — with off-target contamination ω: solver recovers h when ω is given
  T5 — ω alternate-fit: starting from ω=0, recover both h and ω
  T6 — coverage scaling: changing λ shouldn't change the recovered h
  T7 — IRLS vs Anscombe agreement on clean data
"""
from __future__ import annotations
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import numpy as np
from block_solver import solve_block_wls, solve_block_irls, update_omega


def gen_synth(F=20, K=200, n_samples=1, coverage=30, seed=0,
              contamination_frac=0.0, contamination_rate=0.0):
    """Simulate one block.

    Returns:
        h_true:   n_samples × F (true founder freqs)
        counts:   n_samples × K (Poisson-sampled counts)
        cn:       F × K (random binary copy numbers, ~50% allele freq)
        omega:    K-vector of true contamination rates
    """
    rng = np.random.default_rng(seed)
    # Random binary cn matrix; mean carrier freq ~0.3 to 0.5
    af = rng.uniform(0.1, 0.5, K)
    cn = (rng.random((F, K)) < af[None, :]).astype(np.int8)

    # Random h on simplex per sample
    h_true = rng.dirichlet(np.ones(F) * 0.5, size=n_samples)

    # Optional contamination
    omega = np.zeros(K)
    if contamination_frac > 0:
        idx = rng.choice(K, int(K * contamination_frac), replace=False)
        omega[idx] = contamination_rate

    # Expected counts
    mu = coverage * (h_true @ cn) + coverage * omega[None, :]
    mu = np.maximum(mu, 1e-6)
    counts = rng.poisson(mu)
    return h_true, counts, cn, omega


def test_clean_overdetermined():
    """T1+T2: clean k-mers, K >> F. Should recover h tightly."""
    h_true, counts, cn, omega = gen_synth(F=20, K=500, coverage=50, seed=1)
    h_hat, _ = solve_block_wls(counts[0], cn, coverage=50.0)
    err = np.linalg.norm(h_hat - h_true[0])
    print(f"  T1/T2: F=20, K=500, cov=50 → ||h_hat - h_true|| = {err:.4f}")
    print(f"         max |Δh|: {np.abs(h_hat - h_true[0]).max():.4f}")
    print(f"         h_true[0:5] = {h_true[0,:5].round(3)}")
    print(f"         h_hat [0:5] = {h_hat[:5].round(3)}")
    assert err < 0.05, f"recovery error {err:.4f} too high"


def test_under_determined():
    """T3: K < F, system under-determined. Solver should still produce valid simplex output."""
    h_true, counts, cn, _ = gen_synth(F=80, K=40, coverage=30, seed=2)
    h_hat, _ = solve_block_wls(counts[0], cn, coverage=30.0)
    print(f"\n  T3: F=80, K=40 → simplex check:")
    print(f"         min(h)={h_hat.min():.4f}, sum(h)={h_hat.sum():.4f}")
    assert h_hat.min() >= -1e-6, "h has negative entries"
    assert abs(h_hat.sum() - 1) < 1e-3, f"h doesn't sum to 1 ({h_hat.sum():.4f})"


def test_known_omega():
    """T4: 30% of k-mers have contamination at rate 5 reads/coverage. With omega given, recover h."""
    h_true, counts, cn, omega = gen_synth(F=20, K=500, coverage=50,
                                           contamination_frac=0.3, contamination_rate=5.0, seed=3)
    h_hat, _ = solve_block_wls(counts[0], cn, coverage=50.0, omega=omega)
    err = np.linalg.norm(h_hat - h_true[0])
    print(f"\n  T4: with known ω → ||h_hat - h_true|| = {err:.4f}")
    assert err < 0.08, f"recovery error {err:.4f} too high with known omega"


def test_alternate_omega():
    """T5: alternate-solve recovery. Start from ω=0, fit h, re-estimate ω, refit h, ..."""
    F, K, N = 20, 500, 50  # 50 samples, identifies omega well
    h_true, counts, cn, omega_true = gen_synth(F=F, K=K, n_samples=N, coverage=50,
                                                contamination_frac=0.3, contamination_rate=5.0, seed=4)
    coverage = np.full(N, 50.0)

    # Initial: omega = 0
    omega_hat = np.zeros(K)
    h_hat = np.zeros((N, F))
    for outer in range(3):
        for i in range(N):
            h_hat[i], _ = solve_block_wls(counts[i], cn, coverage=coverage[i], omega=omega_hat)
        omega_hat = update_omega(counts, cn, h_hat, coverage)
        err_h = np.linalg.norm(h_hat - h_true) / N
        err_omega = np.linalg.norm(omega_hat - omega_true)
        print(f"\n  T5 outer iter {outer}: avg ||Δh||={err_h:.4f}, ||Δω||={err_omega:.4f}")
    # After 3 iterations should be close
    assert err_h < 0.1, f"avg h error {err_h:.4f} too high"
    assert err_omega < 5, f"omega error {err_omega:.4f} too high"


def test_coverage_invariance():
    """T6: doubling coverage shouldn't change recovered h."""
    h_true, counts1, cn, _ = gen_synth(F=20, K=500, coverage=20, seed=5)
    h_true2, counts2, _, _ = gen_synth(F=20, K=500, coverage=80, seed=5)  # same seed → same h_true
    assert np.allclose(h_true, h_true2)
    h1, _ = solve_block_wls(counts1[0], cn, coverage=20.0)
    h2, _ = solve_block_wls(counts2[0], cn, coverage=80.0)
    print(f"\n  T6: ||h@cov20 - h@cov80|| = {np.linalg.norm(h1 - h2):.4f}")
    print(f"      ||h@cov20 - h_true||   = {np.linalg.norm(h1 - h_true[0]):.4f}")
    print(f"      ||h@cov80 - h_true||   = {np.linalg.norm(h2 - h_true[0]):.4f}")


def test_irls_vs_anscombe():
    """T7: agreement on clean data."""
    h_true, counts, cn, _ = gen_synth(F=20, K=500, coverage=50, seed=7)
    h_a, _ = solve_block_wls(counts[0], cn, coverage=50.0)
    h_i, _ = solve_block_irls(counts[0], cn, coverage=50.0)
    print(f"\n  T7: ||h_anscombe - h_irls|| = {np.linalg.norm(h_a - h_i):.4f}")
    print(f"      ||h_anscombe - h_true|| = {np.linalg.norm(h_a - h_true[0]):.4f}")
    print(f"      ||h_irls - h_true||     = {np.linalg.norm(h_i - h_true[0]):.4f}")


if __name__ == "__main__":
    print("="*70)
    print("Block solver synthetic tests")
    print("="*70)
    test_clean_overdetermined()
    test_under_determined()
    test_known_omega()
    test_alternate_omega()
    test_coverage_invariance()
    test_irls_vs_anscombe()
    print("\nAll assertions passed.")
