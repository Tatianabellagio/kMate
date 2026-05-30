"""
Try unconstrained NNLS: solve for q ≥ 0 (founder rates), post-normalize to get h.

Standard NNLS is well-behaved and uses scipy.optimize.nnls. No simplex pathology.
"""
from __future__ import annotations
import sys, os, time
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import numpy as np
from scipy.sparse import load_npz
from scipy.optimize import nnls
import cvxpy as cp


DATA = os.path.join(os.path.dirname(__file__), "..", "data")


def main():
    cn_sparse = load_npz(os.path.join(DATA, "test_chr1_first200.kmer_pa.npz"))
    meta = np.load(os.path.join(DATA, "test_chr1_first200.meta.npz"), allow_pickle=True)
    founders = meta["founders"]
    F, K = cn_sparse.shape
    kmer_pa = np.asarray(cn_sparse.todense()).astype(np.int8)
    ac = kmer_pa.sum(axis=0)
    counts = np.load(os.path.join(DATA, "sim_chr1", "uniform82_counts.npz"))["counts"]
    cov = counts.sum() * F / ac.sum()
    h_true = np.full(F, 1.0/F)

    print(f"kmer_pa: {F} × {K:,}, cov={cov:.1f}×, truth=1/{F}={1/F:.4f}")

    def report(label, h_hat, raw_q=None):
        h_hat = np.asarray(h_hat)
        h_hat = h_hat / h_hat.sum() if h_hat.sum() > 0 else h_hat
        err = np.linalg.norm(h_hat - h_true)
        cv = h_hat.std() / max(1e-9, h_hat.mean())
        nz = (h_hat < 1e-4).sum()
        print(f"  {label:55s}  ||Δh||={err:.4f}  CV={cv:.3f}  n_zero={nz}/{F}")
        if raw_q is not None:
            print(f"    raw q stats: min={raw_q.min():.4f}  max={raw_q.max():.4f}  "
                  f"mean={raw_q.mean():.4f}  median={np.median(raw_q):.4f}")
        return err

    print(f"\n{'-'*72}\nUnconstrained methods\n{'-'*72}")

    # U1: scipy NNLS (no upper limit, just q ≥ 0). Solves Σ (c - kmer_pa^T q)^2.
    # Note: solving for q (founder absolute rate, includes coverage)
    # Predicted: kmer_pa^T q. Observed: c.
    print(f"\n  Solving NNLS (full system, may be slow)...")
    t = time.time()
    q_nnls, residual = nnls(kmer_pa.T.astype(float), counts.astype(float), maxiter=5000)
    print(f"  [{time.time()-t:.0f}s, residual={residual:.0f}]")
    h_u1 = q_nnls / q_nnls.sum() if q_nnls.sum() > 0 else q_nnls
    report("U1  scipy NNLS", h_u1, raw_q=q_nnls)

    # U2: CVXPY without simplex constraint (just q ≥ 0)
    print(f"\n  Solving CVXPY q ≥ 0 only (no sum constraint)...")
    t = time.time()
    q = cp.Variable(F, nonneg=True)
    pred = kmer_pa.T @ q
    obj = cp.Minimize(cp.sum_squares(counts - pred))
    cp.Problem(obj).solve(solver="SCS", verbose=False)
    q_u2 = q.value
    print(f"  [{time.time()-t:.0f}s]")
    h_u2 = q_u2 / q_u2.sum() if q_u2.sum() > 0 else q_u2
    report("U2  CVXPY q ≥ 0 (no simplex)", h_u2, raw_q=q_u2)

    # U3: Poisson NLL without simplex
    print(f"\n  Solving Poisson NLL, q ≥ 0...")
    t = time.time()
    q = cp.Variable(F, nonneg=True)
    mu = kmer_pa.T @ q + 1e-3
    obj = cp.Minimize(cp.sum(mu) - cp.sum(cp.multiply(counts, cp.log(mu))))
    cp.Problem(obj).solve(solver="SCS", verbose=False)
    q_u3 = q.value
    print(f"  [{time.time()-t:.0f}s]")
    h_u3 = q_u3 / q_u3.sum() if q_u3.sum() > 0 else q_u3
    report("U3  Poisson NLL, q ≥ 0", h_u3, raw_q=q_u3)

    # U4: KL divergence, q ≥ 0
    print(f"\n  Solving KL divergence, q ≥ 0...")
    t = time.time()
    q = cp.Variable(F, nonneg=True)
    mu = kmer_pa.T @ q + 1e-3
    obj = cp.Minimize(cp.sum(cp.kl_div(counts, mu)))
    cp.Problem(obj).solve(solver="SCS", verbose=False)
    q_u4 = q.value
    print(f"  [{time.time()-t:.0f}s]")
    h_u4 = q_u4 / q_u4.sum() if q_u4.sum() > 0 else q_u4
    report("U4  KL div, q ≥ 0", h_u4, raw_q=q_u4)

    # U5: Identifiability check. The "true" q under uniform truth is (cov/F) for each founder.
    # Check that the system at h=truth fits well.
    q_true = np.full(F, cov/F)
    pred_true = kmer_pa.T @ q_true
    rmse_true = np.sqrt(np.mean((counts - pred_true)**2))
    print(f"\n  At TRUE q (uniform cov/F={cov/F:.3f}):")
    print(f"    RMSE: {rmse_true:.2f}")
    print(f"    expected Poisson noise: sqrt(mean count) = {np.sqrt(counts.mean()):.2f}")
    print(f"    sum(predicted) = {pred_true.sum():.0f}  vs sum(obs) = {counts.sum()}")

    # U6: Show top vs bottom for U1 (scipy NNLS, ground truth with non-neg constraint)
    print(f"\n  U1 (scipy NNLS) top vs bottom founders:")
    for i in np.argsort(-h_u1)[:5]:
        print(f"    top {founders[i]}: q={q_nnls[i]:.4f}  h={h_u1[i]:.4f}")
    for i in np.argsort(h_u1)[:5]:
        print(f"    bot {founders[i]}: q={q_nnls[i]:.4f}  h={h_u1[i]:.4f}")


if __name__ == "__main__":
    main()
