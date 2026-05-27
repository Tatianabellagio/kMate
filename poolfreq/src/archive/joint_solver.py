"""
Joint hapFIRE-style solver: per-bubble haplotype reduction + ONE global founder vector.

For each bubble v:
  - Reduce cn[:, kmers_v] (F × K_v) to unique haplotype patterns (H_v × K_v)
  - Estimate haplotype freqs p_v via EM (well-conditioned: typically H_v ~5-6, K_v ~50-300)
  - Build founder→haplotype dosage matrix M_v (H_v × F): M_v[h, f] = 1 iff founder f's
    haplotype at bubble v is h

JOINT projection:
  Variables:  h ∈ R^F (single global founder freq vector)
  Objective:  minimize Σ_v λ_v ‖M_v · h - p_v‖²
  Constraints: h ≥ 0, sum h = 1

This pools evidence from ALL bubbles to constrain ONE founder vector. Vastly
over-determined (each bubble adds ~5-6 constraints, hundreds of bubbles).

This is the actual hapFIRE-style architecture — what their CVXPY block does,
adapted for our pangenome data structure.
"""
from __future__ import annotations
import numpy as np
import cvxpy as cp
from hapfire_solver import reduce_to_unique_haplotypes, solve_em_simple


def stage1_block(counts_b, cn_b, max_iter=100, tol=1e-7):
    """Stage 1: reduce one block's cn matrix and run EM on haplotype freqs.

    Returns:
        p_haplotype:  H_v-vector of haplotype frequencies (sum=1)
        M:            H_v × F matrix; M[h, f] = 1 iff founder f's haplotype is h
    """
    F = cn_b.shape[0]
    cn_uniq, founder_to_hap, hap_to_founders = reduce_to_unique_haplotypes(cn_b)
    H = cn_uniq.shape[0]
    p, _ = solve_em_simple(counts_b, cn_uniq, max_iter=max_iter, tol=tol)
    M = np.zeros((H, F), dtype=np.float32)
    for h_idx, founders in enumerate(hap_to_founders):
        for f in founders:
            M[h_idx, f] = 1
    return p, M


def solve_joint(counts, cn, bubble_id, weights_per_block=None, verbose=False):
    """Joint solve: stage 1 per bubble, then JOINT CVXPY for global h.

    Parameters:
        counts:       K-vector of k-mer counts (full)
        cn:           F × K binary copy-number matrix
        bubble_id:    K-vector mapping each k-mer to its bubble (0..n_bubbles-1)
        weights_per_block: optional weights to scale each block's residual contribution

    Returns:
        h:           F-vector of founder frequencies
        diagnostics: dict with per-block info
    """
    F = cn.shape[0]
    n_bubbles = len(set(bubble_id))

    # Stage 1: per-bubble haplotype reduction + EM
    p_per_block = []   # H_v vectors
    M_per_block = []   # H_v × F matrices
    bubble_indices = sorted(set(bubble_id))
    for b in bubble_indices:
        mask = bubble_id == b
        if mask.sum() < 3:
            continue
        cn_b = cn[:, mask]
        c_b = counts[mask]
        try:
            p, M = stage1_block(c_b, cn_b)
            if p.size < 2:
                continue  # no haplotype variation, skip
            p_per_block.append(p)
            M_per_block.append(M)
        except Exception as e:
            if verbose:
                print(f"  block {b}: stage1 failed ({e})")
            continue

    n_blocks_solved = len(p_per_block)
    if verbose:
        print(f"  Stage 1: solved {n_blocks_solved} blocks")
        print(f"  H_v distribution: min={min(p.size for p in p_per_block)}, "
              f"max={max(p.size for p in p_per_block)}, "
              f"median={int(np.median([p.size for p in p_per_block]))}")

    # Stage 2: JOINT CVXPY
    h = cp.Variable(F, nonneg=True)
    if weights_per_block is None:
        weights_per_block = [1.0] * n_blocks_solved

    losses = []
    for w, p, M in zip(weights_per_block, p_per_block, M_per_block):
        # M is H_v × F, so M @ h is H_v-vector predicting per-haplotype freq
        diff = M @ h - p
        losses.append(w * cp.sum_squares(diff))

    obj = cp.Minimize(sum(losses))
    cons = [cp.sum(h) == 1]
    prob = cp.Problem(obj, cons)
    prob.solve(solver="SCS", verbose=verbose)

    if h.value is None:
        if verbose:
            print(f"  CVXPY solver failed: status={prob.status}")
        return None, {"failure": True}

    return np.asarray(h.value), {
        "n_blocks_solved": n_blocks_solved,
        "objective": prob.value,
        "H_v_distribution": {
            "min": min(p.size for p in p_per_block),
            "max": max(p.size for p in p_per_block),
            "median": int(np.median([p.size for p in p_per_block])),
        },
    }


def solve_joint_l1(counts, cn, bubble_id, verbose=False):
    """Same as solve_joint but with L1 loss (more robust to outlier blocks)."""
    F = cn.shape[0]
    p_per_block, M_per_block = [], []
    for b in sorted(set(bubble_id)):
        mask = bubble_id == b
        if mask.sum() < 3:
            continue
        try:
            p, M = stage1_block(counts[mask], cn[:, mask])
            if p.size < 2:
                continue
            p_per_block.append(p)
            M_per_block.append(M)
        except Exception:
            continue
    h = cp.Variable(F, nonneg=True)
    losses = [cp.norm1(M @ h - p) for p, M in zip(p_per_block, M_per_block)]
    prob = cp.Problem(cp.Minimize(sum(losses)), [cp.sum(h) == 1])
    prob.solve(solver="SCS", verbose=verbose)
    return np.asarray(h.value), {"n_blocks": len(p_per_block)}
