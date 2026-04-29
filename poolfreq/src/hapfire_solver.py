"""
Two-stage hapFIRE-style solver:
  Stage 1: identify unique haplotype patterns within a block (reduce F × K → H × K)
  Stage 2: EM on the reduced (H × K) matrix → haplotype frequencies (length H)
  Stage 3: project haplotype frequencies → founder frequencies via simplex CVXPY
            (or simple equal-share for a single block)

Within an LD block, two founders that share identical k-mer patterns are
unidentifiable from this block's data. hapFIRE's insight is to reduce the state
space to just the UNIQUE haplotypes (H ≤ F, often << F), making the per-block
EM problem well-conditioned. Cross-block LD then disambiguates which founders
carry which haplotypes (since they typically diverge across blocks).
"""
from __future__ import annotations
import numpy as np
import cvxpy as cp


def reduce_to_unique_haplotypes(cn: np.ndarray) -> tuple[np.ndarray, np.ndarray, list[list[int]]]:
    """Reduce a F × K cn matrix to its unique row patterns.

    Returns:
        cn_unique:        H × K matrix where each row is a unique haplotype pattern
        founder_to_hap:   F-vector mapping each founder to its haplotype index (0..H-1)
        hap_to_founders:  list of lists; hap_to_founders[h] = list of founder indices
                          carrying haplotype h
    """
    F, K = cn.shape
    # Hash each row by tobytes() for unique identification
    seen: dict[bytes, int] = {}
    founder_to_hap = np.zeros(F, dtype=np.int32)
    cn_unique_rows = []
    hap_to_founders: list[list[int]] = []

    for f in range(F):
        key = cn[f].tobytes()
        if key in seen:
            h_idx = seen[key]
            hap_to_founders[h_idx].append(f)
        else:
            h_idx = len(seen)
            seen[key] = h_idx
            cn_unique_rows.append(cn[f].copy())
            hap_to_founders.append([f])
        founder_to_hap[f] = h_idx

    cn_unique = np.array(cn_unique_rows, dtype=cn.dtype)
    return cn_unique, founder_to_hap, hap_to_founders


def solve_em_simple(counts: np.ndarray, cn: np.ndarray, max_iter: int = 100,
                    tol: float = 1e-7) -> tuple[np.ndarray, dict]:
    """Plain EM on counts vs cn matrix (no coverage scaling needed in the M-step)."""
    H, K = cn.shape
    cn = cn.astype(np.float32)
    counts = counts.astype(np.float64)

    p = np.full(H, 1.0 / H)
    for it in range(max_iter):
        denom = p @ cn  # K-vector
        denom = np.maximum(denom, 1e-15)
        cw = counts / denom
        p_new = p * (cn @ cw)
        total = p_new.sum()
        if total > 0:
            p_new = p_new / total
        else:
            p_new = np.full(H, 1.0 / H)
        delta = np.linalg.norm(p_new - p)
        p = p_new
        if delta < tol:
            break
    return p, {"iterations": it + 1, "delta": delta}


def project_haplotypes_to_founders(
    p_haplotype: np.ndarray,        # H-vector of haplotype frequencies (sum=1)
    hap_to_founders: list[list[int]], # for each h, list of founder indices carrying it
    F: int,
    method: str = "equal_share",
) -> np.ndarray:
    """Project haplotype frequencies to founder frequencies.

    Methods:
        equal_share: distribute each haplotype's freq uniformly among its founders
        cvxpy:       hapFIRE-style constrained least-squares projection
    """
    if method == "equal_share":
        h_founder = np.zeros(F)
        for h_idx, founders in enumerate(hap_to_founders):
            if not founders:
                continue
            share = p_haplotype[h_idx] / len(founders)
            for f in founders:
                h_founder[f] += share
        return h_founder

    elif method == "cvxpy":
        H = len(hap_to_founders)
        # Build M: H × F where M[h, f] = 1 if founder f carries haplotype h
        M = np.zeros((H, F), dtype=np.float32)
        for h_idx, founders in enumerate(hap_to_founders):
            for f in founders:
                M[h_idx, f] = 1
        h = cp.Variable(F, nonneg=True)
        diff = M @ h - p_haplotype
        prob = cp.Problem(cp.Minimize(cp.norm(diff)), [cp.sum(h) == 1])
        prob.solve(solver="SCS", verbose=False)
        if h.value is None:
            return project_haplotypes_to_founders(p_haplotype, hap_to_founders, F, "equal_share")
        return np.asarray(h.value)

    else:
        raise ValueError(method)


def solve_block_two_stage(
    counts: np.ndarray,
    cn: np.ndarray,
    coverage: float = None,  # not used internally; kept for API compatibility
    project_method: str = "equal_share",
    em_max_iter: int = 100,
    em_tol: float = 1e-7,
    return_diagnostics: bool = False,
):
    """Full two-stage solve: reduce → EM → project.

    Returns:
        h_founder:  F-vector of estimated founder frequencies (sum=1)
        diagnostics dict (if return_diagnostics)
    """
    F, K = cn.shape
    cn_unique, founder_to_hap, hap_to_founders = reduce_to_unique_haplotypes(cn)
    H = cn_unique.shape[0]

    p_hap, info = solve_em_simple(counts, cn_unique,
                                   max_iter=em_max_iter, tol=em_tol)

    h_founder = project_haplotypes_to_founders(
        p_hap, hap_to_founders, F, method=project_method
    )

    if return_diagnostics:
        return h_founder, {
            "n_unique_haplotypes": H,
            "n_founders": F,
            "haplotype_freqs": p_hap,
            "hap_to_founders": hap_to_founders,
            "em_iterations": info["iterations"],
        }
    return h_founder


def solve_multi_block(
    counts_per_block: list[np.ndarray],
    cn_per_block: list[np.ndarray],
    coverage: float = None,
    project_method: str = "average",
    em_max_iter: int = 100,
):
    """Run two-stage solve on multiple blocks and aggregate to a single h_founder.

    Methods to aggregate:
        average:   per-block solve then average h_founder
        joint_cvxpy:  build joint CVXPY across blocks (TODO)
    """
    F = cn_per_block[0].shape[0]
    h_per_block = []
    for counts, cn in zip(counts_per_block, cn_per_block):
        h_b = solve_block_two_stage(counts, cn, coverage, project_method="equal_share",
                                    em_max_iter=em_max_iter)
        h_per_block.append(h_b)
    return np.mean(h_per_block, axis=0)
