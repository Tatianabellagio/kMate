"""Li–Stephens-style smoothing of per-window founder-frequency (h) estimates.

`smooth_h_across_blocks` is the only production entry point in this module — it is
imported by `per_sample_per_chrom.py` for window (`star2`) mode to smooth per-window
h vectors across adjacent windows.

The experimental BigLD per-block-haplotype EM driver that used to live here
(`load_block_hap_cn`, `project_h_class_to_founders`, `run_chrom_bigld_haplotype`)
was a *tried* approach, not production; it now lives in
`src/archive/block_haplotype_bigld.py`. See `docs/METHODS_TRIED_AND_RESULTS.md`.
"""
from __future__ import annotations
import numpy as np


def smooth_h_across_blocks(
    h_founder_per_block: dict[int, np.ndarray],
    block_pos_start: np.ndarray,
    block_pos_end: np.ndarray,
    n_eco: int,
    recomb_rate: float = 4e-8,
    alpha: float = 0.7,
    n_passes: int = 3,
):
    """Li-Stephens-style HMM smoothing pass on per-block h_founder estimates.

    For each block b, blend its h with a recomb-distance-weighted average of
    its left and right valid neighbors:

        h_b' = alpha * h_b + (1 - alpha) * weighted_avg(h_neighbors)

    Apply n_passes times. Keeps blocks with sparse evidence (single-founder
    dominant) consistent with neighbors that share ancestry, addressing the
    "winner-takes-all" failure mode in heavy-recombination regimes.

    h_founder_per_block: dict {block_idx: (n_eco,) array}, only for blocks
        with successful EM (others are skipped).
    Returns: a new dict with the same keys, smoothed h vectors.
    """
    if not h_founder_per_block or n_passes <= 0 or alpha >= 1.0:
        return h_founder_per_block

    valid_blocks = sorted(h_founder_per_block.keys())
    n_blocks_total = len(block_pos_start)
    mid_pos = (block_pos_start + block_pos_end) / 2.0

    current = dict(h_founder_per_block)
    for _it in range(n_passes):
        new_h: dict[int, np.ndarray] = {}
        for b in valid_blocks:
            h_b = current[b]
            # Find left/right valid neighbors
            left_b = right_b = None
            for j in range(b - 1, -1, -1):
                if j in current:
                    left_b = j
                    break
            for j in range(b + 1, n_blocks_total):
                if j in current:
                    right_b = j
                    break
            neighbor_avg = None
            if left_b is not None and right_b is not None:
                gap_l = abs(mid_pos[b] - mid_pos[left_b])
                gap_r = abs(mid_pos[b] - mid_pos[right_b])
                w_l = float(np.exp(-recomb_rate * gap_l))
                w_r = float(np.exp(-recomb_rate * gap_r))
                norm = w_l + w_r
                if norm > 0:
                    neighbor_avg = (w_l * current[left_b] + w_r * current[right_b]) / norm
            elif left_b is not None:
                neighbor_avg = current[left_b]
            elif right_b is not None:
                neighbor_avg = current[right_b]
            if neighbor_avg is not None:
                h_smoothed = alpha * h_b + (1.0 - alpha) * neighbor_avg
                # Renormalize to simplex
                s = h_smoothed.sum()
                if s > 0:
                    h_smoothed = h_smoothed / s
                new_h[b] = h_smoothed.astype(np.float32)
            else:
                new_h[b] = h_b
        current = new_h
    return current
