"""
Per-LD-block founder-frequency estimation.

Architecture (HAFpipe / hapFIRE style):
    For evolved pool samples, chromosomes are mosaic — different genomic
    regions carry different ancestry mixtures because of recombination +
    selection. A single genome-wide h discards this signal. Per-block h
    captures local ancestry, and per-record alt_freq becomes:

        alt_freq(record) = h_local[block(record)] @ cn_var[:, record]

Block definitions:
    - 'window': fixed bp windows (default 200 kb). Simple; ignores LD.
    - 'ld' (TODO): LD-cluster based on founder-genotype correlation
      between adjacent bubbles. Variable-length, biologically correct.
      hapFIRE uses BigLD.

For a window of W kb in Arabidopsis (selfing ~97%), expected founder
haplotypes per window is small (10–50 distinct), so EM is well-conditioned
even with K_window ~ 100 k-mers/founder.
"""
from __future__ import annotations
import time
from dataclasses import dataclass

import numpy as np
from em_solver import solve_em


@dataclass
class BlockSpec:
    """Genomic block: (chrom, start, end) in 1-based coordinates."""
    chrom: str
    start: int
    end: int

    def __repr__(self):
        return f"{self.chrom}:{self.start}-{self.end}"


def define_windows(bubble_chrom, bubble_start, bubble_end,
                   window_bp=200_000, window_step: int | None = None):
    """Define fixed-bp windows covering all bubbles.

    Args:
        window_bp: width of each window in bp.
        window_step: step between window starts. If None or == window_bp, windows
            are disjoint (legacy behaviour). If < window_bp, windows overlap by
            (window_bp - window_step). E.g. window_bp=10000, window_step=3000 →
            10 kb windows stepping every 3 kb, each k-mer in ~3-4 covering
            windows. Used for Route 2 (overlapping-window smoothing).

    Returns list of BlockSpec tuples, ordered by (chrom, start).
    """
    if window_step is None:
        window_step = window_bp
    if window_step <= 0:
        raise ValueError(f"window_step must be > 0, got {window_step}")
    blocks = []
    for chrom in sorted(set(bubble_chrom)):
        m = bubble_chrom == chrom
        if not m.any():
            continue
        chrom_max = int(bubble_end[m].max())
        n_w = (chrom_max + window_step - 1) // window_step
        for w in range(n_w):
            start = w * window_step + 1
            end = start + window_bp - 1
            blocks.append(BlockSpec(chrom=chrom, start=start, end=end))
    return blocks


def assign_kmers_to_blocks(kmer_bubble_id, bubble_chrom, bubble_start,
                           bubble_end, blocks):
    """For each k-mer, return its block index (or -1 if no block contains it).

    A k-mer is assigned to the block whose [start, end] range overlaps its
    bubble centroid. Bubbles spanning a block boundary are assigned to the
    block containing their start.

    Returns a 1D K-vec assuming disjoint windows (each k-mer in exactly one
    block).
    """
    K = len(kmer_bubble_id)
    bubble_centroid = (bubble_start + bubble_end) // 2

    # Precompute per-chrom block index ranges for fast lookup
    chrom_blocks = {}
    for i, b in enumerate(blocks):
        chrom_blocks.setdefault(b.chrom, []).append((b.start, b.end, i))
    for c in chrom_blocks:
        chrom_blocks[c].sort()

    # For each bubble, find its block
    bubble_block = np.full(len(bubble_chrom), -1, dtype=np.int32)
    for b_idx in range(len(bubble_chrom)):
        c = str(bubble_chrom[b_idx])
        p = int(bubble_centroid[b_idx])
        if c not in chrom_blocks:
            continue
        # Linear scan within chromosome — could binary search, but block
        # count per chrom is small (~150 for 200kb on Chr1)
        for s, e, idx in chrom_blocks[c]:
            if s <= p <= e:
                bubble_block[b_idx] = idx
                break

    # Map k-mers via their bubble_id
    kmer_block = bubble_block[kmer_bubble_id]
    return kmer_block


def assign_records_to_blocks(record_chrom, record_pos, blocks):
    """Per-cn_var-record block index. Same scheme as assign_kmers_to_blocks."""
    chrom_blocks = {}
    for i, b in enumerate(blocks):
        chrom_blocks.setdefault(b.chrom, []).append((b.start, b.end, i))
    for c in chrom_blocks:
        chrom_blocks[c].sort()

    rec_block = np.full(len(record_chrom), -1, dtype=np.int32)
    # Group records by chrom and use searchsorted for vectorization
    for c, items in chrom_blocks.items():
        starts = np.array([s for s, e, idx in items])
        ends = np.array([e for s, e, idx in items])
        idxs = np.array([idx for s, e, idx in items])
        mask = record_chrom == c
        if not mask.any():
            continue
        pos_c = record_pos[mask]
        # find the block whose start <= pos <= end. Since blocks are
        # contiguous and sorted, binary search by start.
        i_starts = np.searchsorted(starts, pos_c, side="right") - 1
        i_starts = np.clip(i_starts, 0, len(starts) - 1)
        in_range = (pos_c >= starts[i_starts]) & (pos_c <= ends[i_starts])
        block_for_rec = np.where(in_range, idxs[i_starts], -1)
        rec_block[mask] = block_for_rec
    return rec_block


def solve_em_per_block(counts, cn_kmer_dense, kmer_block, n_blocks,
                       coverage, em_max_iter=200, tol=1e-7,
                       min_kmers_per_block=200, verbose=False,
                       n_workers=4,
                       global_anchor_weight: float = 0.0,
                       omega=None):
    """Run EM independently per block.

    omega: optional K-vec of per-k-mer weights ω_k (e.g. 1/m_b). Sliced per block
    and passed to solve_em (omega=None → unweighted MLE; identical to old behavior).

    Args:
        counts: K-vec of observed counts (already filtered to block coverage)
        cn_kmer_dense: F × K float32 cn matrix
        kmer_block: K-vec of block index per k-mer (-1 = no block)
        n_blocks: number of blocks
        coverage: scalar coverage estimate
        min_kmers_per_block: blocks with fewer NONZERO-count k-mers fall
            back to a global-h estimate (their h is undefined locally).
        n_workers: thread workers for per-block EM (numpy releases GIL during
            BLAS, so threading shares memory. Default 4. Set to 1 for serial.)
        global_anchor_weight: λ ≥ 0. If > 0, per-block EM is solved with a
            Dirichlet pseudocount centered on `global_h` (the chrom-wide EM
            solution that's already computed for fallback). λ = 0 → pure
            per-window MLE (legacy behaviour). 0.05–0.5 = mild anchor.

    Returns:
        h_blocks: (n_blocks, F) float32 — per-block ancestry vectors. For
            blocks with too few k-mers, falls back to global-h (computed
            once over all k-mers).
        block_status: (n_blocks,) int — 0=local fit, 1=global fallback,
            2=excluded (no k-mers).
        global_h: (F,) — fallback value
    """
    from concurrent.futures import ThreadPoolExecutor
    from threadpoolctl import threadpool_limits

    F = cn_kmer_dense.shape[0]
    h_blocks = np.zeros((n_blocks, F), dtype=np.float32)
    block_status = np.full(n_blocks, 2, dtype=np.int8)

    # Global fallback: EM on ALL k-mers
    if verbose:
        print(f"  computing global-h fallback...", flush=True)
    t = time.time()
    global_h, info = solve_em(counts, cn_kmer_dense, coverage,
                              max_iter=em_max_iter, tol=tol,
                              omega=omega)
    global_h = global_h.astype(np.float32)
    if verbose:
        print(f"    global EM: {info['iterations']} iters, {time.time()-t:.0f}s",
              flush=True)
        if global_anchor_weight > 0:
            print(f"    global_anchor_weight λ={global_anchor_weight} — "
                  f"per-window EM will be anchored toward h_global",
                  flush=True)

    nz = counts > 0
    t0 = time.time()

    # Pre-compute per-block k-mer index lists once.
    block_kmer_idx = [np.flatnonzero(kmer_block == b) for b in range(n_blocks)]

    def _fit_one(b):
        idxs = block_kmer_idx[b]
        if len(idxs) == 0:
            return b, global_h, 2
        idxs_nz = idxs[nz[idxs]]
        if len(idxs_nz) < min_kmers_per_block:
            return b, global_h, 1
        cn_b = np.ascontiguousarray(cn_kmer_dense[:, idxs_nz])
        c_b = counts[idxs_nz]
        omega_b = None if omega is None else omega[idxs_nz]
        if global_anchor_weight > 0:
            h_b, _ = solve_em(c_b, cn_b, coverage,
                              max_iter=em_max_iter, tol=tol,
                              prior_h=global_h,
                              prior_weight=global_anchor_weight,
                              omega=omega_b)
        else:
            h_b, _ = solve_em(c_b, cn_b, coverage,
                              max_iter=em_max_iter, tol=tol,
                              omega=omega_b)
        return b, h_b.astype(np.float32), 0

    # Limit per-thread BLAS to avoid oversubscription. n_workers × inner_threads
    # should ≈ total cores.
    inner_threads = max(1, 8 // n_workers)
    if n_workers == 1:
        results = [_fit_one(b) for b in range(n_blocks)]
    else:
        with threadpool_limits(limits=inner_threads):
            with ThreadPoolExecutor(max_workers=n_workers) as ex:
                results = list(ex.map(_fit_one, range(n_blocks)))

    n_local = n_global = n_excluded = 0
    for b, h_b, st in results:
        h_blocks[b] = h_b
        block_status[b] = st
        if st == 0: n_local += 1
        elif st == 1: n_global += 1
        else: n_excluded += 1

    if verbose:
        print(f"  per-block EM: {n_local} local fits, {n_global} global "
              f"fallbacks, {n_excluded} excluded; {time.time()-t0:.0f}s "
              f"(workers={n_workers}, inner_threads={inner_threads})",
              flush=True)
    return h_blocks, block_status, global_h


def project_blocks_to_records(h_blocks, global_h, cn_var, record_block,
                              cn_var_called=None):
    """Project per-block h to per-record alt freqs by hard window assignment.

    Each record uses the h of the window containing it, h_blocks[record_block[r]]
    (already set to global_h for fallback/excluded windows by solve_em_per_block;
    records with no window, record_block == -1, use global_h). With HMM smoothing
    applied upstream, this hard assignment is the production projection.

    With cn_var_called, AF is the missing-aware (h@cn_var)/(h@cn_var_called);
    without it, ./. is treated as REF (→ SV AF under-call at high missingness).

    Returns (alt_freq, info), both length-N float64. info[r] is the projection
    denominator — the h-weighted called mass at r (1.0 when no called mask).
    """
    N = cn_var.shape[1]
    out = np.zeros(N, dtype=np.float64)
    info = np.ones(N, dtype=np.float64)
    rb = np.asarray(record_block)
    cn_var_csc = cn_var.tocsc()
    cvc_csc = cn_var_called.tocsc() if cn_var_called is not None else None
    for b in np.unique(rb):
        mask = (rb == b)
        if not mask.any():
            continue
        h = (h_blocks[b] if b >= 0 else global_h).astype(np.float32)
        num = (h @ cn_var_csc[:, mask].toarray()).astype(np.float64)
        if cvc_csc is not None:
            den = (h @ cvc_csc[:, mask].toarray()).astype(np.float64)
            out[mask] = num / np.maximum(den, 1e-12)
            info[mask] = den
        else:
            out[mask] = num
    return out, info
