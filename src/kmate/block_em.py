"""
Per-LD-block founder-frequency estimation.

Architecture (HAFpipe / hapFIRE style):
    For evolved pool samples, chromosomes are mosaic — different genomic
    regions carry different ancestry mixtures because of recombination +
    selection. A single genome-wide h discards this signal. Per-block h
    captures local ancestry, and per-record alt_freq becomes:

        alt_freq(record) = h_local[block(record)] @ var_pa[:, record]

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
from .em_solver import solve_em, haploblock_collapse_indices


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
    """Per-var_pa-record block index. Same scheme as assign_kmers_to_blocks."""
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


def solve_em_per_block(counts, kmer_pa_dense, kmer_block, n_blocks,
                       coverage, em_max_iter=200, tol=1e-7,
                       min_kmers_per_block=200, verbose=False,
                       n_workers=4,
                       global_anchor_weight: float = 0.0,
                       omega=None,
                       local_only: bool = False,
                       normalize: str = "per_founder",
                       haploblock_collapse: bool = True,
                       haploblock_eps: float = 0.0):
    """Run EM independently per block.

    omega: optional K-vec of per-k-mer weights ω_k (e.g. 1/m_b). Sliced per block
    and passed to solve_em (omega=None → unweighted MLE; identical to old behavior).

    haploblock_collapse (default True): before fitting each block, collapse the F
    founders into the DISTINCT HAPLOTYPES present in that block (founders with a
    byte-identical k-mer presence pattern over the block's k-mers — exact identity,
    eps=0). The EM is then solved over the K_b ≤ F distinct haplotypes, not the full
    founder set, and each class frequency is split equally back to its member
    founders. This is the general kMate design (block → haploblock → EM, hapFIRE-like):
    it removes the flat-likelihood ridge from k-mer-indistinguishable founders and
    fits only what the data can resolve. When every founder is distinct in a block
    (K_b == F, e.g. coarse r2=0.1 windows) it is an EXACT no-op vs fitting F directly
    (the Kb==F branch fits founders in native order, no permutation), so
    coarse-block/global-style runs are byte-unchanged. Set False for the legacy
    always-fit-F behaviour.

    haploblock_eps (default 0.0): haplotype-merge tolerance passed to
    haploblock_collapse_indices — 0.0 = exact byte-identical presence patterns
    (lossless); >0 = merge founders whose presence differs in ≤ eps·(block k-mers),
    approximate, see haploblock_collapse_indices.

    Args:
        counts: K-vec of observed counts (already filtered to block coverage)
        kmer_pa_dense: F × K float32 kmer_pa matrix
        kmer_block: K-vec of block index per k-mer (-1 = no block)
        n_blocks: number of blocks
        coverage: scalar coverage estimate
        min_kmers_per_block: blocks with fewer NONZERO-count k-mers fall
            back to a global-h estimate (their h is undefined locally).
        local_only: if True, the model is GLOBAL-FREE. Low-evidence
            (status 1) and empty (status 2) blocks are filled with NaN
            instead of global_h, and the per-block EM is never anchored
            (caller must pass global_anchor_weight=0). global_h is still
            computed and returned for diagnostics/storage, but it never
            enters the per-block h or (via the caller) the AF projection.
            Records in NaN blocks project to NaN AF. Use for recombinant
            pools where the chrom-wide mixture is the wrong prior.
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
    try:
        from threadpoolctl import threadpool_limits
    except ImportError:
        # threadpoolctl only bounds BLAS oversubscription across the parallel
        # per-block fits; the EM is still correct without it (just potentially
        # thread-greedy). Degrade gracefully rather than hard-crash a cohort run
        # if the env is missing the package.
        from contextlib import contextmanager
        @contextmanager
        def threadpool_limits(limits=None):
            yield

    F = kmer_pa_dense.shape[0]
    h_blocks = np.zeros((n_blocks, F), dtype=np.float32)
    block_status = np.full(n_blocks, 2, dtype=np.int8)

    # Global fallback: EM on ALL k-mers
    if verbose:
        print(f"  computing global-h fallback...", flush=True)
    t = time.time()
    global_h, info = solve_em(counts, kmer_pa_dense, coverage,
                              max_iter=em_max_iter, tol=tol,
                              omega=omega, normalize=normalize)
    global_h = global_h.astype(np.float32)
    # Global-free mode: low-evidence/empty blocks get NaN, not global_h.
    fallback_h = np.full(F, np.nan, dtype=np.float32) if local_only else global_h
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
            return b, fallback_h, 2
        idxs_nz = idxs[nz[idxs]]
        if len(idxs_nz) < min_kmers_per_block:
            return b, fallback_h, 1
        c_b = counts[idxs_nz]
        omega_b = None if omega is None else omega[idxs_nz]
        # per_founder normalizer over the FULL window (all its k-mers, incl. c_k=0),
        # not just the observed idxs_nz — else it's conditioned on this run's zero draws.
        w_full = (np.ones(len(idxs), np.float32) if omega is None
                  else omega[idxs].astype(np.float32))
        kfw_full = (kmer_pa_dense[:, idxs] @ w_full).astype(np.float32)   # F-vec

        if haploblock_collapse:
            # block -> distinct haplotypes present -> EM over K_b (not F). Members of a
            # class match over idxs (⊇ idxs_nz), so any representative row stands in.
            lab, reps, csize, Kb = haploblock_collapse_indices(kmer_pa_dense[:, idxs],
                                                               eps=haploblock_eps)
            if Kb < F:                                        # genuine collapse
                cn_c = np.ascontiguousarray(kmer_pa_dense[reps][:, idxs_nz])   # K_b × n_nz
                kfw_c = kfw_full[reps]
                if global_anchor_weight > 0:
                    # NOTE: anchor is summed per class then split EQUALLY back to members
                    # (h_b below), so any INTRA-class prior asymmetry in global_h is lost.
                    # Acceptable because members are (≤eps) k-mer-indistinguishable here so
                    # the data can't separate them anyway; but this is why the production
                    # window recipe is LOCAL-ONLY (global_anchor_weight=0) — see
                    # per_sample_per_chrom.run_one_chrom_window.
                    prior_c = np.bincount(lab, weights=global_h, minlength=Kb).astype(np.float32)
                    h_c, _ = solve_em(c_b, cn_c, coverage, max_iter=em_max_iter, tol=tol,
                                      prior_h=prior_c, prior_weight=global_anchor_weight,
                                      omega=omega_b, normalize=normalize, kfw=kfw_c)
                else:
                    h_c, _ = solve_em(c_b, cn_c, coverage, max_iter=em_max_iter, tol=tol,
                                      omega=omega_b, normalize=normalize, kfw=kfw_c)
                h_b = (h_c / csize)[lab]                      # split class freq equally to members
                return b, h_b.astype(np.float32), 0
            # Kb == F: no distinct-haplotype merging — fall through to the direct
            # F-founder fit (EXACT no-op, avoids a needless row permutation + reorder).

        # direct fit over the full founder set (legacy path, and the Kb==F no-op above)
        cn_b = np.ascontiguousarray(kmer_pa_dense[:, idxs_nz])
        if global_anchor_weight > 0:
            h_b, _ = solve_em(c_b, cn_b, coverage, max_iter=em_max_iter, tol=tol,
                              prior_h=global_h, prior_weight=global_anchor_weight,
                              omega=omega_b, normalize=normalize, kfw=kfw_full)
        else:
            h_b, _ = solve_em(c_b, cn_b, coverage, max_iter=em_max_iter, tol=tol,
                              omega=omega_b, normalize=normalize, kfw=kfw_full)
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


def project_blocks_to_records(h_blocks, global_h, var_pa, record_block,
                              var_called=None):
    """Project per-block h to per-record alt freqs by hard window assignment.

    Each record uses the h of the window containing it, h_blocks[record_block[r]]
    (already set to global_h for fallback/excluded windows by solve_em_per_block;
    records with no window, record_block == -1, use global_h). With HMM smoothing
    applied upstream, this hard assignment is the production projection.

    With var_called, AF is the missing-aware (h@var_pa)/(h@var_called);
    without it, ./. is treated as REF (→ SV AF under-call at high missingness).

    Returns (alt_freq, info), both length-N float64. info[r] is the projection
    denominator — the h-weighted called mass at r (1.0 when no called mask).
    """
    N = var_pa.shape[1]
    out = np.zeros(N, dtype=np.float64)
    info = np.ones(N, dtype=np.float64)
    rb = np.asarray(record_block)
    var_pa_csc = var_pa.tocsc()
    cvc_csc = var_called.tocsc() if var_called is not None else None
    for b in np.unique(rb):
        mask = (rb == b)
        if not mask.any():
            continue
        h = (h_blocks[b] if b >= 0 else global_h).astype(np.float32)
        num = (h @ var_pa_csc[:, mask].toarray()).astype(np.float64)
        if cvc_csc is not None:
            den = (h @ cvc_csc[:, mask].toarray()).astype(np.float64)
            out[mask] = num / np.maximum(den, 1e-12)
            info[mask] = den
        else:
            out[mask] = num
    return out, info
