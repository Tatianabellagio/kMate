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

    NOTE: this returns a 1D K-vec assuming disjoint windows (each k-mer in
    exactly one block). For OVERLAPPING windows use
    `assign_kmers_to_blocks_multi` instead, which returns a per-block list.
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


def assign_kmers_to_blocks_multi(kmer_bubble_id, bubble_chrom, bubble_start,
                                  bubble_end, blocks):
    """For each k-mer, return ALL block indices containing its bubble centroid.

    Used for overlapping windows where each k-mer participates in 2+ EM solves.

    Returns:
        block_kmer_idx: list of length n_blocks; entry b is an int64 array of
            k-mer indices that belong to block b.
    """
    bubble_centroid = (bubble_start + bubble_end) // 2

    # Per-chrom: sorted (start, end, block_idx) triples + arrays for searchsorted
    chrom_blocks: dict[str, list[tuple[int, int, int]]] = {}
    for i, b in enumerate(blocks):
        chrom_blocks.setdefault(b.chrom, []).append((b.start, b.end, i))
    for c in chrom_blocks:
        chrom_blocks[c].sort()

    n_blocks = len(blocks)
    # Per-bubble: list of block indices
    bubble_blocks: list[list[int]] = [[] for _ in range(len(bubble_chrom))]
    for b_idx in range(len(bubble_chrom)):
        c = str(bubble_chrom[b_idx])
        p = int(bubble_centroid[b_idx])
        if c not in chrom_blocks:
            continue
        # Walk blocks (sorted by start) — they're contiguous in start order.
        # Worst-case O(blocks/chrom) per bubble; still cheap for ~10k windows.
        for s, e, idx in chrom_blocks[c]:
            if s > p:
                break  # sorted by start; nothing later contains p
            if p <= e:
                bubble_blocks[b_idx].append(idx)

    # Build block_kmer_idx by inverting bubble_blocks via kmer_bubble_id
    # First, accumulate (kmer_id, block_id) pairs efficiently:
    # for each k-mer, all of its bubble's blocks.
    block_kmer_lists: list[list[int]] = [[] for _ in range(n_blocks)]
    for k in range(len(kmer_bubble_id)):
        b_idx = int(kmer_bubble_id[k])
        for blk in bubble_blocks[b_idx]:
            block_kmer_lists[blk].append(k)
    block_kmer_idx = [np.asarray(lst, dtype=np.int64) for lst in block_kmer_lists]
    return block_kmer_idx


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


def find_haplotype_clusters(cn_block, min_cluster_overlap=1):
    """Group founders that share an identical haplotype within a block.

    Args:
        cn_block: F × K_block dense int8/float32 (founders × k-mers in block)
        min_cluster_overlap: minimum non-zero entries to consider a founder
            informative within this block (founders with all-zero rows are
            placed in a single "uninformative" cluster)

    Returns:
        cluster_id: F-vec of cluster index (0..C-1)
        n_clusters: number of distinct clusters
    """
    F = cn_block.shape[0]
    # Hash each founder's row by its bytes
    cn_int = (cn_block > 0).astype(np.uint8)  # binarize
    row_active = cn_int.sum(axis=1) >= min_cluster_overlap
    cluster_id = np.full(F, -1, dtype=np.int32)
    seen = {}
    next_id = 0
    for f in range(F):
        if not row_active[f]:
            continue
        h = cn_int[f].tobytes()
        if h in seen:
            cluster_id[f] = seen[h]
        else:
            seen[h] = next_id
            cluster_id[f] = next_id
            next_id += 1
    # All-zero founders go to a shared uninformative cluster
    if (~row_active).any():
        cluster_id[~row_active] = next_id
        next_id += 1
    return cluster_id, next_id


def solve_em_clustered_per_block(counts, cn_kmer_dense, kmer_block, n_blocks,
                                 coverage, em_max_iter=200, tol=1e-7,
                                 min_kmers_per_block=200, verbose=False):
    """Per-block cluster-EM — fixes rank deficiency in 231-dim direct EM.

    For each block:
      1. Cluster founders by exact-match haplotype within the block
         (typically reduces 231→10–50 clusters per 200 kb in selfing pop)
      2. Run multinomial EM in cluster-frequency space (well-conditioned)
      3. Distribute cluster freqs uniformly within cluster → per-founder h
         (mathematically equivalent for the alt_freq projection: alt_freq(r)
          is invariant to redistributing mass among founders sharing the
          same local haplotype)

    Returns:
        h_blocks_founder: (n_blocks, F) — per-founder h per block (uniform
            within each cluster, by construction; useful for downstream
            comparison with non-clustered results)
        block_status: (n_blocks,) — 0=local, 1=global-fallback, 2=excluded
        global_h: F-vec global fallback (uniform-founder EM)
        cluster_info: list of dicts, one per block, with
            {n_clusters, cluster_id (F-vec), h_cluster (C-vec)}
    """
    F = cn_kmer_dense.shape[0]
    h_blocks = np.zeros((n_blocks, F), dtype=np.float32)
    block_status = np.full(n_blocks, 2, dtype=np.int8)
    cluster_info = [None] * n_blocks

    if verbose:
        print(f"  computing global-h fallback...", flush=True)
    t = time.time()
    global_h, info = solve_em(counts, cn_kmer_dense, coverage,
                              max_iter=em_max_iter, tol=tol)
    global_h = global_h.astype(np.float32)
    if verbose:
        print(f"    global EM: {info['iterations']} iters, {time.time()-t:.0f}s")

    nz = counts > 0
    cluster_sizes_log = []
    n_local = n_global = n_excl = 0
    t0 = time.time()
    for b in range(n_blocks):
        in_block = (kmer_block == b)
        if not in_block.any():
            h_blocks[b] = global_h
            block_status[b] = 2
            n_excl += 1
            continue
        in_block_nz = in_block & nz
        if in_block_nz.sum() < min_kmers_per_block:
            h_blocks[b] = global_h
            block_status[b] = 1
            n_global += 1
            continue
        cn_b = np.ascontiguousarray(cn_kmer_dense[:, in_block_nz])
        c_b = counts[in_block_nz]

        # Cluster founders within the block
        cluster_id, n_clust = find_haplotype_clusters(cn_b)
        cluster_sizes_log.append(n_clust)
        # Build cluster representative cn_kmer (C × K_block)
        # — each cluster's row is the OR (or any representative) of its founders
        # — since exact-match clustering, all members are identical
        cn_cluster = np.zeros((n_clust, cn_b.shape[1]), dtype=np.float32)
        cluster_size = np.zeros(n_clust, dtype=np.int32)
        for f in range(F):
            cid = cluster_id[f]
            if cid < 0:
                continue
            if cluster_size[cid] == 0:
                cn_cluster[cid] = cn_b[f]
            cluster_size[cid] += 1

        # Cluster-EM (well-conditioned: dim = n_clust, typically 10–50)
        h_cluster, _ = solve_em(c_b, cn_cluster, coverage,
                                max_iter=em_max_iter, tol=tol)
        h_cluster = h_cluster.astype(np.float32)

        # Distribute uniformly within cluster → per-founder h for this block
        h_founder = np.zeros(F, dtype=np.float32)
        for f in range(F):
            cid = cluster_id[f]
            if cid >= 0 and cluster_size[cid] > 0:
                h_founder[f] = h_cluster[cid] / cluster_size[cid]

        h_blocks[b] = h_founder
        block_status[b] = 0
        cluster_info[b] = {
            "n_clusters": int(n_clust),
            "cluster_id": cluster_id,
            "h_cluster": h_cluster,
        }
        n_local += 1
    if verbose:
        if cluster_sizes_log:
            print(f"  cluster dim per block: median={int(np.median(cluster_sizes_log))}, "
                  f"min={min(cluster_sizes_log)}, max={max(cluster_sizes_log)} "
                  f"(out of F={F})")
        print(f"  per-block cluster-EM: {n_local} local, {n_global} fallback, "
              f"{n_excl} excluded; {time.time()-t0:.0f}s")
    return h_blocks, block_status, global_h, cluster_info


def solve_em_per_block(counts, cn_kmer_dense, kmer_block, n_blocks,
                       coverage, em_max_iter=200, tol=1e-7,
                       min_kmers_per_block=200, verbose=False,
                       n_workers=4,
                       global_anchor_weight: float = 0.0,
                       block_kmer_idx_override: list | None = None):
    """Run EM independently per block.

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
                              max_iter=em_max_iter, tol=tol)
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

    # Pre-compute per-block masks once. Use override (overlapping windows)
    # if provided; otherwise build from the legacy 1D kmer_block scheme.
    if block_kmer_idx_override is not None:
        block_kmer_idx = block_kmer_idx_override
        if verbose:
            n_kmers_in_some_block = sum(len(idx) for idx in block_kmer_idx)
            print(f"  using overlapping-window block_kmer_idx: "
                  f"{n_kmers_in_some_block:,} (kmer, block) pairs across "
                  f"{n_blocks} blocks "
                  f"(avg {n_kmers_in_some_block/max(1,n_blocks):.0f} kmers/block)",
                  flush=True)
    else:
        block_kmer_idx = [None] * n_blocks
        for b in range(n_blocks):
            block_kmer_idx[b] = np.flatnonzero(kmer_block == b)

    def _fit_one(b):
        idxs = block_kmer_idx[b]
        if len(idxs) == 0:
            return b, global_h, 2
        idxs_nz = idxs[nz[idxs]]
        if len(idxs_nz) < min_kmers_per_block:
            return b, global_h, 1
        cn_b = np.ascontiguousarray(cn_kmer_dense[:, idxs_nz])
        c_b = counts[idxs_nz]
        if global_anchor_weight > 0:
            h_b, _ = solve_em(c_b, cn_b, coverage,
                              max_iter=em_max_iter, tol=tol,
                              prior_h=global_h,
                              prior_weight=global_anchor_weight)
        else:
            h_b, _ = solve_em(c_b, cn_b, coverage,
                              max_iter=em_max_iter, tol=tol)
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


def project_blocks_to_records_overlap(h_blocks, block_status, global_h,
                                       cn_var, blocks,
                                       record_chrom, record_pos,
                                       eps_bp: float = 100.0,
                                       cn_var_called=None):
    """Project per-block h to per-record AFs using INVERSE-DISTANCE AVERAGING
    across all blocks that COVER the record.

    For overlapping windows (Route 2): each record at position p is inside ~3
    windows; compute h_r as the inverse-distance-weighted average of those
    windows' h vectors, then project alt_freq = h_r @ cn_var[:, r].

    Args:
        h_blocks, block_status, global_h: as in project_blocks_to_records
        blocks: list of BlockSpec
        record_chrom, record_pos: per-record metadata
        eps_bp: small offset to inverse-distance weights so a record exactly
            at a window center doesn't blow up; ε=100 bp matches the typical
            window step granularity.
    """
    N = cn_var.shape[1]
    out = np.zeros(N, dtype=np.float64)
    record_chrom = np.asarray(record_chrom)
    record_pos = np.asarray(record_pos).astype(np.int64)

    # Effective per-block h with fallback for excluded/global windows
    F = h_blocks.shape[1]
    h_eff = h_blocks.copy().astype(np.float32)
    fallback_mask = (block_status != 0)
    h_eff[fallback_mask] = global_h.astype(np.float32)

    # Per-chrom: sorted blocks for fast cover-search
    chrom_to_block_idxs: dict[str, list[int]] = {}
    for i, b in enumerate(blocks):
        chrom_to_block_idxs.setdefault(b.chrom, []).append(i)
    block_centers = np.array([(b.start + b.end) / 2.0 for b in blocks])
    block_starts_arr = np.array([b.start for b in blocks], dtype=np.int64)
    block_ends_arr = np.array([b.end for b in blocks], dtype=np.int64)

    cn_var_csc = cn_var.tocsc()
    for chrom, idxs in chrom_to_block_idxs.items():
        idxs_sorted = np.array(sorted(idxs, key=lambda i: blocks[i].start),
                               dtype=np.int64)
        starts_c = block_starts_arr[idxs_sorted]
        ends_c = block_ends_arr[idxs_sorted]
        centers_c = block_centers[idxs_sorted]
        chr_mask = record_chrom == chrom
        if not chr_mask.any():
            continue
        rec_idx = np.flatnonzero(chr_mask)
        positions = record_pos[rec_idx]

        # For each record p, the covering windows are those with start <= p <= end.
        # Since windows are sorted by start and have fixed width W, we can find
        # the leftmost candidate (rightmost block whose start <= p) and walk
        # back/forward while end >= p / start <= p.
        # Vectorized: for each record, find candidate left index via searchsorted
        # on starts_c (rightmost start <= p). Then expand to the cover set.
        n_blocks_c = len(idxs_sorted)

        # Compute the maximum possible window width to know how far back to look
        widths_c = ends_c - starts_c + 1
        max_w = int(widths_c.max()) if len(widths_c) else 0
        # The "step" between adjacent starts gives # of overlapping windows
        if n_blocks_c > 1:
            steps = np.diff(starts_c)
            min_step = int(steps.min())
        else:
            min_step = max_w
        max_overlap = max(1, max_w // max(1, min_step))

        # Per-record cover sets: walk from candidate index +/- max_overlap
        right = np.searchsorted(starts_c, positions, side='right') - 1
        right = np.clip(right, 0, n_blocks_c - 1)

        rec_block_lists: list[np.ndarray] = []
        for ri, p in enumerate(positions):
            r = int(right[ri])
            cover = []
            # walk back: while starts_c[k] <= p and ends_c[k] >= p
            k = r
            while k >= 0 and starts_c[k] <= p:
                if ends_c[k] >= p:
                    cover.append(int(idxs_sorted[k]))
                k -= 1
                if (r - k) > max_overlap + 2:
                    break  # bounded walk
            # walk forward (in case starts_c[r] > p missed some)
            k = r + 1
            while k < n_blocks_c and starts_c[k] <= p:
                if ends_c[k] >= p:
                    cover.append(int(idxs_sorted[k]))
                k += 1
                if (k - r) > max_overlap + 2:
                    break
            if not cover:
                # No block covers this record — fall back to nearest center
                d = np.abs(centers_c - p)
                cover = [int(idxs_sorted[int(np.argmin(d))])]
            rec_block_lists.append(np.asarray(cover, dtype=np.int64))

        # Compute weights and aggregate. MAR projection (2026-05-21): divide
        # by called-mask weighted h to match global-mode semantics. Without
        # the called mask we'd silently treat ./. as REF → SV AF under-call.
        cv = cn_var_csc[:, rec_idx].toarray().astype(np.float32)  # F × n_rec
        if cn_var_called is not None:
            cvc = cn_var_called.tocsc()[:, rec_idx].toarray().astype(np.float32)
        else:
            cvc = None
        for ri, p in enumerate(positions):
            cover = rec_block_lists[ri]
            ctrs = block_centers[cover]
            d = np.abs(ctrs - p) + eps_bp
            w = (1.0 / d).astype(np.float32)
            w /= w.sum()
            h_avg = (h_eff[cover].T @ w).astype(np.float32)
            num = float(h_avg @ cv[:, ri])
            if cvc is not None:
                den = float(h_avg @ cvc[:, ri])
                out[rec_idx[ri]] = num / max(den, 1e-12)
            else:
                out[rec_idx[ri]] = num
    return out


def project_blocks_to_records(h_blocks, block_status, global_h,
                              cn_var, record_block,
                              record_chrom=None, record_pos=None,
                              blocks=None, smooth=True,
                              cn_var_called=None):
    """Project per-block h to per-record alt freqs.

    smooth=False (legacy): each record uses h_blocks[record_block[r]] directly
        — produces hard step discontinuities at block boundaries.

    smooth=True (default, HAFpipe-like): for each record, linearly
        interpolate between the h of the block containing it and the h of
        the geometrically nearest neighbor block on the same chromosome.
        Weighted by inverse distance to block centers. Removes the step
        discontinuities while preserving local mosaic structure.

    For smooth=True, record_chrom, record_pos, and blocks must be provided.

    cn_var_called (MAR projection, 2026-05-21): if provided, divide each
    record's AF by `h_avg @ cn_var_called[:, r]` (h-weighted called mass).
    Matches global-mode semantics. Without it, ./. is silently treated as REF
    → SV AF under-call at high-missingness records.
    """
    N = cn_var.shape[1]
    out = np.zeros(N, dtype=np.float64)

    if not smooth:
        rb = np.asarray(record_block)
        unique_b = np.unique(rb)
        cn_var_csc = cn_var.tocsc()
        cvc_csc = cn_var_called.tocsc() if cn_var_called is not None else None
        for b in unique_b:
            mask = (rb == b)
            if not mask.any():
                continue
            h = h_blocks[b] if b >= 0 else global_h
            cv = cn_var_csc[:, mask]
            num = (h.astype(np.float32) @ cv.toarray()).astype(np.float64)
            if cvc_csc is not None:
                cvc = cvc_csc[:, mask]
                den = (h.astype(np.float32) @ cvc.toarray()).astype(np.float64)
                out[mask] = num / np.maximum(den, 1e-12)
            else:
                out[mask] = num
        return out

    # Smooth: linear interpolation between adjacent blocks on same chrom.
    if record_chrom is None or record_pos is None or blocks is None:
        raise ValueError("smooth=True requires record_chrom, record_pos, blocks")

    record_chrom = np.asarray(record_chrom)
    record_pos = np.asarray(record_pos)

    # Pre-compute h_for_block including the global fallback for missing/excluded
    F = h_blocks.shape[1]
    h_eff = h_blocks.copy().astype(np.float32)
    fallback_mask = (block_status != 0)
    h_eff[fallback_mask] = global_h.astype(np.float32)

    # Per-chromosome processing (within-chrom interpolation only)
    chrom_to_block_idxs = {}
    for i, b in enumerate(blocks):
        chrom_to_block_idxs.setdefault(b.chrom, []).append(i)
    block_centers = np.array([(b.start + b.end) / 2.0 for b in blocks])

    cn_var_csc = cn_var.tocsc()
    cvc_csc = cn_var_called.tocsc() if cn_var_called is not None else None
    for chrom, idxs in chrom_to_block_idxs.items():
        idxs = np.array(sorted(idxs, key=lambda i: blocks[i].start))
        centers_c = block_centers[idxs]
        chr_mask = record_chrom == chrom
        if not chr_mask.any():
            continue
        positions = record_pos[chr_mask]
        # For each record on this chrom, find the two flanking block centers.
        # right_idx in [0..len(idxs)] is where pos would insert
        right_idx = np.searchsorted(centers_c, positions)
        # Left and right block indices (clip at edges)
        right_idx_clip = np.clip(right_idx, 1, len(idxs))
        left_idx_clip = right_idx_clip - 1
        left_centers = centers_c[left_idx_clip]
        right_centers = centers_c[np.clip(right_idx_clip, 0, len(idxs)-1)]
        # When pos is outside both edges, weight collapses to nearest block
        denom = np.maximum(right_centers - left_centers, 1.0)
        w_right = np.clip((positions - left_centers) / denom, 0.0, 1.0)
        w_left = 1.0 - w_right
        # Map back to global block ids
        gb_left = idxs[left_idx_clip]
        gb_right = idxs[np.clip(right_idx_clip, 0, len(idxs)-1)]

        # Subset cn_var to these records
        rec_idx = np.flatnonzero(chr_mask)
        cv = cn_var_csc[:, rec_idx].toarray()  # F × n_rec dense
        # h_left @ cv and h_right @ cv (vectorized over records)
        Hl = h_eff[gb_left]   # (n_rec, F)
        Hr = h_eff[gb_right]  # (n_rec, F)
        # alt_freq[i] = w_left[i] * (Hl[i] @ cv[:,i]) + w_right[i] * (Hr[i] @ cv[:,i])
        # = sum_f cv[f,i] * (w_left[i] * Hl[i,f] + w_right[i] * Hr[i,f])
        # Compute as einsum
        af_num = (w_left * np.einsum("if,fi->i", Hl, cv)
                  + w_right * np.einsum("if,fi->i", Hr, cv))
        if cvc_csc is not None:
            # MAR projection: divide by interpolated h-weighted called mass.
            cvc = cvc_csc[:, rec_idx].toarray()
            af_den = (w_left * np.einsum("if,fi->i", Hl, cvc)
                      + w_right * np.einsum("if,fi->i", Hr, cvc))
            out[rec_idx] = af_num / np.maximum(af_den, 1e-12)
        else:
            out[rec_idx] = af_num
    return out
