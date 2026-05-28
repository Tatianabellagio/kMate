"""
LD-defined blocks for per-block EM, replacing fixed-bp windows.

Algorithm (hapFIRE's CompleteLDPartition):
    For each variant i, find the rightmost variant j (within a search window)
    such that r²(i, j) > threshold. Compute right_max[i] = cummax of these.
    Blocks are split at indices i where right_max[i-1] < i — meaning no
    earlier variant has any LD partner reaching variant i.

This produces *fully-independent* blocks: no variant inside a block has any
LD partner outside the block. Block sizes are variable, adapted to the
recombination landscape and selfing-induced LD structure.

Reference: Czech, Spence, Bellagio, Exposito-Alonso et al. 2022 bioRxiv,
hapFIRE/haplotype_generation.py:CompleteLDPartition.

For Arabidopsis (selfing ~97%, LD spans 10–100 kb between related accessions)
the typical block size at r²=0.1 is 1–10 kb of contiguous variants — much
finer than 200 kb fixed windows.
"""
from __future__ import annotations
import time
from dataclasses import dataclass
from typing import Iterable

import numpy as np
from scipy.sparse import issparse

from block_em import BlockSpec


def compute_ld_blocks_gabriel(
    cn_var,                    # F × N sparse 0/1
    record_chrom,
    record_pos,
    r2_threshold: float = 0.5,
    smooth_window: int = 50,
    min_block_records: int = 100,
    max_block_bp: int = 500_000,
    verbose: bool = True,
) -> list[BlockSpec]:
    """LD-block boundaries from Gabriel-style adjacent-r² drops.

    For each pair of adjacent variants (i, i+1) on a chromosome, compute the
    mean r² between record i and the prior `smooth_window` records (the
    "current block centroid"). When this drops below r2_threshold, a new
    block begins.

    Compared to CompleteLDPartition, this captures *local* LD breaks rather
    than chromosome-wide reach. For highly-inbred Arabidopsis (where every
    record has SOME long-range LD partner, making CompleteLDPartition
    degenerate), the Gabriel-style criterion gives 100s of blocks per chrom
    at biologically relevant scales.

    Soft cap at max_block_bp prevents runaway blocks across centromeres.
    """
    record_chrom = np.asarray(record_chrom)
    record_pos = np.asarray(record_pos)
    F, N = cn_var.shape
    blocks: list[BlockSpec] = []

    for chrom in sorted(set(record_chrom)):
        m = record_chrom == chrom
        idx_chrom = np.flatnonzero(m)
        if len(idx_chrom) == 0:
            continue
        order = np.argsort(record_pos[idx_chrom])
        idx_chrom = idx_chrom[order]
        positions = record_pos[idx_chrom]
        N_c = len(idx_chrom)

        if verbose:
            print(f"  {chrom}: {N_c:,} records ...", flush=True)
        t = time.time()

        G = cn_var[:, idx_chrom]
        if issparse(G):
            G = G.toarray()
        G = G.astype(np.float32)
        G -= G.mean(axis=0, keepdims=True)
        norm = np.linalg.norm(G, axis=0)
        nz_var = norm > 1e-9
        G /= np.maximum(norm, 1e-12)

        # For each i ≥ 1, compute mean |r| between G[:,i] and G[:,i-W..i-1]
        # Done in batches for memory.
        block_start = 0
        chrom_blocks_pos = []  # list of (i_start, i_end) record-index ranges
        BATCH = 8192
        prev_window_end = 0
        for b0 in range(0, N_c, BATCH):
            b1 = min(b0 + BATCH, N_c)
            cur = G[:, b0:b1]                          # F × batch
            ws = max(0, b0 - smooth_window)
            ctx = G[:, ws:b1]                          # F × (batch + W or less)
            r = ctx.T @ cur                            # (batch+W) × batch
            absr = np.abs(r)
            # For each variant i in [b0, b1), avg |r| over its prior W records
            for k in range(b1 - b0):
                i = b0 + k
                if i == 0:
                    continue  # first record always starts block 0
                # context is G[:, max(0, i-W):i]  → indices in `ctx` are
                # [max(0, i - W) - ws .. i - ws)
                left = max(0, i - smooth_window)
                col_in_ctx_lo = left - ws
                col_in_ctx_hi = i - ws
                if col_in_ctx_hi <= col_in_ctx_lo:
                    continue
                # current variant is column `k` in `cur`, also column (i - ws) in ctx
                # We want absr[col_in_ctx_lo:col_in_ctx_hi, k]
                vals = absr[col_in_ctx_lo:col_in_ctx_hi, k]
                # If current variant has zero variance, skip (treat as continuation)
                if not nz_var[i]:
                    continue
                # Take average; LD threshold is on mean |r| → equivalent r²
                mean_absr = float(vals.mean())
                if (mean_absr * mean_absr) < r2_threshold or \
                   (positions[i] - positions[block_start] > max_block_bp):
                    # close prior block [block_start, i)
                    chrom_blocks_pos.append((block_start, i))
                    block_start = i
        # final block
        if block_start < N_c:
            chrom_blocks_pos.append((block_start, N_c))

        # Merge tiny blocks
        merged_record_blocks = []
        for (i0, i1) in chrom_blocks_pos:
            if (i1 - i0) < min_block_records and merged_record_blocks:
                p_i0, _ = merged_record_blocks[-1]
                merged_record_blocks[-1] = (p_i0, i1)
            else:
                merged_record_blocks.append((i0, i1))

        for (i0, i1) in merged_record_blocks:
            blocks.append(BlockSpec(
                chrom=str(chrom),
                start=int(positions[i0]),
                end=int(positions[i1 - 1]),
            ))

        if verbose:
            sizes_bp = np.array([blocks[-1].end - blocks[-1].start + 1
                                 for _ in range(1)])  # last only
            sizes = [b.end - b.start + 1
                     for b in blocks if b.chrom == str(chrom)]
            print(f"    {len(merged_record_blocks)} blocks; "
                  f"sizes bp: median={int(np.median(sizes)):,}, "
                  f"max={max(sizes):,}; [{time.time()-t:.0f}s]", flush=True)

    return blocks


def compute_ld_blocks(
    cn_var,                 # F × N sparse 0/1 (founder genotype at biallelic record)
    record_chrom,           # N-vec
    record_pos,             # N-vec
    r2_threshold: float = 0.1,
    search_window: int = 100,
    min_block_records: int = 50,
    verbose: bool = True,
) -> list[BlockSpec]:
    """Compute LD-coherent blocks via CompleteLDPartition.

    For each variant i, scan up to `search_window` neighbors and record the
    rightmost variant j with r²(i, j) > r2_threshold. Block boundaries are
    placed where this "right reach" cummax breaks (i.e., the current variant
    has no LD partner anywhere to the left or right beyond itself).

    Tiny blocks (< min_block_records) are merged into the previous block.
    """
    record_chrom = np.asarray(record_chrom)
    record_pos = np.asarray(record_pos)
    F, N = cn_var.shape
    blocks: list[BlockSpec] = []

    for chrom in sorted(set(record_chrom)):
        m = record_chrom == chrom
        idx_chrom = np.flatnonzero(m)
        if len(idx_chrom) == 0:
            continue
        # Sort by position (records should already be sorted, but ensure)
        order = np.argsort(record_pos[idx_chrom])
        idx_chrom = idx_chrom[order]
        positions = record_pos[idx_chrom]
        N_c = len(idx_chrom)

        if verbose:
            print(f"  {chrom}: {N_c:,} records, computing LD ...", flush=True)
        t = time.time()

        # Pull dense F × N_c chunk and standardize for fast r² via dot products
        G = cn_var[:, idx_chrom]
        if issparse(G):
            G = G.toarray()
        G = G.astype(np.float32)
        # Subtract per-column mean and divide by per-column L2 → r = G[:,i]·G[:,j]
        G -= G.mean(axis=0, keepdims=True)
        norm = np.linalg.norm(G, axis=0)
        G /= np.maximum(norm, 1e-12)
        # Records with zero variance (all-0 or all-1 rows) → norm=0 → division
        # gave 0; their dot products with anything are 0; r²=0 → never a partner.

        # right_partner[i] = max j in [i+1, i+W) with r²(i,j) > threshold
        right_partner = np.arange(N_c, dtype=np.int64)
        thresh = np.float32(np.sqrt(r2_threshold))  # |r| ≥ √r² gives r² ≥ thresh

        # Vectorize per-variant: compute G[:,i]^T @ G[:,i+1:i+W]
        # Process in mini-batches to amortize memory
        BATCH = 4096
        for b0 in range(0, N_c, BATCH):
            b1 = min(b0 + BATCH, N_c)
            # Right-neighbors window for batch [b0, b1):
            # j ranges in [b0+1, b1+search_window). Compute full
            # cross-product G[:, b0:b1].T @ G[:, b0:b1+search_window]
            j_end = min(b1 + search_window, N_c)
            left = G[:, b0:b1]              # F × batch
            right = G[:, b0:j_end]          # F × (batch + W)
            r = left.T @ right              # batch × (batch + W)
            absr = np.abs(r)
            # For each variant i in [b0, b1), find max j in (i, i+W] with |r|>=thresh
            for k in range(b1 - b0):
                i = b0 + k
                jstart = k + 1
                jend = min(k + 1 + search_window, j_end - b0)
                row = absr[k, jstart:jend]
                # rightmost above threshold
                idx_above = np.where(row >= thresh)[0]
                if len(idx_above) > 0:
                    right_partner[i] = b0 + jstart + int(idx_above.max())

        # cummax(right_partner) gives the maximum reach of any variant ≤ i
        right_max = np.maximum.accumulate(right_partner)

        # Block boundaries: i where right_max[i-1] < i (variant i is unreachable
        # from any earlier variant)
        boundary_idxs = [0]
        for i in range(1, N_c):
            if right_max[i - 1] < i:
                boundary_idxs.append(i)
        boundary_idxs.append(N_c)

        # Build BlockSpecs and merge tiny ones
        chrom_blocks = []
        for k in range(len(boundary_idxs) - 1):
            i0, i1 = boundary_idxs[k], boundary_idxs[k + 1]
            chrom_blocks.append(BlockSpec(
                chrom=str(chrom),
                start=int(positions[i0]),
                end=int(positions[i1 - 1]),
            ))
        # Merge tiny blocks into prior
        merged: list[BlockSpec] = []
        carry = None
        for k in range(len(boundary_idxs) - 1):
            i0, i1 = boundary_idxs[k], boundary_idxs[k + 1]
            n_recs = i1 - i0
            if n_recs < min_block_records and merged:
                # Extend prior block's end
                merged[-1] = BlockSpec(
                    chrom=merged[-1].chrom,
                    start=merged[-1].start,
                    end=int(positions[i1 - 1]),
                )
            elif n_recs < min_block_records:
                # First block tiny — keep, will absorb next
                merged.append(BlockSpec(
                    chrom=chrom_blocks[k].chrom,
                    start=chrom_blocks[k].start,
                    end=chrom_blocks[k].end,
                ))
            else:
                merged.append(chrom_blocks[k])

        if verbose:
            sizes_bp = np.array([b.end - b.start + 1 for b in merged])
            print(f"    {len(merged)} blocks (raw {len(chrom_blocks)}); "
                  f"sizes bp: median={int(np.median(sizes_bp)):,}, "
                  f"max={sizes_bp.max():,}, [{time.time()-t:.0f}s]",
                  flush=True)

        blocks.extend(merged)
    return blocks
