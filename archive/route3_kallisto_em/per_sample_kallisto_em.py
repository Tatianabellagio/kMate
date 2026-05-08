"""Kallisto-style per-window pseudoalignment EM for founder pool-seq frequencies.

Architecture (per Bray et al. 2016 kallisto, applied to founder panels):

1. INDEX: cn_full[f, k] (founder × k-mer copy-number) is bit-packed once. Each
   k-mer column becomes a 231-bit founder compatibility set.
2. PSEUDOALIGN per read: canonical 31-merize the read, look up each k-mer's
   founder set, and INTERSECT (bitwise AND) the sets within each 10 kb window
   the read touches. Output: per (read, window) "founder compatibility class"
   S — the set of founders that could explain ALL the read's k-mers in that
   window.
3. EQUIVALENCE-CLASS COUNTS: hash (window_id, S) and accumulate counts. Reads
   that produce identical compatibility sets in the same window collapse to
   one EC. The EC counts are the sufficient statistic for the EM (this is the
   kallisto trick — per-iter cost is O(unique_ec) not O(reads)).
4. PER-WINDOW EM on EC counts:
        h_new[f] ∝ h[f] · Σ_c n_c · I[f ∈ S_c] / (Σ_{f' ∈ S_c} h[f'])
   Same multiplicative simplex update as cactus_em, but with EC observations
   instead of per-k-mer Poisson counts. This recovers the within-read phase
   information that per-k-mer counting throws away.
5. PROJECT: per-record AF = h_window @ cn_var[:, record], same as cactus_em.

Compared to cactus_em --block-mode window:
  - Same cn_full index. Same cn_var projection. Same windows.
  - DIFFERENT EM: counts EC observations (read-level, phase-aware), not k-mer
    counts (per-k-mer, marginalized).

This is a prototype: pure-Python read parsing, in-memory bit-packed cn_full,
sequential per-window EM. Optimization goal AFTER verifying the approach
beats window_10kb on the n50_g3 recomb sim.
"""
from __future__ import annotations
import argparse
import gc
import os
import sys
import time
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy.sparse import load_npz

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from block_em import define_windows


_COMPLEMENT_TBL = bytes.maketrans(b"ACGTNacgtn", b"TGCANtgcan")

# 2-bit encoding for ACGT. N maps to 0xFF as a sentinel; any k-mer containing
# 0xFF is rejected (skip).
_BASE2BIT = np.full(256, 0xFF, dtype=np.uint8)
for c, v in zip(b"ACGTacgt", [0, 1, 2, 3, 0, 1, 2, 3]):
    _BASE2BIT[c] = v


def revcomp_bytes(s: bytes) -> bytes:
    return s.translate(_COMPLEMENT_TBL)[::-1]


def encode_kmer_u64(kmer: bytes, k: int = 31) -> int:
    """Encode a k-mer string into a uint64 via 2-bit packing.
    Returns -1 (0xFFFFFFFFFFFFFFFF) if any base is N or invalid.
    """
    val = 0
    for c in kmer:
        b = _BASE2BIT[c]
        if b == 0xFF:
            return np.uint64(0xFFFFFFFFFFFFFFFF)
        val = (val << 2) | int(b)
    return np.uint64(val)


def encode_panel_kmers_u64(kmer_strings: np.ndarray, k: int = 31):
    """Vectorized encoding of N panel k-mer strings to a (N,) uint64 array.
    K-mers containing N (or any non-ACGT) get 0xFFFFFFFFFFFFFFFF as a sentinel.
    """
    N = len(kmer_strings)
    out = np.zeros(N, dtype=np.uint64)
    for i in range(N):
        s = kmer_strings[i]
        if isinstance(s, str):
            s = s.encode("ascii")
        out[i] = encode_kmer_u64(s, k)
    return out


def canonical_kmer_u64_iter(seq: bytes, k: int = 31):
    """Yield canonical (min of fwd/revcomp) uint64 encodings of k-mers in seq.

    Uses a rolling 2-bit buffer for both strands. Skips windows containing N
    or any non-ACGT base. Two integers maintained: fwd, rev. Canonical = min.
    """
    L = len(seq)
    if L < k:
        return
    mask = (np.uint64(1) << np.uint64(2 * k)) - np.uint64(1)
    rc_shift = np.uint64(2 * (k - 1))
    fwd = np.uint64(0)
    rev = np.uint64(0)
    valid = 0  # bases since last N
    for i in range(L):
        b = _BASE2BIT[seq[i]]
        if b == 0xFF:
            valid = 0
            fwd = np.uint64(0)
            rev = np.uint64(0)
            continue
        b_u = np.uint64(b)
        b_rc = np.uint64(3 - b)
        fwd = ((fwd << np.uint64(2)) | b_u) & mask
        rev = (rev >> np.uint64(2)) | (b_rc << rc_shift)
        valid += 1
        if valid >= k:
            yield fwd if fwd < rev else rev


def pack_cn_full_as_bitmasks(cn_full, n_founders: int):
    """Convert cn_full (F × K CSR sparse) to a packed (n_words, K) uint64 array.

    Each k-mer column becomes a bitmask: bit f set iff founder f carries the k-mer.
    n_words = ceil(F / 64). For F=231, n_words=4.

    Iterates per-founder (231 rows) rather than per-k-mer (~20M cols) — 100×
    faster. Vectorized OR via numpy fancy indexing on column positions.
    """
    n_words = (n_founders + 63) // 64
    K = cn_full.shape[1]
    packed = np.zeros((n_words, K), dtype=np.uint64)
    cn_csr = cn_full.tocsr()
    for f in range(n_founders):
        cols = cn_csr.getrow(f).indices
        if len(cols) == 0:
            continue
        w = f >> 6
        b = f & 0x3F
        bit = np.uint64(1) << np.uint64(b)
        packed[w, cols] |= bit
    return packed


def all_ones_mask(n_founders: int) -> np.ndarray:
    """Bitmask with bits 0..n_founders-1 set, packed into uint64 words."""
    n_words = (n_founders + 63) // 64
    mask = np.zeros(n_words, dtype=np.uint64)
    for f in range(n_founders):
        w = f >> 6
        b = f & 0x3F
        mask[w] |= np.uint64(1) << np.uint64(b)
    return mask


def mask_members(mask: np.ndarray, n_founders: int) -> np.ndarray:
    """Return indices of founders set in the bitmask."""
    out = []
    for f in range(n_founders):
        w = f >> 6
        b = f & 0x3F
        if mask[w] & (np.uint64(1) << np.uint64(b)):
            out.append(f)
    return np.asarray(out, dtype=np.int32)


def stream_fastq_seqs(fq_path: str, max_reads: int | None = None):
    """Yield read sequences (bytes) from a FASTQ file. Supports plain or .gz."""
    import gzip
    opener = gzip.open if fq_path.endswith(".gz") else open
    n = 0
    with opener(fq_path, "rb") as fh:
        line_idx = 0
        for line in fh:
            if line_idx % 4 == 1:
                yield line.rstrip(b"\n\r")
                n += 1
                if max_reads is not None and n >= max_reads:
                    return
            line_idx += 1


def lookup_kmers_u64(query_u64: np.ndarray,
                      sorted_panel_u64: np.ndarray,
                      sort_idx: np.ndarray) -> np.ndarray:
    """Vectorized lookup: for each query k-mer (uint64), return its column id
    in the original panel layout, or -1 if not in panel.
    """
    pos = np.searchsorted(sorted_panel_u64, query_u64)
    pos_clip = np.clip(pos, 0, len(sorted_panel_u64) - 1)
    found = (pos < len(sorted_panel_u64)) & (sorted_panel_u64[pos_clip] == query_u64)
    out = np.where(found, sort_idx[pos_clip], -1).astype(np.int64)
    return out


def encode_reads_batch_u64(seqs: list[bytes], k: int = 31):
    """Vectorized canonical k-mer encoding for a batch of equal-length reads.

    Returns:
        kmers_u64: (n_reads, n_kmers_per_read) uint64 — canonical (min of fwd, rev)
        valid:     (n_reads, n_kmers_per_read) bool — False if k-mer contains N
        n_per_row: int — n_kmers_per_read = L - k + 1 (assumes all reads same length)

    Reads not of length L are skipped (caller can filter or pad).
    """
    if not seqs:
        return (np.empty((0, 0), dtype=np.uint64),
                np.empty((0, 0), dtype=bool), 0)
    L = len(seqs[0])
    # Filter to reads of consistent length
    seqs = [s for s in seqs if len(s) == L]
    if not seqs:
        return (np.empty((0, 0), dtype=np.uint64),
                np.empty((0, 0), dtype=bool), 0)
    n_reads = len(seqs)
    arr = np.frombuffer(b"".join(seqs), dtype=np.uint8).reshape(n_reads, L)
    codes = _BASE2BIT[arr]                            # (N, L) uint8, 0xFF for N
    valid_base = codes != 0xFF                        # (N, L)
    # codes_safe: replace 0xFF with 0 to avoid corrupting the shift; we mask later
    codes_safe = np.where(valid_base, codes, 0).astype(np.uint64)
    n_kmers = L - k + 1
    if n_kmers <= 0:
        return (np.empty((0, 0), dtype=np.uint64),
                np.empty((0, 0), dtype=bool), 0)

    # Forward: fwd[r, i] = packed bits codes[r, i..i+k-1]
    fwd = np.zeros((n_reads, n_kmers), dtype=np.uint64)
    for j in range(k):
        fwd = (fwd << np.uint64(2)) | codes_safe[:, j:j + n_kmers]
    # Reverse-complement: rev[r, i] = revcomp of codes[r, i..i+k-1]
    rev = np.zeros((n_reads, n_kmers), dtype=np.uint64)
    # rev_codes_safe: complement = 3 - code (for valid bases)
    cmp_safe = (np.uint64(3) - codes_safe) & np.uint64(0b11)
    for j in range(k):
        # bit position in rev for original index j: shifted left by 2*(k-1-j)
        rev |= cmp_safe[:, j:j + n_kmers] << np.uint64(2 * (k - 1 - j))
    canon = np.minimum(fwd, rev)

    # Validity: a k-mer at pos i is valid iff all bases [i, i+k-1] are valid.
    # Use a sliding-window AND via cumulative count of invalid bases.
    invalid = (~valid_base).astype(np.int32)   # 1 where N
    csum = np.zeros((n_reads, L + 1), dtype=np.int32)
    np.cumsum(invalid, axis=1, out=csum[:, 1:])
    n_invalid_in_window = csum[:, k:] - csum[:, :L - k + 1]
    valid_kmer = n_invalid_in_window == 0
    return canon, valid_kmer, n_kmers


def pseudoalign_reads_to_ec(
    reads_paths: list[str],
    sorted_panel_u64: np.ndarray,
    sort_idx: np.ndarray,
    kmer_to_window: np.ndarray,
    packed_cn: np.ndarray,
    n_founders: int,
    n_windows: int,
    k: int = 31,
    max_reads_per_file: int | None = None,
    log_every: int = 1_000_000,
    batch_size: int = 50_000,
):
    """Pseudoalign reads → per-(window, founder-mask) equivalence-class counts.

    Vectorized batch processing: encode many reads at once, look up all their
    k-mers in one searchsorted call, then group/intersect per (read, window).

    Returns:
        ec_counts: dict {(window_id, mask_bytes): count}
        n_reads_used: int (reads with at least one informative k-mer)
        n_reads_total: int
    """
    ec_counts: dict[tuple[int, bytes], int] = defaultdict(int)
    n_reads_total = 0
    n_reads_used = 0
    t0 = time.time()

    def flush_batch(seqs):
        nonlocal n_reads_total, n_reads_used
        canon, valid, n_per = encode_reads_batch_u64(seqs, k=k)
        if canon.size == 0:
            n_reads_total += len(seqs)
            return
        n_reads_in_batch = canon.shape[0]
        n_reads_total += len(seqs)
        # Lookup all k-mers in panel at once
        flat = canon.ravel()
        kid_flat = lookup_kmers_u64(flat, sorted_panel_u64, sort_idx)
        kid = kid_flat.reshape(n_reads_in_batch, n_per)
        # Mark invalid k-mers as -1
        kid[~valid] = -1
        # Window assignment per k-mer
        win = np.where(kid >= 0, kmer_to_window[np.maximum(kid, 0)], -1)
        # Process per-read
        for r in range(n_reads_in_batch):
            kids_r = kid[r]
            wins_r = win[r]
            mask_keep = wins_r >= 0
            if not mask_keep.any():
                continue
            kids_v = kids_r[mask_keep]
            wins_v = wins_r[mask_keep]
            order = np.argsort(wins_v, kind="stable")
            wins_s = wins_v[order]
            kids_s = kids_v[order]
            cps = np.concatenate([[0],
                np.where(np.diff(wins_s) != 0)[0] + 1,
                [len(wins_s)]])
            informative = False
            for j in range(len(cps) - 1):
                s, e = cps[j], cps[j + 1]
                w_id = int(wins_s[s])
                cols = packed_cn[:, kids_s[s:e]]
                m = np.bitwise_and.reduce(cols, axis=1)
                if not m.any():
                    continue
                ec_counts[(w_id, m.tobytes())] += 1
                informative = True
            if informative:
                n_reads_used += 1

    for fq_path in reads_paths:
        batch: list[bytes] = []
        for seq in stream_fastq_seqs(fq_path, max_reads=max_reads_per_file):
            batch.append(bytes(seq))
            if len(batch) >= batch_size:
                flush_batch(batch)
                batch = []
                if n_reads_total > 0 and n_reads_total % log_every < batch_size:
                    el = time.time() - t0
                    print(f"  pseudoalign: {n_reads_total:,} reads "
                          f"({n_reads_used:,} informative, "
                          f"{len(ec_counts):,} EC) "
                          f"[{el:.0f}s, "
                          f"{n_reads_total/max(1, el):.0f} reads/s]",
                          flush=True)
        if batch:
            flush_batch(batch)
    return ec_counts, n_reads_used, n_reads_total


def em_one_simplex(
    M: np.ndarray,    # (F, n_ec) membership
    nc: np.ndarray,   # (n_ec,) counts
    n_founders: int,
    h_init: np.ndarray | None = None,
    prior_h: np.ndarray | None = None,
    prior_weight: float = 0.0,
    max_iter: int = 200,
    tol: float = 1e-6,
):
    """One simplex EM solve. Optional Dirichlet-style prior centered on prior_h.

    h_new[f] ∝ h[f] * (M @ cw) + prior_weight * N * prior_h[f]
              normalized to simplex.

    prior_weight=0 → pure MLE. prior_weight=1 → equal influence of prior and data
    when N total observations.
    """
    h = (h_init.copy() if h_init is not None
         else np.full(n_founders, 1.0 / n_founders, dtype=np.float32))
    N = float(nc.sum())
    if N <= 0:
        return h, 0
    pseudo = (None if prior_weight <= 0 or prior_h is None
              else (np.float32(prior_weight) * np.float32(N)
                    * prior_h.astype(np.float32)))
    last_delta = 1.0
    for it in range(max_iter):
        denom = np.maximum(h @ M, np.float32(1e-12))
        cw = nc / denom
        em_term = h * (M @ cw)
        if pseudo is not None:
            em_term = em_term + pseudo
        s = em_term.sum()
        h_new = em_term / max(s, 1e-30)
        last_delta = float(np.linalg.norm(h_new - h))
        h = h_new
        if last_delta < tol:
            break
    return h, it + 1


def build_window_M_nc(ecs: list[tuple[bytes, int]],
                       n_founders: int,
                       n_words: int) -> tuple[np.ndarray, np.ndarray]:
    """Build (F, n_ec) membership matrix M and (n_ec,) counts vector nc."""
    n_ec = len(ecs)
    M = np.zeros((n_founders, n_ec), dtype=np.float32)
    nc = np.zeros(n_ec, dtype=np.float32)
    for c, (mask_words, count) in enumerate(ecs):
        mask_arr = np.frombuffer(mask_words, dtype=np.uint64).copy()
        members = mask_members(mask_arr, n_founders)
        if members.size:
            M[members, c] = 1.0
        nc[c] = float(count)
    return M, nc


def em_solve_per_window_ec(
    ec_by_window: dict[int, list[tuple[bytes, int]]],
    n_windows: int,
    n_founders: int,
    max_iter: int = 200,
    tol: float = 1e-6,
    global_anchor_weight: float = 0.0,
):
    """Run kallisto-style EM per window on EC counts.

    If global_anchor_weight > 0, first solve a global EM over the union of all
    ECs (window-blind) to get h_global, then re-solve per window with a
    Dirichlet pseudocount toward h_global.

    Returns:
        h_blocks:     (n_windows, n_founders) float32
        block_status: (n_windows,) int — 0 = solved, 1 = empty
        h_global:     (n_founders,) float32 — chrom-wide EM result (uniform if
                      global_anchor_weight == 0)
    """
    n_words = (n_founders + 63) // 64
    h_blocks = np.zeros((n_windows, n_founders), dtype=np.float32)
    status = np.full(n_windows, 1, dtype=np.int8)
    fallback_h = np.full(n_founders, 1.0 / n_founders, dtype=np.float32)

    # Optional: solve a global EM first by collapsing all ECs into one big problem.
    # Note: a single EC's compatibility set is window-specific, but founder
    # membership is window-agnostic (cn_full[:, k] doesn't depend on the window).
    # So the global solve is just an EM over all ECs treated as a single mixture.
    if global_anchor_weight > 0:
        all_ecs: list[tuple[bytes, int]] = []
        for ecs in ec_by_window.values():
            all_ecs.extend(ecs)
        if all_ecs:
            print(f"  [anchor] solving global EM over {len(all_ecs):,} ECs...",
                  flush=True)
            t = time.time()
            M_all, nc_all = build_window_M_nc(all_ecs, n_founders, n_words)
            h_global, n_it = em_one_simplex(
                M_all, nc_all, n_founders,
                max_iter=max_iter, tol=tol)
            print(f"  [anchor] global EM done ({n_it} iter, "
                  f"{time.time()-t:.0f}s); top h: {np.argsort(-h_global)[:5]}",
                  flush=True)
            del M_all, nc_all
        else:
            h_global = fallback_h.copy()
    else:
        h_global = fallback_h.copy()

    for w, ecs in ec_by_window.items():
        if not ecs:
            h_blocks[w] = h_global if global_anchor_weight > 0 else fallback_h
            continue
        M, nc = build_window_M_nc(ecs, n_founders, n_words)
        if nc.sum() == 0:
            h_blocks[w] = h_global if global_anchor_weight > 0 else fallback_h
            continue
        h_w, _ = em_one_simplex(
            M, nc, n_founders,
            prior_h=h_global if global_anchor_weight > 0 else None,
            prior_weight=global_anchor_weight,
            max_iter=max_iter, tol=tol)
        h_blocks[w] = h_w
        status[w] = 0
    return h_blocks, status, h_global


def project_blocks_to_records_simple(
    h_blocks: np.ndarray,
    cn_var,
    record_block: np.ndarray,
    fallback_h: np.ndarray,
    chunk: int = 50_000,
):
    """Per-record alt freq using the h of the record's window (no smoothing for prototype).

    Chunked: avoids materializing very large dense slices for windows
    containing many records (the b=-1 fallback bucket can hold millions of
    records on multi-chrom panels and would otherwise allocate GB-scale
    dense matrices).
    """
    N = cn_var.shape[1]
    out = np.zeros(N, dtype=np.float64)
    cv_csc = cn_var.tocsc()
    unique_b = np.unique(record_block)
    for b in unique_b:
        rec_idx = np.where(record_block == b)[0]
        if rec_idx.size == 0:
            continue
        h = (h_blocks[b] if b >= 0 else fallback_h).astype(np.float32)
        for s in range(0, rec_idx.size, chunk):
            sub = rec_idx[s:s + chunk]
            cv = cv_csc[:, sub]
            out[sub] = (h @ cv.toarray()).astype(np.float64)
    return out


def assign_kmers_to_windows(
    bubble_id_per_kmer: np.ndarray,
    bubble_chrom: np.ndarray,
    bubble_start: np.ndarray,
    bubble_end: np.ndarray,
    blocks,
    chrom_filter: str | None = None,
):
    """Assign each k-mer to its containing window (-1 if no window)."""
    bubble_centroid = (bubble_start + bubble_end) // 2
    chrom_blocks: dict[str, list[tuple[int, int, int]]] = {}
    for i, b in enumerate(blocks):
        chrom_blocks.setdefault(b.chrom, []).append((b.start, b.end, i))
    for c in chrom_blocks:
        chrom_blocks[c].sort()
    bubble_window = np.full(len(bubble_chrom), -1, dtype=np.int32)
    for b_idx in range(len(bubble_chrom)):
        c = str(bubble_chrom[b_idx])
        if chrom_filter is not None and c != chrom_filter:
            continue
        p = int(bubble_centroid[b_idx])
        if c not in chrom_blocks:
            continue
        # Linear scan — fine for ~3k windows; binary search would be faster for many calls
        for s, e, idx in chrom_blocks[c]:
            if s <= p <= e:
                bubble_window[b_idx] = idx
                break
    kmer_window = bubble_window[bubble_id_per_kmer]
    return kmer_window


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cn-kmer-prefix", required=True,
                    help="prefix for cn_full (e.g. data/cn_full_231_v2/cn)")
    ap.add_argument("--cn-var", required=True)
    ap.add_argument("--cn-var-meta", required=True)
    ap.add_argument("--reads", required=True, nargs="+")
    ap.add_argument("--sample", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--chroms", nargs="+", default=["Chr1"])
    ap.add_argument("--window-bp", type=int, default=10_000)
    ap.add_argument("--em-max-iter", type=int, default=200)
    ap.add_argument("--max-reads-per-file", type=int, default=None,
                    help="Cap reads per FASTQ file (for prototyping/smoke)")
    ap.add_argument("--save-h-blocks", default=None)
    ap.add_argument("--save-ec-cache", default=None,
                    help="optional: cache pseudoalignment EC counts to npz so "
                         "anchor-weight sweeps can skip re-pseudoaligning")
    ap.add_argument("--load-ec-cache", default=None,
                    help="optional: load EC counts from a previous --save-ec-cache "
                         "run; skips read processing entirely")
    ap.add_argument("--global-anchor-weight", type=float, default=0.0,
                    help="KL prior weight λ toward the chrom-wide EM solution. "
                         "0 = pure per-window MLE (default). "
                         "0.05-0.3 = mild anchor, helps low-evidence windows. "
                         "1.0 = anchor as influential as data.")
    args = ap.parse_args()

    print(f"=== {args.sample} (kallisto-EM mode, window={args.window_bp} bp) ===",
          flush=True)
    t0 = time.time()

    # 1. Load cn_full + bubble metadata + uint64 panel cache (no string array)
    print("loading cn_full + bubble metadata + uint64 panel cache...", flush=True)
    import resource
    def _mem():
        return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6
    cn_blocks, bubble_ids = [], []
    bubble_chrom_l, bubble_start_l, bubble_end_l = [], [], []
    sorted_panel_blocks, sort_idx_blocks = [], []
    founders = None
    bubble_offset = 0
    kmer_offset = 0
    for c in args.chroms:
        cn_path = args.cn_kmer_prefix + f"_{c}.cn.npz"
        meta_path = args.cn_kmer_prefix + f"_{c}.meta.npz"
        u64_path = args.cn_kmer_prefix + f"_{c}_kmers_u64.npz"
        if not os.path.exists(cn_path):
            print(f"  WARN: {cn_path} missing, skipping", flush=True)
            continue
        if not os.path.exists(u64_path):
            print(f"  ERROR: {u64_path} missing — run "
                  f"scripts/precompute_panel_u64_index.py first", flush=True)
            sys.exit(1)
        cn_b = load_npz(cn_path)
        # Load only bubble fields from meta — skip the 2 GB kmer_index array
        meta = np.load(meta_path, allow_pickle=True)
        cn_blocks.append(cn_b)
        bubble_ids.append(meta["bubble_id"] + bubble_offset)
        bubble_chrom_l.append(meta["bubble_chrom"])
        bubble_start_l.append(meta["bubble_start"])
        bubble_end_l.append(meta["bubble_end"])
        bubble_offset += int(meta["bubble_id"].max()) + 1
        if founders is None:
            founders = meta["founders"]
        del meta
        u64 = np.load(u64_path)
        sorted_panel_blocks.append(u64["sorted_u64"])
        # sort_idx is original column ID per chrom — offset by kmer_offset for
        # the genome-wide concatenated layout
        sort_idx_blocks.append(u64["sort_idx"] + kmer_offset)
        kmer_offset += cn_b.shape[1]
        del u64
        print(f"  [{c}] cn_b loaded, RSS={_mem():.2f} GB", flush=True)

    if len(cn_blocks) == 1:
        cn_full = cn_blocks[0]   # avoid scipy.hstack doubling memory
    else:
        from scipy.sparse import hstack as sp_hstack
        cn_full = sp_hstack(cn_blocks, format="csr")
    cn_blocks.clear()
    print(f"  cn_full assembled, RSS={_mem():.2f} GB", flush=True)
    bubble_id_per_kmer = np.concatenate(bubble_ids)
    bubble_chrom = np.concatenate(bubble_chrom_l)
    bubble_start = np.concatenate(bubble_start_l)
    bubble_end = np.concatenate(bubble_end_l)

    # If multiple chroms, we need to re-sort the concatenated u64 panel
    # because each chrom's sort_idx was sorted within-chrom.
    if len(sorted_panel_blocks) > 1:
        all_u64 = np.concatenate(sorted_panel_blocks)
        all_idx = np.concatenate(sort_idx_blocks)
        order = np.argsort(all_u64)
        sorted_panel_u64 = all_u64[order]
        sort_idx = all_idx[order].astype(np.int64)
        del all_u64, all_idx, order
    else:
        sorted_panel_u64 = sorted_panel_blocks[0]
        sort_idx = sort_idx_blocks[0].astype(np.int64)
    del sorted_panel_blocks, sort_idx_blocks

    n_founders = len(founders)
    K = cn_full.shape[1]
    print(f"  cn_full: {cn_full.shape} ({cn_full.nnz:,} nnz), "
          f"founders: {n_founders}, k-mers: {K:,}, "
          f"sorted_panel_u64: {sorted_panel_u64.nbytes/1e9:.2f} GB",
          flush=True)

    # 2. Define windows and assign each k-mer to its window
    print("defining windows + assigning k-mers...", flush=True)
    blocks = define_windows(bubble_chrom, bubble_start, bubble_end,
                            window_bp=args.window_bp)
    n_windows = len(blocks)
    print(f"  defined {n_windows:,} windows of {args.window_bp:,} bp", flush=True)
    kmer_window = assign_kmers_to_windows(
        bubble_id_per_kmer, bubble_chrom, bubble_start, bubble_end, blocks)
    n_assigned = (kmer_window >= 0).sum()
    print(f"  k-mers assigned to a window: {n_assigned:,}/{K:,} "
          f"({n_assigned/K:.1%})", flush=True)

    # 3. (panel u64 already loaded from cache above)

    # 4. Pack cn_full into bit-masked uint64 columns; free sparse cn_full
    print(f"bit-packing cn_full into ({(n_founders+63)//64}, {K}) uint64...",
          flush=True)
    t = time.time()
    packed_cn = pack_cn_full_as_bitmasks(cn_full, n_founders)
    print(f"  packed shape={packed_cn.shape}, {packed_cn.nbytes/1e9:.2f} GB "
          f"[{time.time()-t:.0f}s]", flush=True)
    del cn_full; gc.collect()

    # 5. Pseudoalign reads → per-(window, mask) EC counts (or load from cache)
    if args.load_ec_cache and os.path.exists(args.load_ec_cache):
        print(f"loading EC counts from {args.load_ec_cache}...", flush=True)
        t = time.time()
        cache = np.load(args.load_ec_cache, allow_pickle=False)
        ec_w = cache["ec_w"]
        ec_mask = cache["ec_mask"]   # (n_ec, n_words) uint64
        ec_count = cache["ec_count"]
        n_used = int(cache["n_used"])
        n_total = int(cache["n_total"])
        print(f"  loaded {len(ec_w):,} EC ({n_used:,}/{n_total:,} reads "
              f"informative) [{time.time()-t:.0f}s]", flush=True)
        ec_by_window: dict[int, list[tuple[bytes, int]]] = defaultdict(list)
        for i in range(len(ec_w)):
            ec_by_window[int(ec_w[i])].append(
                (ec_mask[i].tobytes(), int(ec_count[i])))
    else:
        print(f"pseudoaligning reads from {args.reads}...", flush=True)
        t = time.time()
        ec_counts, n_used, n_total = pseudoalign_reads_to_ec(
            args.reads, sorted_panel_u64, sort_idx, kmer_window, packed_cn,
            n_founders=n_founders, n_windows=n_windows,
            k=31, max_reads_per_file=args.max_reads_per_file)
        print(f"  done in {time.time()-t:.0f}s. {n_used:,}/{n_total:,} reads "
              f"informative. {len(ec_counts):,} unique EC.", flush=True)
        ec_by_window: dict[int, list[tuple[bytes, int]]] = defaultdict(list)
        for (w, mask_bytes), count in ec_counts.items():
            ec_by_window[w].append((mask_bytes, count))

        if args.save_ec_cache:
            print(f"saving EC cache → {args.save_ec_cache}...", flush=True)
            n_words = packed_cn.shape[0]
            ec_w_arr, ec_mask_arr, ec_count_arr = [], [], []
            for w, ecs in ec_by_window.items():
                for mb, cnt in ecs:
                    ec_w_arr.append(w)
                    ec_mask_arr.append(np.frombuffer(mb, dtype=np.uint64))
                    ec_count_arr.append(cnt)
            np.savez_compressed(
                args.save_ec_cache,
                ec_w=np.asarray(ec_w_arr, dtype=np.int32),
                ec_mask=np.stack(ec_mask_arr).astype(np.uint64),
                ec_count=np.asarray(ec_count_arr, dtype=np.int64),
                n_used=np.int64(n_used),
                n_total=np.int64(n_total),
                n_windows=np.int32(n_windows),
                n_founders=np.int32(n_founders))

    n_w_with_ec = sum(1 for w in ec_by_window if ec_by_window[w])
    print(f"  windows with EC: {n_w_with_ec}/{n_windows}", flush=True)

    # 6. Per-window EM (optionally with global-anchor KL prior)
    print(f"running per-window EM on EC counts "
          f"(global_anchor_weight={args.global_anchor_weight})...", flush=True)
    t = time.time()
    h_blocks, block_status, h_global = em_solve_per_window_ec(
        ec_by_window, n_windows, n_founders,
        max_iter=args.em_max_iter, tol=1e-6,
        global_anchor_weight=args.global_anchor_weight)
    print(f"  EM done in {time.time()-t:.0f}s. "
          f"{(block_status==0).sum()}/{n_windows} windows solved.", flush=True)

    # 7. Project per-window h to per-record AFs
    # Free panel arrays no longer needed (packed_cn, sorted_panel_u64, sort_idx)
    print(f"freeing panel arrays before projection (RSS pre-free: {_mem():.2f} GB)...",
          flush=True)
    del packed_cn, sorted_panel_u64, sort_idx, kmer_window
    gc.collect()
    print(f"  RSS post-free: {_mem():.2f} GB", flush=True)

    print("projecting to per-record AFs...", flush=True)
    cn_var = load_npz(args.cn_var)
    print(f"  cn_var loaded ({cn_var.shape}, {cn_var.nnz:,} nnz), "
          f"RSS={_mem():.2f} GB", flush=True)
    var_meta = np.load(args.cn_var_meta, allow_pickle=True)
    rec_chrom = np.asarray(var_meta["chrom"]).astype(str)
    rec_pos = np.asarray(var_meta["pos"]).astype(np.int64)

    # Vectorized record → window assignment (per chromosome via searchsorted)
    rec_block = np.full(len(rec_chrom), -1, dtype=np.int32)
    chrom_to_block_idxs: dict[str, list[int]] = {}
    for i, b in enumerate(blocks):
        chrom_to_block_idxs.setdefault(b.chrom, []).append(i)
    for c, idxs in chrom_to_block_idxs.items():
        idxs = np.array(sorted(idxs, key=lambda i: blocks[i].start),
                        dtype=np.int64)
        starts = np.array([blocks[i].start for i in idxs], dtype=np.int64)
        ends = np.array([blocks[i].end for i in idxs], dtype=np.int64)
        chr_mask = (rec_chrom == c)
        if not chr_mask.any():
            continue
        positions = rec_pos[chr_mask]
        # find rightmost block whose start <= position
        right = np.searchsorted(starts, positions, side="right") - 1
        valid = (right >= 0) & (right < len(starts))
        # require position <= ends[right]
        candidate = np.where(valid, right, 0)
        in_window = valid & (positions <= ends[candidate])
        assigned = np.where(in_window, idxs[candidate], -1)
        rec_block[chr_mask] = assigned.astype(np.int32)
    print(f"  records assigned to a window: "
          f"{(rec_block >= 0).sum():,}/{len(rec_block):,}",
          flush=True)
    fallback_h = np.full(n_founders, 1.0 / n_founders, dtype=np.float32)
    freqs = project_blocks_to_records_simple(h_blocks, cn_var, rec_block,
                                             fallback_h)

    # 8. Write TSV
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    print(f"writing TSV to {args.out}...", flush=True)
    ref_arr = np.asarray(var_meta["ref_len"]).astype(np.int64)
    alt_arr = np.asarray(var_meta["alt_len"]).astype(np.int64)
    with open(args.out, "w") as f:
        f.write("chrom\tpos\tref_len\talt_len\talt_freq\n")
        for i in range(len(rec_chrom)):
            f.write(f"{rec_chrom[i]}\t{rec_pos[i]}\t{ref_arr[i]}\t"
                    f"{alt_arr[i]}\t{freqs[i]:.6f}\n")

    if args.save_h_blocks:
        np.savez_compressed(args.save_h_blocks,
                            h_blocks=h_blocks,
                            h_global=h_global,
                            block_status=block_status,
                            block_chrom=np.array([b.chrom for b in blocks]),
                            block_start=np.array([b.start for b in blocks]),
                            block_end=np.array([b.end for b in blocks]),
                            founders=founders,
                            global_anchor_weight=np.float32(args.global_anchor_weight))
        print(f"saved h_blocks → {args.save_h_blocks}", flush=True)
    print(f"=== done. wall: {time.time()-t0:.0f}s ===", flush=True)


if __name__ == "__main__":
    main()
