"""Build per-block haplotype-level k-mer cn matrices for one chromosome.

We use the BigLD block index ONLY for block POSITIONS (chrom, start, end).
For each block we then dedup founders independently by their actual k-mer
content (full FASTA slice → set of canonical k-mers). Founders with
identical block sequences — same SNPs AND same SVs — collapse into one
class. Founders that share SNPs but differ on SVs become DIFFERENT classes,
so SV-distinguishing k-mers can drive their separate per-class freqs in
the EM.

This is strictly more informative than SNP-only dedup (which is what
hapFIRE's hap_idx_d1 carries): unique-hap classes here are sequence-unique,
so cn[class, kmer] is naturally 0/1, no rep-founder bias.

We further drop k-mers that appear discriminative in >1 block so per-block
EMs operate on disjoint k-mer sets — read counts are attributable to
exactly one block.

Output (single .npz per chromosome): per-block sequence-unique class
matrices, founder→class mapping per block, and the kept k-mer strings.
"""
from __future__ import annotations
import argparse
import gc
import os
import sys
import time
from collections import defaultdict
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pysam
from scipy.sparse import csr_matrix, csc_matrix, save_npz


_RC = str.maketrans("ACGT", "TGCA")

# 2-bit base encoding lookup. Byte values are A=65,C=67,G=71,T=83 (ASCII).
# Build a 256-byte LUT mapping ASCII byte → 2-bit code (0/1/2/3) or 255 for non-ACGT.
_BASE_LUT = np.full(256, 255, dtype=np.uint8)
for c, code in zip(b"ACGT", [0, 1, 2, 3]):
    _BASE_LUT[c] = code
for c, code in zip(b"acgt", [0, 1, 2, 3]):
    _BASE_LUT[c] = code


def canonical_kmers(seq: str, k: int) -> set[str]:
    """Pure-Python fallback (slow for large blocks). Use canonical_kmer_hashes
    + decode for production.
    """
    seq = seq.upper()
    out = set()
    for i in range(len(seq) - k + 1):
        kmer = seq[i:i + k]
        if any(c not in "ACGT" for c in kmer):
            continue
        rc = kmer.translate(_RC)[::-1]
        out.add(min(kmer, rc))
    return out


def canonical_kmer_hashes(seq: str, k: int) -> np.ndarray:
    """Numpy-accelerated canonical k-mer hashing.

    Returns a 1D uint64 array of UNIQUE canonical k-mer hashes from seq
    (k <= 32, so the 2-bit packing fits in uint64). N-containing k-mers
    are skipped.

    For k=31 each hash uses 62 bits (top 2 bits zero). Canonical = min(
    fwd_hash, rc_hash). The hash is bijective with the canonical k-mer
    sequence (modulo bit-twiddling), so np.unique on the hash array is
    equivalent to set-of-canonical-kmer-strings.
    """
    if k > 32:
        raise ValueError(f"k must be <= 32 for uint64 packing; got {k}")
    if len(seq) < k:
        return np.empty(0, dtype=np.uint64)
    # Encode: ASCII bytes → 2-bit codes, with 255 for non-ACGT
    s = np.frombuffer(seq.encode("ascii"), dtype=np.uint8)
    codes = _BASE_LUT[s]  # (L,) uint8 in {0,1,2,3,255}

    L = codes.size
    n_kmers = L - k + 1
    if n_kmers <= 0:
        return np.empty(0, dtype=np.uint64)

    # Rolling-window (n_kmers, k) uint8
    win = np.lib.stride_tricks.sliding_window_view(codes, k)  # (n_kmers, k)

    # Mark k-mers containing any 255 (non-ACGT base)
    valid = (win != 255).all(axis=1)

    # Pack 2-bit codes into uint64 forward hash. shifts = [2*(k-1), 2*(k-2), ..., 0]
    shifts = (np.arange(k - 1, -1, -1, dtype=np.int64) * 2).astype(np.uint64)
    win64 = win.astype(np.uint64)
    fwd = (win64 << shifts[None, :]).sum(axis=1)  # (n_kmers,) uint64

    # Reverse complement: complement each base (XOR 3), then reverse position-order.
    # rc_hash = sum over j of (3^code[k-1-j]) << (2*(k-1-j))
    #         = sum over i of (3^code[i]) << (2*i)
    rc_codes = (win ^ 3).astype(np.uint64)  # (n_kmers, k) — complement
    # Reversed shift order: the i-th position from the end gets shift 2*i
    rc_shifts = (np.arange(k, dtype=np.int64) * 2).astype(np.uint64)
    # In rc, original position k-1-j contributes to bit position j; so for
    # window position i (0..k-1), rc bit position is k-1-i. Equivalent:
    # rc = sum over i of complement(win[i]) << (2 * (k-1-i))
    # which equals fwd computed on reversed-and-complemented sequence.
    rc_win = rc_codes[:, ::-1]
    rc = (rc_win << shifts[None, :]).sum(axis=1)  # same shifts as fwd but applied to rev-complemented
    # Actually we want: rc_hash for each kmer is the fwd_hash of (reverse(complement(kmer)))
    # rev(complement(kmer)) per row = rc_win above (already reversed and complemented)
    # So rc = (rc_win << shifts).sum() — same packing as fwd. Already done above.

    canonical = np.minimum(fwd, rc)
    canonical = canonical[valid]
    return np.unique(canonical)


def hash_to_kmer(h: int, k: int) -> str:
    """Decode a 2-bit canonical k-mer hash back to its canonical string."""
    bases = "ACGT"
    out = []
    for i in range(k - 1, -1, -1):
        out.append(bases[(h >> (2 * i)) & 3])
    return "".join(out)


def canonical_kmer_hashes_with_counts(seq: str, k: int):
    """Like canonical_kmer_hashes but returns (sorted_unique_hashes, counts).

    counts[i] = number of times sorted_unique_hashes[i] appears in seq.
    Used for genomic-uniqueness filtering: any k-mer with count > 1 in
    a single founder's full Chr1 is "multi-mapping" within that founder.
    """
    if k > 32:
        raise ValueError(f"k must be <= 32 for uint64 packing; got {k}")
    if len(seq) < k:
        return (np.empty(0, dtype=np.uint64), np.empty(0, dtype=np.int32))
    s = np.frombuffer(seq.encode("ascii"), dtype=np.uint8)
    codes = _BASE_LUT[s]
    L = codes.size
    n_kmers = L - k + 1
    if n_kmers <= 0:
        return (np.empty(0, dtype=np.uint64), np.empty(0, dtype=np.int32))
    win = np.lib.stride_tricks.sliding_window_view(codes, k)
    valid = (win != 255).all(axis=1)
    shifts = (np.arange(k - 1, -1, -1, dtype=np.int64) * 2).astype(np.uint64)
    win64 = win.astype(np.uint64)
    fwd = (win64 << shifts[None, :]).sum(axis=1)
    rc_codes = (win ^ 3).astype(np.uint64)
    rc_win = rc_codes[:, ::-1]
    rc = (rc_win << shifts[None, :]).sum(axis=1)
    canonical = np.minimum(fwd, rc)
    canonical = canonical[valid]
    return np.unique(canonical, return_counts=True)


# Module-level cache for FASTAs in worker processes
_FA_CACHE: dict[str, pysam.FastaFile] = {}


def _open_fa(path: str) -> pysam.FastaFile:
    fa = _FA_CACHE.get(path)
    if fa is None:
        fa = pysam.FastaFile(path)
        _FA_CACHE[path] = fa
    return fa


def collapse_classes_by_jaccard(class_hashes: list[np.ndarray],
                                founder_to_class: np.ndarray,
                                max_jaccard_dist: float):
    """PHG-style allele grouping: merge classes within max_jaccard_dist
    (1 - Jaccard) of each other via connected components.

    Vectorized implementation: build sparse class×kmer carrier matrix once,
    compute pairwise intersections via cn @ cn.T (single sparse matmul),
    and apply the threshold. ~1000× faster than a Python pairwise loop for
    n>50.

    class_hashes: list of (k_i,) uint64 arrays (sorted, unique k-mer hashes
                  per class).
    founder_to_class: (N_eco,) int32 mapping founders to original class idx.
    max_jaccard_dist: 0..1, maximum 1-Jaccard distance to merge. mxDiv in PHG.
                      0 = no merging (strict identity); 1 = merge everything.

    Returns: (new_class_hashes, new_founder_to_class) — merged.
    """
    n = len(class_hashes)
    if n <= 1 or max_jaccard_dist <= 0:
        return class_hashes, founder_to_class

    # Build the union-of-hashes index (sorted) for vectorized class membership
    union = np.unique(np.concatenate(class_hashes)) if class_hashes else \
        np.empty(0, dtype=np.uint64)
    if union.size == 0:
        return class_hashes, founder_to_class

    # Build sparse cn[class, kmer-in-union]
    rows: list[np.ndarray] = []
    cols: list[np.ndarray] = []
    for c, h in enumerate(class_hashes):
        if h.size == 0:
            continue
        idx = np.searchsorted(union, h)
        cols.append(idx)
        rows.append(np.full(idx.size, c, dtype=np.int32))
    if not rows:
        return class_hashes, founder_to_class
    rows_arr = np.concatenate(rows)
    cols_arr = np.concatenate(cols)
    data = np.ones(rows_arr.size, dtype=np.int32)
    cn = csr_matrix((data, (rows_arr, cols_arr)),
                    shape=(n, union.size), dtype=np.int32)

    # Pairwise intersection: cn @ cn.T → (n, n) int32
    inter = (cn @ cn.T).toarray()  # |class_i ∩ class_j|
    sizes = np.asarray(cn.sum(axis=1)).flatten()  # |class_i|
    # |class_i ∪ class_j| = sizes[i] + sizes[j] - inter[i, j]
    union_pw = sizes[:, None] + sizes[None, :] - inter
    # Avoid div-by-zero; classes with size 0 stay singletons
    safe_union = np.where(union_pw > 0, union_pw, 1)
    jaccard = inter / safe_union
    # Merge edges: (1 - jaccard) <= max_jaccard_dist, i.e., jaccard >= 1 - threshold
    merge_thresh = 1.0 - max_jaccard_dist
    # Build edges (upper triangular only)
    iu, ju = np.where(np.triu(jaccard >= merge_thresh, k=1))

    # Union-find over the n classes
    parent = list(range(n))
    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    for a, b in zip(iu, ju):
        ra, rb = find(int(a)), find(int(b))
        if ra != rb:
            parent[ra] = rb

    # Build merged classes by component
    comp_to_members: dict[int, list[int]] = defaultdict(list)
    for i in range(n):
        comp_to_members[find(i)].append(i)
    comps_sorted = sorted(comp_to_members.values(), key=lambda m: min(m))
    old_to_new = np.empty(n, dtype=np.int32)
    new_class_hashes: list[np.ndarray] = []
    for new_id, members in enumerate(comps_sorted):
        for m in members:
            old_to_new[m] = new_id
        merged = np.unique(np.concatenate([class_hashes[m] for m in members])) \
            if members else np.empty(0, dtype=np.uint64)
        new_class_hashes.append(merged)
    new_founder_to_class = old_to_new[founder_to_class].astype(np.int32)
    return new_class_hashes, new_founder_to_class


def process_block(args_tuple):
    """Build per-block sequence-unique class × kmer cn matrix using
    numpy-vectorized canonical k-mer hashing.

    For each of the 231 founders, slice its FASTA in the block window,
    compute the sorted unique uint64 canonical k-mer hashes via
    canonical_kmer_hashes(). Then dedup founders by their k-mer-hash array:
    founders with identical block sequences (same SNPs AND same SVs)
    collapse into one class.

    If jaccard_mxdiv > 0, additionally collapse near-identical classes via
    pairwise k-mer Jaccard (PHG-style). Two classes are merged when
    1 - Jaccard(set_i, set_j) <= jaccard_mxdiv.

    Returns: (block_idx, kmer_hashes_kept (uint64 array),
              cn_kept_uint8 (n_uniq, K_kept), founder_to_class int32 (N_eco,))

    Discriminative filter applied: k-mers in [1, n_uniq-1] carriers within
    block. Cross-block dedup is applied later, globally.
    """
    (block_idx, chrom, pos_start, pos_end,
     ecotypes, fastas_dir, k, jaccard_mxdiv, mask_positions_for_chrom) = args_tuple
    n_eco = len(ecotypes)

    # Build a per-base relative mask for this block window.
    # mask_positions_for_chrom is a sorted np.int64 array of 1-based positions
    # to N-out for the current chromosome (low-MAF SNP positions). We slice
    # to this block's window and convert to 0-based offsets within the slice.
    if mask_positions_for_chrom is not None and mask_positions_for_chrom.size > 0:
        lo = np.searchsorted(mask_positions_for_chrom, pos_start, side="left")
        hi = np.searchsorted(mask_positions_for_chrom, pos_end, side="right")
        mask_in_block = mask_positions_for_chrom[lo:hi]
        mask_offsets = mask_in_block - pos_start  # 0-based offsets in the slice
    else:
        mask_offsets = None

    # Hash each founder's full FASTA slice → sorted uint64 hash array
    founder_hashes: list[np.ndarray] = []
    for f in range(n_eco):
        eco = ecotypes[f]
        fa = _open_fa(os.path.join(fastas_dir, f"{eco}.chr.fa"))
        # 1-based pos_start, pysam fetch is 0-based half-open
        seq = fa.fetch(chrom, pos_start - 1, pos_end)
        if mask_offsets is not None and mask_offsets.size > 0:
            # N-out low-MAF SNP positions; canonical_kmer_hashes will skip
            # k-mers overlapping any N.
            seq_arr = np.frombuffer(seq.encode("ascii"), dtype=np.uint8).copy()
            valid = (mask_offsets >= 0) & (mask_offsets < seq_arr.size)
            seq_arr[mask_offsets[valid]] = ord("N")
            seq = seq_arr.tobytes().decode("ascii")
        founder_hashes.append(canonical_kmer_hashes(seq, k=k))

    # Dedup founders by hash array: convert each array to bytes for hashing
    # in a dict (numpy arrays aren't hashable directly; tobytes() is).
    seen: dict[bytes, int] = {}
    founder_to_class = np.empty(n_eco, dtype=np.int32)
    class_hashes: list[np.ndarray] = []
    for f, h in enumerate(founder_hashes):
        key = h.tobytes()
        cid = seen.get(key)
        if cid is None:
            cid = len(class_hashes)
            seen[key] = cid
            class_hashes.append(h)
        founder_to_class[f] = cid

    # Optional PHG-style allele grouping by k-mer Jaccard distance
    if jaccard_mxdiv > 0 and len(class_hashes) > 1:
        class_hashes, founder_to_class = collapse_classes_by_jaccard(
            class_hashes, founder_to_class, jaccard_mxdiv)

    n_uniq = len(class_hashes)

    # Build union of k-mer hashes across all classes
    if not class_hashes:
        return (block_idx, np.empty(0, dtype=np.uint64),
                np.zeros((n_uniq, 0), dtype=np.uint8), founder_to_class)
    union_hashes = np.unique(np.concatenate(class_hashes)) if class_hashes else \
        np.empty(0, dtype=np.uint64)
    if union_hashes.size == 0:
        return (block_idx, union_hashes,
                np.zeros((n_uniq, 0), dtype=np.uint8), founder_to_class)

    # Build cn[class, kmer]: for each class, mark which positions in
    # union_hashes it carries (np.searchsorted is fast vs the sorted union).
    # We accumulate (row, col) coords for sparse construction — density ~0.01
    # so dense (n_uniq × n_kmer) is ~50× larger than sparse equivalent.
    rows: list[np.ndarray] = []
    cols: list[np.ndarray] = []
    for c, h in enumerate(class_hashes):
        if h.size == 0:
            continue
        idx = np.searchsorted(union_hashes, h)
        cols.append(idx)
        rows.append(np.full(idx.size, c, dtype=np.int32))
    if not rows:
        return (block_idx, np.empty(0, dtype=np.uint64),
                csr_matrix((n_uniq, 0), dtype=np.uint8), founder_to_class)
    rows_arr = np.concatenate(rows)
    cols_arr = np.concatenate(cols)
    data = np.ones(rows_arr.size, dtype=np.uint8)
    cn = csr_matrix((data, (rows_arr, cols_arr)),
                    shape=(n_uniq, union_hashes.size), dtype=np.uint8)

    # Discriminative filter: keep cols where 1 <= sum < n_uniq
    col_sum = np.asarray(cn.sum(axis=0)).flatten()
    discrim = (col_sum > 0) & (col_sum < n_uniq)
    if not discrim.any():
        return (block_idx, np.empty(0, dtype=np.uint64),
                csr_matrix((n_uniq, 0), dtype=np.uint8), founder_to_class)
    keep_idx = np.where(discrim)[0]
    kmer_hashes_kept = union_hashes[keep_idx]
    cn_kept = cn[:, keep_idx].tocsr()
    return (block_idx, kmer_hashes_kept, cn_kept, founder_to_class)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--block-index", required=True,
                    help="hapfire_block_index_<chrom>.npz")
    ap.add_argument("--fastas-dir", required=True,
                    help="dir of <eco>.chr.fa per ecotype")
    ap.add_argument("--out", required=True,
                    help="output .npz path")
    ap.add_argument("--k", type=int, default=31)
    ap.add_argument("--chrom-prefix", default="Chr",
                    help="prefix to add to block_chrom values from index")
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--chrom-filter", default=None,
                    help="if given (e.g. 'Chr1'), only build blocks on this chrom")
    ap.add_argument("--max-blocks", type=int, default=None,
                    help="for testing: only process the first N blocks")
    ap.add_argument("--jaccard-mxdiv", type=float, default=0.0,
                    help="PHG-style allele grouping: merge classes whose "
                         "1-Jaccard k-mer distance is <= this value. "
                         "0.0 (default) = strict identity (no merging); "
                         "0.05 = merge classes within ~5%% k-mer difference.")
    ap.add_argument("--snp-maf-threshold", type=float, default=0.0,
                    help="Mask SNP positions in panel cn_var with MAF below "
                         "this threshold (replace with 'N' in founder FASTA "
                         "slice; k-mers overlapping N are dropped). 0.0 "
                         "(default) = no masking; 0.05 = drop SNPs at <5%% MAF.")
    ap.add_argument("--cn-var", default=None,
                    help="path to cn_var .npz for MAF computation; required "
                         "if --snp-maf-threshold > 0")
    ap.add_argument("--cn-var-meta", default=None,
                    help="path to cn_var .meta.npz; required if "
                         "--snp-maf-threshold > 0")
    args = ap.parse_args()

    z = np.load(args.block_index, allow_pickle=True)
    ecotypes = np.asarray(z["ecotypes"]).astype(str)
    block_chrom = np.asarray(z["block_chrom"]).astype(str)
    block_pos_start = np.asarray(z["block_pos_start"]).astype(np.int64)
    block_pos_end = np.asarray(z["block_pos_end"]).astype(np.int64)
    n_eco = len(ecotypes)
    # We only use the (chrom, pos_start, pos_end) positions from the block
    # index — NOT the SNP-derived hap_idx_d1. Per-block dedup is recomputed
    # from the full FASTA sequence (captures SVs as well as SNPs).

    # Optional MAF mask: from cn_var, identify SNP records (ref_len==alt_len==1)
    # with MAF below threshold and collect their positions per chromosome.
    # FASTA slicing per founder will N-out these positions, dropping k-mers
    # that overlap them and effectively merging founders that differed only
    # at low-MAF positions.
    mask_positions_per_chrom: dict[str, np.ndarray] = {}
    if args.snp_maf_threshold > 0:
        if not args.cn_var or not args.cn_var_meta:
            raise SystemExit("--snp-maf-threshold > 0 requires --cn-var and --cn-var-meta")
        from scipy.sparse import load_npz
        print(f"Loading cn_var for MAF mask: {args.cn_var}", flush=True)
        cn_var = load_npz(args.cn_var)
        cv_meta = np.load(args.cn_var_meta, allow_pickle=True)
        rec_chrom = np.asarray(cv_meta["chrom"]).astype(str)
        rec_pos = np.asarray(cv_meta["pos"]).astype(np.int64)
        rec_ref = np.asarray(cv_meta["ref_len"]).astype(np.int64)
        rec_alt = np.asarray(cv_meta["alt_len"]).astype(np.int64)
        ac = np.asarray(cn_var.sum(axis=0)).flatten()
        n_founders = cn_var.shape[0]
        af = ac / n_founders
        maf = np.minimum(af, 1.0 - af)
        is_snp = (rec_ref == 1) & (rec_alt == 1)
        to_mask = is_snp & (maf < args.snp_maf_threshold)
        n_to_mask = int(to_mask.sum())
        print(f"  SNP records below MAF {args.snp_maf_threshold}: {n_to_mask:,} / {is_snp.sum():,}", flush=True)
        for c in np.unique(rec_chrom):
            sel = to_mask & (rec_chrom == c)
            if sel.any():
                mask_positions_per_chrom[c] = np.sort(rec_pos[sel])
            else:
                mask_positions_per_chrom[c] = np.empty(0, dtype=np.int64)

    # Build chrom labels (e.g. '1' -> 'Chr1')
    block_chrom_full = np.array(
        [f"{args.chrom_prefix}{c}" if not c.startswith(args.chrom_prefix) else c
         for c in block_chrom], dtype=object)

    # Filter
    n_blocks_total = len(block_chrom_full)
    sel = np.ones(n_blocks_total, dtype=bool)
    if args.chrom_filter:
        sel &= (block_chrom_full == args.chrom_filter)
    if args.max_blocks:
        idx_sel = np.where(sel)[0][:args.max_blocks]
        sel[:] = False
        sel[idx_sel] = True
    keep_idx = np.where(sel)[0]
    print(f"Processing {len(keep_idx)} / {n_blocks_total} blocks "
          f"(chrom_filter={args.chrom_filter})", flush=True)

    # Build job list
    jobs = []
    for b in keep_idx:
        chrom_b = str(block_chrom_full[b])
        # mask_positions_per_chrom is keyed by the chrom labels in cn_var meta
        # (could be "1" or "Chr1"); accept either.
        mp = mask_positions_per_chrom.get(chrom_b)
        if mp is None:
            mp = mask_positions_per_chrom.get(chrom_b.replace("Chr", ""))
        jobs.append((
            int(b),
            chrom_b,
            int(block_pos_start[b]),
            int(block_pos_end[b]),
            ecotypes,
            args.fastas_dir,
            args.k,
            float(args.jaccard_mxdiv),
            mp,
        ))

    # Worker pool
    t0 = time.time()
    results: dict[int, tuple[list[str], np.ndarray, np.ndarray]] = {}
    if args.threads > 1:
        with Pool(args.threads) as pool:
            for i, res in enumerate(pool.imap_unordered(process_block, jobs, chunksize=4)):
                bidx, kmers, cn, f2c = res
                results[bidx] = (kmers, cn, f2c)
                if (i + 1) % 200 == 0:
                    print(f"  {i+1}/{len(jobs)} blocks "
                          f"[{time.time()-t0:.0f}s]", flush=True)
    else:
        for i, j in enumerate(jobs):
            bidx, kmers, cn, f2c = process_block(j)
            results[bidx] = (kmers, cn, f2c)
            if (i + 1) % 200 == 0:
                print(f"  {i+1}/{len(jobs)} blocks "
                      f"[{time.time()-t0:.0f}s]", flush=True)
    print(f"per-block hashing done [{time.time()-t0:.0f}s]", flush=True)

    # Cross-block global dedup: drop k-mer hashes present in >1 block's discriminative set
    print("global cross-block dedup...", flush=True)
    t = time.time()
    # Concatenate all blocks' kept hashes with a "block_id" tag, then count
    # unique-hash → number of blocks. Use np.unique with return_counts on
    # the concatenated hash array tagged so duplicates within a single block
    # don't count (they shouldn't exist anyway since each block returns
    # np.unique).
    all_hashes = []
    for b in keep_idx:
        h = results[int(b)][0]
        if h.size:
            all_hashes.append(h)
    if all_hashes:
        concat = np.concatenate(all_hashes)
        uniq, counts = np.unique(concat, return_counts=True)
        # Build set of "kmers in only one block" → globally-unique
        unique_set = uniq[counts == 1]
        unique_set_view = unique_set  # already sorted by np.unique
    else:
        unique_set_view = np.empty(0, dtype=np.uint64)
        uniq = unique_set_view
    n_global = uniq.size if all_hashes else 0
    n_unique_global = int(unique_set_view.size)
    n_dropped = n_global - n_unique_global
    print(f"  global k-mer hashes: {n_global:,}; "
          f"unique to one block: {n_unique_global:,}; "
          f"dropped (in >1 block): {n_dropped:,} "
          f"({n_dropped/max(1,n_global)*100:.1f}%) "
          f"[{time.time()-t:.0f}s]", flush=True)

    # Build the unified output
    print("building output structures...", flush=True)
    t = time.time()
    block_idx_arr = []
    block_kmer_offsets = [0]
    kmer_hashes_all: list[np.ndarray] = []
    block_n_haps = []
    block_n_kmers = []
    block_chrom_kept = []
    block_pos_start_kept = []
    block_pos_end_kept = []
    founder_to_class_kept = []  # per block: (n_eco,) int32
    cn_local_blocks: list[csr_matrix] = []  # per-block sparse CSR (n_uniq, n_kmer_kept)

    for b in keep_idx:
        kmer_hashes, cn_full, f2c = results[int(b)]
        n_uniq_b = cn_full.shape[0]
        # cn_full is now a sparse CSR matrix from the worker.
        # Filter to globally-unique hashes (np.searchsorted on sorted unique_set_view)
        if kmer_hashes.size == 0:
            keep_local_idx = np.empty(0, dtype=np.int64)
        else:
            pos = np.searchsorted(unique_set_view, kmer_hashes)
            in_set = (pos < unique_set_view.size) & \
                     (unique_set_view[np.clip(pos, 0, unique_set_view.size - 1)] == kmer_hashes)
            keep_local_idx = np.where(in_set)[0]
        kept_hashes = kmer_hashes[keep_local_idx] if keep_local_idx.size else \
            np.empty(0, dtype=np.uint64)
        if cn_full.shape[1] > 0 and keep_local_idx.size > 0:
            cn_kept = cn_full[:, keep_local_idx].tocsr()
        else:
            cn_kept = csr_matrix((n_uniq_b, 0), dtype=np.uint8)

        block_idx_arr.append(int(b))
        kmer_hashes_all.append(kept_hashes)
        block_kmer_offsets.append(block_kmer_offsets[-1] + kept_hashes.size)
        block_n_haps.append(int(n_uniq_b))
        block_n_kmers.append(int(kept_hashes.size))
        block_chrom_kept.append(str(block_chrom_full[b]))
        block_pos_start_kept.append(int(block_pos_start[b]))
        block_pos_end_kept.append(int(block_pos_end[b]))
        founder_to_class_kept.append(f2c)
        cn_local_blocks.append(cn_kept)

    kmer_hashes_concat = (np.concatenate(kmer_hashes_all)
                          if kmer_hashes_all else np.empty(0, dtype=np.uint64))
    # Decode hashes to canonical k-mer strings (needed for jellyfish query later)
    print(f"  decoding {kmer_hashes_concat.size:,} k-mer hashes to strings "
          f"(k={args.k})...", flush=True)
    kmer_strings_all = [hash_to_kmer(int(h), args.k) for h in kmer_hashes_concat]

    # Concatenate per-block sparse cn matrices into one big sparse CSR.
    # Rows: per-block hap-rows stacked block by block; cols: per-block kmers
    # mapped to consecutive ranges in kmer_strings_all per block_kmer_offsets.
    print(f"  concatenating cn rows for {len(cn_local_blocks)} blocks...", flush=True)
    rows_per_block = block_n_haps
    total_rows = sum(rows_per_block)
    total_cols = block_kmer_offsets[-1]

    # Build big CSR by collecting (row_global, col_global) coords from each
    # block's existing sparse cn. Avoids reconstructing dense.
    indptr_list = [0]
    indices_chunks: list[np.ndarray] = []
    nnz_running = 0
    for bi, cn_b in enumerate(cn_local_blocks):
        n_haps_b = cn_b.shape[0]
        col_offset = block_kmer_offsets[bi]
        # cn_b is CSR — its indptr gives row boundaries; indices are cols.
        cb_indptr = cn_b.indptr
        cb_indices = cn_b.indices
        # Append indices with col_offset
        if cb_indices.size:
            indices_chunks.append(cb_indices + col_offset)
        # Append row indptr boundaries (relative offsets, then add nnz_running)
        for h in range(1, n_haps_b + 1):
            indptr_list.append(int(cb_indptr[h]) + nnz_running)
        nnz_running += int(cb_indptr[-1])
    cn_indptr = np.asarray(indptr_list, dtype=np.int64)
    cn_indices = (np.concatenate(indices_chunks).astype(np.int64)
                  if indices_chunks else np.empty(0, dtype=np.int64))
    cn_data = np.ones(cn_indices.size, dtype=np.uint8)

    # Per-block hap-row offset (cumulative, into the big CSR)
    block_hap_offsets = np.cumsum([0] + rows_per_block).astype(np.int64)

    print(f"  shape: rows={total_rows}, cols={total_cols}, nnz={cn_indices.size}; "
          f"density={cn_indices.size/(total_rows*max(1,total_cols)):.4f} "
          f"[{time.time()-t:.0f}s]", flush=True)

    # Save
    print(f"saving to {args.out}...", flush=True)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        args.out,
        ecotypes=ecotypes,
        block_idx=np.asarray(block_idx_arr, dtype=np.int64),
        block_chrom=np.asarray(block_chrom_kept, dtype=object),
        block_pos_start=np.asarray(block_pos_start_kept, dtype=np.int64),
        block_pos_end=np.asarray(block_pos_end_kept, dtype=np.int64),
        block_n_haps=np.asarray(block_n_haps, dtype=np.int32),
        block_n_kmers=np.asarray(block_n_kmers, dtype=np.int32),
        block_kmer_offsets=np.asarray(block_kmer_offsets, dtype=np.int64),
        block_hap_offsets=block_hap_offsets,
        kmer_strings=np.asarray(kmer_strings_all, dtype=object),
        cn_indptr=cn_indptr,
        cn_indices=cn_indices,
        cn_data=cn_data,
        # per-block: which sequence-unique class each founder maps to
        # shape (n_blocks, n_eco) int32. Replaces SNP-derived hap_idx_d1.
        founder_to_class=np.asarray(founder_to_class_kept, dtype=np.int32),
        k=args.k,
    )
    print(f"done. total wall: {time.time()-t0:.0f}s", flush=True)
    print(f"  total k-mers: {len(kmer_strings_all):,}")
    print(f"  total haps:   {total_rows:,}")
    print(f"  median kmers/block: {int(np.median(block_n_kmers))}")
    print(f"  median haps/block:  {int(np.median(block_n_haps))}")


if __name__ == "__main__":
    main()
