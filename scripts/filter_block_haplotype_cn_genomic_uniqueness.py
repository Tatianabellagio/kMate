"""PHG/PanGenie-style genomic-uniqueness filter for block_haplotype_cn npz.

Post-process step on a chr1_full*.npz output of build_block_haplotype_cn.py.

For each k-mer in the kept set, computes its **global occurrence count**
across all 231 founders' full Chr1 sequences. A k-mer is dropped if:
  - global_occurrence > expected_carriers (the discriminative carrier count
    within the k-mer's assigned block)

In other words: kept k-mers must appear at exactly one position per
*carrier* founder's Chr1, with NO additional off-target positions in the
genome. This mirrors PanGenie's "genomic_count == local_count" filter
applied at the block level instead of bubble level.

The cn matrix and kmer_strings are filtered to keep only good k-mers;
founder_to_class is unchanged (founder partitioning by block sequence
remains identical — we just drop noisy k-mers from each class's signature).

Usage:
    python filter_block_haplotype_cn_genomic_uniqueness.py \\
        --in chr1_full.npz --out chr1_full_genuniq.npz \\
        --fastas-dir founder_fastas_231 --threads 8
"""
from __future__ import annotations
import argparse
import os
import sys
import time
from pathlib import Path

import numpy as np
import pysam
from scipy.sparse import csr_matrix

# Reuse the hashing helpers from the builder
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from build_block_haplotype_cn import canonical_kmer_hashes_with_counts


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="in_npz", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--fastas-dir", required=True)
    ap.add_argument("--chrom", default="Chr1")
    ap.add_argument("--threads", type=int, default=1,
                    help="not used yet — sequential per-founder hashing")
    args = ap.parse_args()

    print(f"Loading {args.in_npz}...", flush=True)
    z = np.load(args.in_npz, allow_pickle=True)
    ecotypes = np.asarray(z["ecotypes"]).astype(str)
    block_idx = np.asarray(z["block_idx"]).astype(np.int64)
    block_chrom = np.asarray(z["block_chrom"]).astype(str)
    block_pos_start = np.asarray(z["block_pos_start"]).astype(np.int64)
    block_pos_end = np.asarray(z["block_pos_end"]).astype(np.int64)
    block_n_haps = np.asarray(z["block_n_haps"]).astype(np.int32)
    block_n_kmers = np.asarray(z["block_n_kmers"]).astype(np.int32)
    block_kmer_offsets = np.asarray(z["block_kmer_offsets"]).astype(np.int64)
    block_hap_offsets = np.asarray(z["block_hap_offsets"]).astype(np.int64)
    kmer_strings = z["kmer_strings"]  # keep as object array for now
    cn_indptr = np.asarray(z["cn_indptr"]).astype(np.int64)
    cn_indices = np.asarray(z["cn_indices"]).astype(np.int64)
    cn_data = np.asarray(z["cn_data"]).astype(np.uint8)
    founder_to_class = np.asarray(z["founder_to_class"]).astype(np.int32)
    k = int(z["k"])

    n_eco = len(ecotypes)
    n_blocks = len(block_idx)
    n_kmers_total = int(block_kmer_offsets[-1])
    total_haps = int(block_hap_offsets[-1])
    print(f"  blocks={n_blocks}, n_eco={n_eco}, total_kmers={n_kmers_total:,}, "
          f"total_haps={total_haps:,}", flush=True)

    # Step A: Vectorized kmer_strings → uint64 hashes.
    # Concatenate all k-mer strings into one big buffer, decode bases to 2-bit
    # codes, reshape, then bit-pack into uint64.
    print("Step A: encoding kmer strings → uint64 hashes (vectorized)...", flush=True)
    t = time.time()
    base_lut = np.full(256, 255, dtype=np.uint8)
    for c, code in zip(b"ACGT", [0, 1, 2, 3]):
        base_lut[c] = code
    # Concat all strings (each is k chars). Total bytes = n × k ≈ 1.1 GB for k=31, n=36M.
    big = "".join(kmer_strings)
    del kmer_strings  # free the object array; we'll reload its underlying ascii from `big`
    kmer_string_bytes = np.frombuffer(big.encode("ascii"), dtype=np.uint8)
    kmer_string_bytes = kmer_string_bytes.reshape(n_kmers_total, k)
    codes = base_lut[kmer_string_bytes].astype(np.uint64)
    shifts = (np.arange(k - 1, -1, -1, dtype=np.int64) * 2).astype(np.uint64)
    kept_hashes = (codes << shifts[None, :]).sum(axis=1)
    del codes, kmer_string_bytes
    print(f"  encoded {n_kmers_total:,} hashes [{time.time()-t:.0f}s]", flush=True)

    # IMPORTANT: kept_hashes is in BLOCK ORDER, not sorted! We need a sorted
    # lookup index. Build a sort-permutation and its inverse.
    sort_idx = np.argsort(kept_hashes, kind="stable")
    sorted_hashes = kept_hashes[sort_idx]
    inv_sort = np.empty_like(sort_idx)
    inv_sort[sort_idx] = np.arange(n_kmers_total, dtype=np.int64)

    # Step B: For each founder, hash full Chr1 with counts, accumulate
    # per-kmer global_count for our kept set
    print(f"Step B: per-founder hashing of full {args.chrom} "
          f"({n_eco} founders)...", flush=True)
    t0 = time.time()
    global_count = np.zeros(n_kmers_total, dtype=np.int32)
    n_founder_with_dups = 0
    n_dup_kmers = 0
    for f, eco in enumerate(ecotypes):
        fa_path = os.path.join(args.fastas_dir, f"{eco}.chr.fa")
        fa = pysam.FastaFile(fa_path)
        seq = fa.fetch(args.chrom)
        fa.close()
        h, cnt = canonical_kmer_hashes_with_counts(seq, k=k)
        # Find which of these hashes are in our kept set (sorted_hashes)
        pos = np.searchsorted(sorted_hashes, h)
        in_set = (pos < sorted_hashes.size) & \
                 (sorted_hashes[np.clip(pos, 0, sorted_hashes.size - 1)] == h)
        keep = np.where(in_set)[0]
        # Accumulate counts. Convert sorted-position to original block-order.
        sorted_pos = pos[keep]
        orig_idx = sort_idx[sorted_pos]
        global_count[orig_idx] += cnt[keep]
        # Track multi-occurrence within founder
        within = (cnt[keep] > 1).sum()
        if within > 0:
            n_founder_with_dups += 1
            n_dup_kmers += int(within)
        if (f + 1) % 30 == 0 or f + 1 == n_eco:
            print(f"  {f+1}/{n_eco} founders [{time.time()-t0:.0f}s]", flush=True)
    print(f"  done [{time.time()-t0:.0f}s]", flush=True)
    print(f"  founders with at least one within-founder duplicate kept-kmer: "
          f"{n_founder_with_dups}/{n_eco}")

    # Step C: Compute carrier_count per kept k-mer (within its block).
    # Each kmer column j belongs to exactly one block (by cross-block dedup),
    # so carrier_count[j] = sum of cn entries in column j (= #class-rows in
    # the assigned block that carry the kmer). Compute directly from the
    # CSR triplets: for each (row, col) entry, increment carrier_count[col].
    print("Step C: computing per-kmer carrier_count...", flush=True)
    t = time.time()
    # cn_indices already lists the (sparse) column for each nnz; counting
    # bincount over cn_indices gives carrier_count per column.
    carrier_count_per_kmer = np.bincount(cn_indices, minlength=n_kmers_total).astype(np.int32)
    print(f"  carrier_count: median={np.median(carrier_count_per_kmer):.1f}, "
          f"min={carrier_count_per_kmer.min()}, "
          f"max={carrier_count_per_kmer.max()} [{time.time()-t:.0f}s]", flush=True)

    # Step D: Apply filter: drop k-mers where global_count > carrier_count.
    # global_count[i] is occurrences of kept-kmer i across all founders' full
    # Chr1; carrier_count[i] is the # carriers within its assigned block.
    # If equal, k-mer is "block-unique" (every carrier has it once, no off-target).
    # If global > carrier: off-target genomic occurrences exist → drop.
    print("Step D: applying genomic-uniqueness filter...", flush=True)
    excess = global_count - carrier_count_per_kmer
    keep_mask = excess <= 0  # global <= carrier (should be exactly == in clean case)
    n_keep = int(keep_mask.sum())
    n_drop = int(n_kmers_total - n_keep)
    print(f"  kept k-mers: {n_keep:,} / {n_kmers_total:,} "
          f"({n_keep/n_kmers_total*100:.1f}%); dropped {n_drop:,}")
    print(f"  excess distribution: median={np.median(excess):.0f} "
          f"p95={np.percentile(excess, 95):.0f} max={excess.max()}")

    # Diagnostic: per-block kmer counts before/after
    new_block_n_kmers = np.zeros(n_blocks, dtype=np.int32)
    for b in range(n_blocks):
        ks = int(block_kmer_offsets[b])
        ke = int(block_kmer_offsets[b + 1])
        new_block_n_kmers[b] = int(keep_mask[ks:ke].sum())
    print(f"\nper-block kmer count BEFORE filter: median="
          f"{np.median(block_n_kmers):.0f}, mean={block_n_kmers.mean():.1f}")
    print(f"per-block kmer count AFTER filter:  median="
          f"{np.median(new_block_n_kmers):.0f}, mean={new_block_n_kmers.mean():.1f}")

    # Step E: Filter cn matrix and kmer_strings (per-block, memory-light).
    # We rebuild cn by remapping column indices: for each nnz entry (r, c),
    # if c is kept (keep_mask[c] = True), write it with new_col[c]. Kept
    # cols get consecutive new indices per-block (so block_kmer_offsets
    # remains contiguous in the new layout).
    print("Step E: filtering cn + kmer_strings (memory-efficient)...", flush=True)
    t = time.time()
    new_block_kmer_offsets = np.concatenate(
        [[0], np.cumsum(new_block_n_kmers)]).astype(np.int64)
    assert new_block_kmer_offsets[-1] == n_keep
    # Build per-block "old col → new col" mapping via per-block within-mask cumsum
    new_col_of_old = np.full(n_kmers_total, -1, dtype=np.int64)
    for b in range(n_blocks):
        ks = int(block_kmer_offsets[b])
        ke = int(block_kmer_offsets[b + 1])
        new_ks = int(new_block_kmer_offsets[b])
        # Within the [ks:ke) span, kept cols get consecutive new indices [new_ks, new_ks + n_keep_b)
        block_keep_mask = keep_mask[ks:ke]
        local_kept = np.where(block_keep_mask)[0]
        new_col_of_old[ks + local_kept] = new_ks + np.arange(local_kept.size, dtype=np.int64)
    # Filter nnz entries: keep where new_col_of_old[col] != -1
    new_cols = new_col_of_old[cn_indices]
    keep_nnz_mask = new_cols >= 0
    new_indices = new_cols[keep_nnz_mask]
    new_data = cn_data[keep_nnz_mask].astype(np.uint8)
    # Rebuild indptr: per-row, count how many nnz survive in that row
    # cn_indptr is current row pointers. For each row r, nnz indices are
    # cn_indices[cn_indptr[r]:cn_indptr[r+1]]. After filtering with
    # keep_nnz_mask, count survivors per row.
    nrows_total = total_haps
    # Per-nnz row-id, then bincount
    nnz_row_id = np.zeros(cn_indices.size, dtype=np.int64)
    for r in range(nrows_total):
        s, e = int(cn_indptr[r]), int(cn_indptr[r + 1])
        nnz_row_id[s:e] = r
    survivors_per_row = np.bincount(nnz_row_id[keep_nnz_mask], minlength=nrows_total)
    new_indptr = np.concatenate([[0], np.cumsum(survivors_per_row)]).astype(np.int64)
    # Decode kept k-mer strings from kept_hashes for the kept-column subset
    # (we're not relying on the original kmer_strings array which we freed)
    print(f"  decoding {n_keep:,} kept k-mer hashes to strings...", flush=True)
    kept_hash_idx = np.where(keep_mask)[0]  # global old-col indices in keep order
    # Reorder kept hashes to match the new layout (per-block, then within-block)
    # Actually new_col_of_old already gives the new layout. We want
    # new_kmer_strings[new_idx] = decoded(kept_hashes[old_idx]) for old_idx
    # where new_col_of_old[old_idx] = new_idx.
    new_kmer_strings = np.empty(n_keep, dtype=object)
    bases_arr = np.array([ord("A"), ord("C"), ord("G"), ord("T")], dtype=np.uint8)
    sel_hashes = kept_hashes[kept_hash_idx]
    new_idx_arr = new_col_of_old[kept_hash_idx]
    # Vectorize hash → string decoding
    h_arr = sel_hashes
    # Extract per-position 2-bit codes: shifts (k-1)*2, (k-2)*2, ..., 0
    decoded = np.zeros((n_keep, k), dtype=np.uint8)
    for i in range(k):
        shift = np.uint64(2 * (k - 1 - i))
        decoded[:, i] = bases_arr[(h_arr >> shift) & np.uint64(3)]
    # Convert each row to ascii string
    decoded_strs = decoded.tobytes()
    for j in range(n_keep):
        new_kmer_strings[new_idx_arr[j]] = decoded_strs[j * k:(j + 1) * k].decode("ascii")
    print(f"  cn shape: ({nrows_total}, {n_keep}), nnz: {new_indices.size:,} "
          f"[{time.time()-t:.0f}s]", flush=True)

    # Step F: Save
    print(f"Step F: saving to {args.out}...", flush=True)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        args.out,
        ecotypes=ecotypes,
        block_idx=block_idx,
        block_chrom=block_chrom,
        block_pos_start=block_pos_start,
        block_pos_end=block_pos_end,
        block_n_haps=block_n_haps,
        block_n_kmers=new_block_n_kmers,
        block_kmer_offsets=new_block_kmer_offsets,
        block_hap_offsets=block_hap_offsets,  # unchanged; rows preserved
        kmer_strings=np.asarray(new_kmer_strings, dtype=object),
        cn_indptr=new_indptr,
        cn_indices=new_indices.astype(np.int64),
        cn_data=new_data,
        founder_to_class=founder_to_class,
        k=k,
    )
    print(f"saved. n_kept={n_keep:,}", flush=True)


if __name__ == "__main__":
    main()
