"""Filter block_haplotype_cn npz to keep only k-mers with carrier_count >= threshold.

Drops singleton fingerprint k-mers (mostly imputation noise + private SNPs)
while keeping multi-carrier "allele-like" k-mers that pool evidence across
founders sharing alleles. Aim is to mimic PanGenie's allele-grouping
structure but starting from FASTA-derived k-mers rather than VCF bubbles.

Usage:
    python filter_block_haplotype_cn_min_carrier.py --in chr1_full_genuniq.npz \\
        --min-cc 5 --out chr1_full_genuniq_cc5.npz
"""
from __future__ import annotations
import argparse
import time
from pathlib import Path

import numpy as np


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="in_npz", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--min-cc", type=int, default=5,
                    help="minimum carrier_count (#class-rows carrying k-mer) to keep")
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
    kmer_strings = np.asarray(z["kmer_strings"])
    cn_indptr = np.asarray(z["cn_indptr"]).astype(np.int64)
    cn_indices = np.asarray(z["cn_indices"]).astype(np.int64)
    cn_data = np.asarray(z["cn_data"]).astype(np.uint8)
    founder_to_class = np.asarray(z["founder_to_class"]).astype(np.int32)
    k = int(z["k"])

    n_blocks = len(block_idx)
    n_kmers_total = int(block_kmer_offsets[-1])
    total_haps = int(block_hap_offsets[-1])
    print(f"  blocks={n_blocks}, total_kmers={n_kmers_total:,}, "
          f"total_haps={total_haps:,}", flush=True)

    # Compute carrier_count per k-mer (# class-rows that carry each k-mer)
    print("computing carrier_count per kmer...", flush=True)
    t = time.time()
    carrier_count = np.bincount(cn_indices, minlength=n_kmers_total).astype(np.int32)
    print(f"  carrier_count: median={np.median(carrier_count):.1f}, "
          f"mean={carrier_count.mean():.1f}, max={carrier_count.max()} "
          f"[{time.time()-t:.0f}s]", flush=True)

    # Filter
    keep_mask = carrier_count >= args.min_cc
    n_keep = int(keep_mask.sum())
    print(f"  keep cc >= {args.min_cc}: {n_keep:,} / {n_kmers_total:,} "
          f"({n_keep/n_kmers_total*100:.1f}%)")

    # Per-block kmer counts
    new_block_n_kmers = np.zeros(n_blocks, dtype=np.int32)
    for b in range(n_blocks):
        ks = int(block_kmer_offsets[b])
        ke = int(block_kmer_offsets[b + 1])
        new_block_n_kmers[b] = int(keep_mask[ks:ke].sum())
    print(f"per-block kmer count BEFORE: median={np.median(block_n_kmers):.0f}, mean={block_n_kmers.mean():.1f}")
    print(f"per-block kmer count AFTER:  median={np.median(new_block_n_kmers):.0f}, mean={new_block_n_kmers.mean():.1f}")
    print(f"blocks emptied (n_kmers=0):  {(new_block_n_kmers == 0).sum()} / {n_blocks}")

    # Build new column index mapping (per-block contiguous)
    print("rebuilding cn (memory-efficient)...", flush=True)
    t = time.time()
    new_block_kmer_offsets = np.concatenate(
        [[0], np.cumsum(new_block_n_kmers)]).astype(np.int64)
    assert new_block_kmer_offsets[-1] == n_keep

    new_col_of_old = np.full(n_kmers_total, -1, dtype=np.int64)
    for b in range(n_blocks):
        ks = int(block_kmer_offsets[b])
        ke = int(block_kmer_offsets[b + 1])
        new_ks = int(new_block_kmer_offsets[b])
        block_keep_mask = keep_mask[ks:ke]
        local_kept = np.where(block_keep_mask)[0]
        new_col_of_old[ks + local_kept] = new_ks + np.arange(local_kept.size, dtype=np.int64)

    # Filter cn nnz entries
    new_cols = new_col_of_old[cn_indices]
    keep_nnz_mask = new_cols >= 0
    new_indices = new_cols[keep_nnz_mask]
    new_data = cn_data[keep_nnz_mask].astype(np.uint8)

    # Rebuild indptr per row
    nrows_total = total_haps
    nnz_row_id = np.zeros(cn_indices.size, dtype=np.int64)
    for r in range(nrows_total):
        s, e = int(cn_indptr[r]), int(cn_indptr[r + 1])
        nnz_row_id[s:e] = r
    survivors_per_row = np.bincount(nnz_row_id[keep_nnz_mask], minlength=nrows_total)
    new_indptr = np.concatenate([[0], np.cumsum(survivors_per_row)]).astype(np.int64)

    # Filter kmer_strings
    keep_global_idx = np.where(keep_mask)[0]
    # Reorder kept kmer_strings to match new layout (per-block, then within-block)
    new_kmer_strings = np.empty(n_keep, dtype=object)
    new_idx_arr = new_col_of_old[keep_global_idx]
    for j, old_i in enumerate(keep_global_idx):
        new_kmer_strings[new_idx_arr[j]] = kmer_strings[old_i]
    print(f"  cn shape: ({nrows_total}, {n_keep}), nnz: {new_indices.size:,} "
          f"[{time.time()-t:.0f}s]", flush=True)

    # Save
    print(f"saving to {args.out}...", flush=True)
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
        block_hap_offsets=block_hap_offsets,
        kmer_strings=new_kmer_strings,
        cn_indptr=new_indptr,
        cn_indices=new_indices.astype(np.int64),
        cn_data=new_data,
        founder_to_class=founder_to_class,
        k=k,
    )
    print(f"saved. n_kept={n_keep:,}", flush=True)


if __name__ == "__main__":
    main()
