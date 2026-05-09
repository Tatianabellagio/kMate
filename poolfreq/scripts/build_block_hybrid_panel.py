"""Build hybrid per-block founder-level cn matrix combining:
  - clean_cc5 panel (haplotype-block class k-mers, scattered class→founder)
  - cn_full_231_v2 (PanGenie bubble-level k-mers, founder-level natively)

Output: a single npz that the existing block_haplotype_em.py can load with
inference_unit="founder". Each per-block cn has rows=231 founders, cols=
(union of clean_cc5 + PanGenie k-mers in that block range).
"""
from __future__ import annotations
import argparse
import time
from pathlib import Path
import numpy as np
from scipy.sparse import csr_matrix, load_npz, vstack as sp_vstack, hstack as sp_hstack

import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from build_block_haplotype_cn import canonical_kmer_hashes, hash_to_kmer


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--clean-npz", required=True,
                    help="clean_cc5 panel (chr1_full_clean_cc5.npz)")
    ap.add_argument("--cn-full-prefix", required=True,
                    help="cn_full prefix; expects <prefix>_<chrom>.cn.npz + meta.npz")
    ap.add_argument("--chrom", default="Chr1")
    ap.add_argument("--out", required=True)
    ap.add_argument("--k", type=int, default=31)
    args = ap.parse_args()

    # Load clean panel
    print(f"Loading {args.clean_npz}...", flush=True)
    z = np.load(args.clean_npz, allow_pickle=True)
    ecotypes = np.asarray(z["ecotypes"]).astype(str)
    block_chrom = np.asarray(z["block_chrom"]).astype(str)
    block_pos_start = np.asarray(z["block_pos_start"]).astype(np.int64)
    block_pos_end = np.asarray(z["block_pos_end"]).astype(np.int64)
    block_idx_arr = np.asarray(z["block_idx"]).astype(np.int64)
    clean_block_n_haps = np.asarray(z["block_n_haps"]).astype(np.int32)
    clean_block_n_kmers = np.asarray(z["block_n_kmers"]).astype(np.int32)
    clean_block_kmer_offsets = np.asarray(z["block_kmer_offsets"]).astype(np.int64)
    clean_block_hap_offsets = np.asarray(z["block_hap_offsets"]).astype(np.int64)
    clean_kmer_strings = z["kmer_strings"]
    clean_cn_indptr = np.asarray(z["cn_indptr"]).astype(np.int64)
    clean_cn_indices = np.asarray(z["cn_indices"]).astype(np.int64)
    clean_cn_data = np.asarray(z["cn_data"]).astype(np.uint8)
    founder_to_class = np.asarray(z["founder_to_class"]).astype(np.int32)
    n_eco = len(ecotypes)
    n_blocks = len(block_idx_arr)
    total_clean_haps = int(clean_block_hap_offsets[-1])
    total_clean_kmers = int(clean_block_kmer_offsets[-1])
    print(f"  blocks={n_blocks}, n_eco={n_eco}, total_clean_kmers={total_clean_kmers:,}",
          flush=True)
    cn_clean = csr_matrix((clean_cn_data, clean_cn_indices, clean_cn_indptr),
                          shape=(total_clean_haps, total_clean_kmers), dtype=np.uint8)

    # Load PanGenie cn_full for this chrom
    print(f"Loading PanGenie cn_full for {args.chrom}...", flush=True)
    cn_full = load_npz(f"{args.cn_full_prefix}_{args.chrom}.cn.npz")
    pmeta = np.load(f"{args.cn_full_prefix}_{args.chrom}.meta.npz", allow_pickle=True)
    pkmer_index = np.asarray(pmeta["kmer_index"]).astype(str)
    bubble_chrom = np.asarray(pmeta["bubble_chrom"]).astype(str)
    bubble_start = np.asarray(pmeta["bubble_start"]).astype(np.int64)
    bubble_end = np.asarray(pmeta["bubble_end"]).astype(np.int64)
    print(f"  cn_full shape: {cn_full.shape}, n_panel_kmers: {len(pkmer_index):,}",
          flush=True)

    # Build per-block lookup of PanGenie k-mer indices: those whose bubble pos
    # falls within block [pos_start, pos_end].
    print("Indexing PanGenie k-mers per block...", flush=True)
    t = time.time()
    # For each PanGenie k-mer, its bubble position
    panel_pos = (bubble_start + bubble_end) / 2.0
    sort_perm = np.argsort(panel_pos)
    sorted_pos = panel_pos[sort_perm]
    print(f"  panel kmers sorted [{time.time()-t:.0f}s]", flush=True)

    # For each block, find panel k-mer indices in [start, end]
    panel_kmers_per_block: list[np.ndarray] = []
    for b in range(n_blocks):
        if str(block_chrom[b]) != args.chrom:
            panel_kmers_per_block.append(np.empty(0, dtype=np.int64))
            continue
        lo = np.searchsorted(sorted_pos, float(block_pos_start[b]), side="left")
        hi = np.searchsorted(sorted_pos, float(block_pos_end[b]), side="right")
        panel_kmers_per_block.append(sort_perm[lo:hi])
    n_panel_total = sum(p.size for p in panel_kmers_per_block)
    print(f"  panel kmers in blocks: {n_panel_total:,} "
          f"(median per block: {int(np.median([p.size for p in panel_kmers_per_block]))}) "
          f"[{time.time()-t:.0f}s]", flush=True)

    # Build hybrid per-block cn at FOUNDER level (rows=231, cols=clean+panel kmers).
    # Per block:
    #   clean cols: scatter cn_clean rows from class to founder via founder_to_class
    #   panel cols: cn_full sliced to panel_kmers_per_block[b] (already founder-level)
    print("Building hybrid per-block cn at founder level...", flush=True)
    t = time.time()
    hybrid_block_n_kmers = []
    hybrid_block_kmer_offsets = [0]
    hybrid_kmer_strings_all: list[str] = []
    cn_blocks_founder: list[csr_matrix] = []

    for b in range(n_blocks):
        # Clean part
        ks = int(clean_block_kmer_offsets[b])
        ke = int(clean_block_kmer_offsets[b + 1])
        n_clean_b = ke - ks
        hs = int(clean_block_hap_offsets[b])
        he = int(clean_block_hap_offsets[b + 1])
        if n_clean_b > 0:
            cn_clean_block = cn_clean[hs:he, ks:ke]      # (n_class_b, n_clean_b)
            f2c_b = founder_to_class[b]                   # (n_eco,) int32 in [0, n_class_b)
            # Scatter: cn_clean_founder[F, K] = cn_clean_block[f2c_b[F], K]
            cn_clean_founder = cn_clean_block[f2c_b]      # (n_eco, n_clean_b) sparse
        else:
            cn_clean_founder = csr_matrix((n_eco, 0), dtype=np.uint8)
        # Panel part
        panel_kmer_idx = panel_kmers_per_block[b]
        if panel_kmer_idx.size > 0:
            cn_panel_block = cn_full[:, panel_kmer_idx].tocsr().astype(np.uint8)
        else:
            cn_panel_block = csr_matrix((n_eco, 0), dtype=np.uint8)
        # Concatenate columns: clean | panel
        cn_combined = sp_hstack([cn_clean_founder, cn_panel_block]).tocsr()
        cn_blocks_founder.append(cn_combined)
        # Update kmer strings
        clean_strs = list(clean_kmer_strings[ks:ke])
        panel_strs = [str(pkmer_index[i]) for i in panel_kmer_idx]
        block_strs = clean_strs + panel_strs
        hybrid_kmer_strings_all.extend(block_strs)
        hybrid_block_n_kmers.append(len(block_strs))
        hybrid_block_kmer_offsets.append(hybrid_block_kmer_offsets[-1] + len(block_strs))
        if (b + 1) % 1000 == 0:
            print(f"  {b+1}/{n_blocks} blocks [{time.time()-t:.0f}s]", flush=True)
    print(f"per-block hybrid construction done [{time.time()-t:.0f}s]", flush=True)

    # Save in the same format as our other panels.
    # Note: at founder level, there's no "class" structure. We set
    # block_n_haps = n_eco (231) for all blocks, founder_to_class = identity,
    # and cn_indptr/indices/data is the founder-level CSR.
    print("concatenating cn matrices into single CSR...", flush=True)
    t = time.time()
    block_hap_offsets = np.arange(n_blocks + 1, dtype=np.int64) * n_eco
    indptr_list = [0]
    indices_chunks: list[np.ndarray] = []
    nnz_running = 0
    for bi, cn_b in enumerate(cn_blocks_founder):
        col_offset = hybrid_block_kmer_offsets[bi]
        cb_indptr = cn_b.indptr
        cb_indices = cn_b.indices
        if cb_indices.size:
            indices_chunks.append(cb_indices + col_offset)
        for h in range(1, n_eco + 1):
            indptr_list.append(int(cb_indptr[h]) + nnz_running)
        nnz_running += int(cb_indptr[-1])
    cn_indptr = np.asarray(indptr_list, dtype=np.int64)
    cn_indices = (np.concatenate(indices_chunks).astype(np.int64)
                  if indices_chunks else np.empty(0, dtype=np.int64))
    cn_data = np.ones(cn_indices.size, dtype=np.uint8)
    print(f"  shape: rows={n_blocks * n_eco}, cols={hybrid_block_kmer_offsets[-1]}, "
          f"nnz={cn_indices.size:,} [{time.time()-t:.0f}s]", flush=True)

    # founder_to_class: identity (each founder is its own "class" at this stage)
    f2c_identity = np.broadcast_to(
        np.arange(n_eco, dtype=np.int32)[None, :], (n_blocks, n_eco)
    ).copy()

    print(f"saving to {args.out}...", flush=True)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        args.out,
        ecotypes=ecotypes,
        block_idx=block_idx_arr,
        block_chrom=block_chrom,
        block_pos_start=block_pos_start,
        block_pos_end=block_pos_end,
        block_n_haps=np.full(n_blocks, n_eco, dtype=np.int32),
        block_n_kmers=np.asarray(hybrid_block_n_kmers, dtype=np.int32),
        block_kmer_offsets=np.asarray(hybrid_block_kmer_offsets, dtype=np.int64),
        block_hap_offsets=block_hap_offsets,
        kmer_strings=np.asarray(hybrid_kmer_strings_all, dtype=object),
        cn_indptr=cn_indptr,
        cn_indices=cn_indices,
        cn_data=cn_data,
        founder_to_class=f2c_identity,
        k=args.k,
    )
    print("DONE", flush=True)
    print(f"  total kmers: {hybrid_block_kmer_offsets[-1]:,}")
    print(f"  median kmers/block: {int(np.median(hybrid_block_n_kmers))}")


if __name__ == "__main__":
    main()
