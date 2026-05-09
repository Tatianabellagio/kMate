"""Build per-block class × kmer cn from CLEAN reconstructed haplotype sequences.

For each BigLD block + each cn_var-derived class, reconstruct the class's
haplotype sequence by walking TAIR10 reference and applying the variants
the class carries (per cn_var). Hash canonical k-mers from this clean
sequence. NO per-founder FASTA noise — `cn` is bias-free by construction.

Inputs:
  --block-index    BigLD block coords
  --cn-var         cn_var .npz (founder × variant)
  --cn-var-meta    cn_var .meta.npz
  --vcf            source VCF mapped to cn_var (one record per cn_var entry,
                   in VCF iteration order; `bcftools view -H ... | wc -l`
                   should equal cn_var.shape[1])
  --ref            TAIR10 reference FASTA (Chr1..5 named)
  --fastas-dir     Not used for k-mers (clean reconstruction); only for
                   compatibility with EM driver (founder_to_class structure)
  --out            output .npz

Pipeline per block:
  1. cn_var-dedup founders by their genotype profile in this block → classes
  2. For each cn_var record in block, get (ref_seq, alt_seq) from VCF
  3. For each class, walk reference applying variants per class genotype →
     clean haplotype sequence
  4. Hash canonical k-mers → class k-mer set
  5. Discriminative filter (cc in [1, n_classes-1])
  6. Cross-block dedup (k-mer in only 1 block's discriminative set)
"""
from __future__ import annotations
import argparse
import os
import sys
import time
from collections import defaultdict
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pysam
from scipy.sparse import csr_matrix, load_npz

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from build_block_haplotype_cn import (
    canonical_kmer_hashes, canonical_kmer_hashes_with_counts, hash_to_kmer,
)


def reconstruct_class_sequence(
    ref_seq: bytes,           # block-range reference, 0-based, length L
    block_pos_start: int,     # 1-based, inclusive
    variants_in_block: list,  # list of (cn_var_idx, pos, ref_seq, alt_seq, class_carries: bool)
) -> str:
    """Walk reference from start to end, applying variants the class carries.

    variants_in_block must be sorted by pos. For each variant:
      - If class carries it: emit alt_seq, advance position by len(ref_seq)
      - If not: emit ref bases until past the variant, then continue
    Overlapping/adjacent variants handled by skipping any whose ref-span starts
    before the current emit position (keeps the first overlapping one).
    """
    out = bytearray()
    cursor = 0  # 0-based offset into ref_seq we've emitted up to (exclusive)
    L = len(ref_seq)
    for cnv_idx, pos, vref, valt, carry in variants_in_block:
        # 0-based offset of variant start in ref_seq
        var_off = pos - block_pos_start
        if var_off < 0 or var_off >= L:
            continue
        if var_off < cursor:
            # Overlap with a previous variant we already applied; skip
            continue
        # Emit reference up to the variant start
        if cursor < var_off:
            out.extend(ref_seq[cursor:var_off])
            cursor = var_off
        # Apply variant
        ref_end = var_off + len(vref)
        if carry:
            out.extend(valt)
            cursor = ref_end
        else:
            # Emit the ref allele (may be shorter than the variant's local span)
            out.extend(ref_seq[var_off:min(ref_end, L)])
            cursor = ref_end
    # Tail
    if cursor < L:
        out.extend(ref_seq[cursor:L])
    return out.decode("ascii", errors="replace")


def process_block_clean(args_tuple):
    """Build per-block class × kmer cn using clean reconstructed sequences.

    With dose_aware=True, cn[class, kmer] = number of times kmer appears in
    class's reconstructed sequence (PanGenie-MODEL_SPEC's "expected k-mer
    copy number" per class) instead of binary 0/1. Helps repeat-affected
    positions where a k-mer has multiple genomic copies in a haplotype.
    """
    args_tuple = list(args_tuple)
    (block_idx, chrom, pos_start, pos_end,
     ecotypes, k,
     cnvar_block,
     variants_in_block,
     ref_seq_bytes) = args_tuple[:9]
    dose_aware = args_tuple[9] if len(args_tuple) > 9 else False
    n_eco = len(ecotypes)

    if cnvar_block.shape[1] == 0:
        founder_to_class = np.zeros(n_eco, dtype=np.int32)
        return (block_idx, np.empty(0, dtype=np.uint64),
                csr_matrix((1, 0), dtype=np.uint8), founder_to_class)

    # cn_var dedup
    seen: dict[bytes, int] = {}
    founder_to_class = np.empty(n_eco, dtype=np.int32)
    members_for_class: list[list[int]] = []
    class_genotypes: list[np.ndarray] = []
    for f in range(n_eco):
        key = cnvar_block[f].tobytes()
        cid = seen.get(key)
        if cid is None:
            cid = len(members_for_class)
            seen[key] = cid
            members_for_class.append([f])
            class_genotypes.append(cnvar_block[f].copy())
        else:
            members_for_class[cid].append(f)
        founder_to_class[f] = cid
    n_classes = len(members_for_class)

    if n_classes <= 1:
        return (block_idx, np.empty(0, dtype=np.uint64),
                csr_matrix((n_classes, 0), dtype=np.uint8), founder_to_class)

    # Reconstruct clean sequence per class and hash k-mers
    class_hashes: list[np.ndarray] = []
    class_counts: list[np.ndarray] = []   # only used in dose_aware mode
    for c in range(n_classes):
        gt = class_genotypes[c]
        vars_list = []
        for v_local_idx, pos, vref, valt in variants_in_block:
            carry = bool(gt[v_local_idx])
            vars_list.append((v_local_idx, pos, vref, valt, carry))
        vars_list.sort(key=lambda x: x[1])
        seq = reconstruct_class_sequence(ref_seq_bytes, pos_start, vars_list)
        if dose_aware:
            h_arr, c_arr = canonical_kmer_hashes_with_counts(seq, k=k)
            class_hashes.append(h_arr)
            class_counts.append(c_arr)
        else:
            h_arr = canonical_kmer_hashes(seq, k=k)
            class_hashes.append(h_arr)

    # Build cn[class, kmer] with discriminative filter
    if all(h.size == 0 for h in class_hashes):
        return (block_idx, np.empty(0, dtype=np.uint64),
                csr_matrix((n_classes, 0), dtype=np.uint8), founder_to_class)
    union_hashes = np.unique(np.concatenate(class_hashes))

    rows: list[np.ndarray] = []
    cols: list[np.ndarray] = []
    data_chunks: list[np.ndarray] = []
    for c, h in enumerate(class_hashes):
        if h.size == 0:
            continue
        idx = np.searchsorted(union_hashes, h)
        cols.append(idx)
        rows.append(np.full(idx.size, c, dtype=np.int32))
        if dose_aware:
            cnts = class_counts[c]
            # Cap doses at 255 to fit in uint8
            data_chunks.append(np.minimum(cnts, 255).astype(np.uint8))
        else:
            data_chunks.append(np.ones(idx.size, dtype=np.uint8))
    rows_arr = np.concatenate(rows) if rows else np.empty(0, dtype=np.int32)
    cols_arr = np.concatenate(cols) if cols else np.empty(0, dtype=np.int64)
    data = (np.concatenate(data_chunks) if data_chunks else
            np.empty(0, dtype=np.uint8))
    cn = csr_matrix((data, (rows_arr, cols_arr)),
                    shape=(n_classes, union_hashes.size), dtype=np.uint8)

    # Discriminative filter: count CARRIERS (binary) not dose, even in dose mode
    if dose_aware:
        cn_binary = (cn > 0).astype(np.int8)
        col_sum = np.asarray(cn_binary.sum(axis=0)).flatten()
    else:
        col_sum = np.asarray(cn.sum(axis=0)).flatten()
    discrim = (col_sum > 0) & (col_sum < n_classes)
    if not discrim.any():
        return (block_idx, np.empty(0, dtype=np.uint64),
                csr_matrix((n_classes, 0), dtype=np.uint8), founder_to_class)
    keep_idx = np.where(discrim)[0]
    return (block_idx, union_hashes[keep_idx],
            cn[:, keep_idx].tocsr(), founder_to_class)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--block-index", required=True)
    ap.add_argument("--cn-var", required=True)
    ap.add_argument("--cn-var-meta", required=True)
    ap.add_argument("--vcf", required=True,
                    help="source VCF for cn_var (e.g. merged_v2_chr.vcf.gz)")
    ap.add_argument("--ref", required=True,
                    help="TAIR10 reference FASTA (Chr1..5 named)")
    ap.add_argument("--out", required=True)
    ap.add_argument("--k", type=int, default=31)
    ap.add_argument("--chrom-prefix", default="Chr")
    ap.add_argument("--chrom-filter", default=None)
    ap.add_argument("--max-blocks", type=int, default=None)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--dose-aware", action="store_true",
                    help="cn[class,kmer] = count of kmer in class sequence (not 0/1)")
    args = ap.parse_args()

    print(f"Loading block index from {args.block_index}...", flush=True)
    z = np.load(args.block_index, allow_pickle=True)
    ecotypes = np.asarray(z["ecotypes"]).astype(str)
    block_chrom = np.asarray(z["block_chrom"]).astype(str)
    block_pos_start = np.asarray(z["block_pos_start"]).astype(np.int64)
    block_pos_end = np.asarray(z["block_pos_end"]).astype(np.int64)

    print(f"Loading cn_var from {args.cn_var}...", flush=True)
    cn_var = load_npz(args.cn_var)
    cv_meta = np.load(args.cn_var_meta, allow_pickle=True)
    cv_founders = np.asarray(cv_meta["founders"]).astype(str)
    cv_chrom = np.asarray(cv_meta["chrom"]).astype(str)
    cv_pos = np.asarray(cv_meta["pos"]).astype(np.int64)
    cv_ref_len = np.asarray(cv_meta["ref_len"]).astype(np.int64)
    cv_alt_len = np.asarray(cv_meta["alt_len"]).astype(np.int64)
    print(f"  cn_var: {cn_var.shape}, n_records: {cn_var.shape[1]:,}", flush=True)

    # Reorder cn_var rows to match block-index ecotype order
    if not np.array_equal(cv_founders, ecotypes):
        print("  reordering cn_var rows to match block-index ecotype order...",
              flush=True)
        rev = {e: i for i, e in enumerate(cv_founders)}
        order = np.array([rev[e] for e in ecotypes])
        cn_var = cn_var[order, :]

    block_chrom_full = np.array(
        [f"{args.chrom_prefix}{c}" if not c.startswith(args.chrom_prefix) else c
         for c in block_chrom], dtype=object)

    n_blocks_total = len(block_chrom_full)
    sel = np.ones(n_blocks_total, dtype=bool)
    if args.chrom_filter:
        sel &= (block_chrom_full == args.chrom_filter)
    if args.max_blocks:
        idx_sel = np.where(sel)[0][:args.max_blocks]
        sel[:] = False
        sel[idx_sel] = True
    keep_idx = np.where(sel)[0]
    print(f"Processing {len(keep_idx)} / {n_blocks_total} blocks", flush=True)

    # Open VCF + reference
    print("Opening VCF + reference for sequence extraction...", flush=True)
    vcf = pysam.VariantFile(args.vcf)
    ref = pysam.FastaFile(args.ref)

    # cn_var meta uses chrom labels matching VCF (Chr1, ..., Chr5).
    # We assume cn_var record order = VCF iteration order.
    # For per-block fetching, we use vcf.fetch(chrom, start-1, end) to get
    # records in range, and align to cn_var via (chrom, pos, ref_len, alt_len)
    # matching. Where multiple alts share a position+ref_len+alt_len, we
    # match in iteration order using a per-position counter.

    # Build a per-(chrom, pos) → list of cn_var indices index for fast lookup
    print("Indexing cn_var by (chrom, pos) for fast block lookup...", flush=True)
    t = time.time()
    cv_pos_index: dict[tuple, list[int]] = defaultdict(list)
    for i in range(len(cv_pos)):
        cv_pos_index[(cv_chrom[i], int(cv_pos[i]))].append(i)
    print(f"  indexed {sum(len(v) for v in cv_pos_index.values()):,} entries "
          f"[{time.time()-t:.0f}s]", flush=True)

    # Pre-fetch per-block reference + variants. Pass to workers.
    print("Pre-fetching per-block ref + variants from VCF...", flush=True)
    t = time.time()
    jobs = []
    for b in keep_idx:
        chrom_b = str(block_chrom_full[b])
        start_b = int(block_pos_start[b])
        end_b = int(block_pos_end[b])
        # Reference for the block
        ref_seq = ref.fetch(chrom_b, start_b - 1, end_b)
        ref_seq_bytes = ref_seq.encode("ascii")
        # Variants in block — fetch via tabix; iterate in VCF order
        variants_in_block = []  # (v_local_idx, pos, ref_bytes, alt_bytes)
        cnvar_records = []      # cn_var indices in same order
        # Per-position counter for matching multi-allelic records
        pos_counter: dict[tuple, int] = defaultdict(int)
        for rec in vcf.fetch(chrom_b, start_b - 1, end_b):
            if not rec.alts:
                continue
            rpos = rec.pos
            r_ref = rec.ref
            r_alt = rec.alts[0]  # build_cn_var.py only used alts[0]
            key = (chrom_b, rpos)
            occ = pos_counter[key]
            pos_counter[key] += 1
            cnvar_indices_at_pos = cv_pos_index.get(key, [])
            if occ >= len(cnvar_indices_at_pos):
                # Mismatch — VCF has more records at this pos than cn_var; skip
                continue
            cv_idx = cnvar_indices_at_pos[occ]
            # Sanity: ref_len/alt_len should match
            if cv_ref_len[cv_idx] != len(r_ref) or cv_alt_len[cv_idx] != len(r_alt):
                continue  # skip if mismatch
            local_idx = len(cnvar_records)
            cnvar_records.append(cv_idx)
            variants_in_block.append(
                (local_idx, rpos, r_ref.encode("ascii"), r_alt.encode("ascii"))
            )
        # Slice cn_var to in-block records, in the SAME order as variants_in_block
        if cnvar_records:
            cnvar_block = cn_var[:, cnvar_records].toarray().astype(np.uint8)
        else:
            cnvar_block = np.zeros((cn_var.shape[0], 0), dtype=np.uint8)
        jobs.append((
            int(b), chrom_b, start_b, end_b,
            ecotypes, args.k,
            cnvar_block, variants_in_block, ref_seq_bytes,
            bool(args.dose_aware),
        ))
    vcf.close()
    ref.close()
    sizes = [j[6].shape[1] for j in jobs]
    print(f"  pre-fetched {len(jobs)} blocks in {time.time()-t:.0f}s; "
          f"per-block n_var: median={int(np.median(sizes))}, max={max(sizes)}",
          flush=True)

    # Worker pool
    t0 = time.time()
    results: dict[int, tuple] = {}
    if args.threads > 1:
        with Pool(args.threads) as pool:
            for i, res in enumerate(pool.imap_unordered(
                    process_block_clean, jobs, chunksize=4)):
                bidx, kmers, cn, f2c = res
                results[bidx] = (kmers, cn, f2c)
                if (i + 1) % 200 == 0:
                    print(f"  {i+1}/{len(jobs)} blocks "
                          f"[{time.time()-t0:.0f}s]", flush=True)
    else:
        for i, j in enumerate(jobs):
            bidx, kmers, cn, f2c = process_block_clean(j)
            results[bidx] = (kmers, cn, f2c)
            if (i + 1) % 200 == 0:
                print(f"  {i+1}/{len(jobs)} blocks "
                      f"[{time.time()-t0:.0f}s]", flush=True)
    print(f"per-block hashing done [{time.time()-t0:.0f}s]", flush=True)

    # Cross-block dedup
    print("global cross-block dedup...", flush=True)
    t = time.time()
    all_hashes = []
    for b in keep_idx:
        h = results[int(b)][0]
        if h.size:
            all_hashes.append(h)
    if all_hashes:
        concat = np.concatenate(all_hashes)
        uniq, counts = np.unique(concat, return_counts=True)
        unique_set_view = uniq[counts == 1]
    else:
        unique_set_view = np.empty(0, dtype=np.uint64)
        uniq = unique_set_view
    n_global = uniq.size
    n_unique_global = int(unique_set_view.size)
    n_dropped = n_global - n_unique_global
    print(f"  global k-mers: {n_global:,}; unique to one block: "
          f"{n_unique_global:,}; dropped: {n_dropped:,} "
          f"({n_dropped/max(1,n_global)*100:.1f}%) [{time.time()-t:.0f}s]",
          flush=True)

    # Build unified output
    print("building output structures...", flush=True)
    t = time.time()
    block_idx_arr = []
    block_kmer_offsets = [0]
    kmer_hashes_all: list[np.ndarray] = []
    block_n_haps: list[int] = []
    block_n_kmers: list[int] = []
    block_chrom_kept: list[str] = []
    block_pos_start_kept: list[int] = []
    block_pos_end_kept: list[int] = []
    founder_to_class_kept: list[np.ndarray] = []
    cn_local_blocks: list[csr_matrix] = []

    for b in keep_idx:
        kmer_hashes, cn_full, f2c = results[int(b)]
        n_classes_b = cn_full.shape[0]
        if kmer_hashes.size == 0:
            keep_local_idx = np.empty(0, dtype=np.int64)
        else:
            pos = np.searchsorted(unique_set_view, kmer_hashes)
            in_set = (pos < unique_set_view.size) & \
                     (unique_set_view[np.clip(pos, 0, unique_set_view.size - 1)] == kmer_hashes)
            keep_local_idx = np.where(in_set)[0]
        kept_hashes = (kmer_hashes[keep_local_idx] if keep_local_idx.size
                       else np.empty(0, dtype=np.uint64))
        if cn_full.shape[1] > 0 and keep_local_idx.size > 0:
            cn_kept = cn_full[:, keep_local_idx].tocsr()
        else:
            cn_kept = csr_matrix((n_classes_b, 0), dtype=np.uint8)
        block_idx_arr.append(int(b))
        kmer_hashes_all.append(kept_hashes)
        block_kmer_offsets.append(block_kmer_offsets[-1] + kept_hashes.size)
        block_n_haps.append(int(n_classes_b))
        block_n_kmers.append(int(kept_hashes.size))
        block_chrom_kept.append(str(block_chrom_full[b]))
        block_pos_start_kept.append(int(block_pos_start[b]))
        block_pos_end_kept.append(int(block_pos_end[b]))
        founder_to_class_kept.append(f2c)
        cn_local_blocks.append(cn_kept)

    kmer_hashes_concat = (np.concatenate(kmer_hashes_all)
                          if kmer_hashes_all else np.empty(0, dtype=np.uint64))
    print(f"  decoding {kmer_hashes_concat.size:,} k-mer hashes to strings...",
          flush=True)
    kmer_strings_all = [hash_to_kmer(int(h), args.k) for h in kmer_hashes_concat]

    # Concat per-block sparse cn matrices into one big CSR
    print(f"  concatenating cn rows for {len(cn_local_blocks)} blocks...",
          flush=True)
    rows_per_block = block_n_haps
    total_rows = sum(rows_per_block)
    total_cols = block_kmer_offsets[-1]
    indptr_list = [0]
    indices_chunks: list[np.ndarray] = []
    nnz_running = 0
    for bi, cn_b in enumerate(cn_local_blocks):
        n_haps_b = cn_b.shape[0]
        col_offset = block_kmer_offsets[bi]
        cb_indptr = cn_b.indptr
        cb_indices = cn_b.indices
        if cb_indices.size:
            indices_chunks.append(cb_indices + col_offset)
        for h in range(1, n_haps_b + 1):
            indptr_list.append(int(cb_indptr[h]) + nnz_running)
        nnz_running += int(cb_indptr[-1])
    cn_indptr = np.asarray(indptr_list, dtype=np.int64)
    cn_indices = (np.concatenate(indices_chunks).astype(np.int64)
                  if indices_chunks else np.empty(0, dtype=np.int64))
    cn_data = np.ones(cn_indices.size, dtype=np.uint8)
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
        founder_to_class=np.asarray(founder_to_class_kept, dtype=np.int32),
        k=args.k,
    )
    print(f"done. total wall: {time.time()-t0:.0f}s", flush=True)
    print(f"  total k-mers: {total_cols:,}")
    print(f"  total class-rows: {total_rows:,}")
    print(f"  median classes/block: {int(np.median(block_n_haps))}")
    print(f"  median kmers/block:   {int(np.median(block_n_kmers))}")


if __name__ == "__main__":
    main()
