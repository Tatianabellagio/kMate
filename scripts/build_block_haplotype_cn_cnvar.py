"""Build per-block haplotype-level k-mer cn matrix using **cn_var-based**
founder dedup (instead of FASTA-slice dedup).

Per BigLD block:
  1. Slice cn_var to records within the block's bp range.
  2. Dedup founders by their variant-genotype profile in this slice
     (SNPs + SVs, panel-VCF only — excludes per-base FASTA noise).
  3. For each class, pick a representative founder and hash its full FASTA
     slice → k-mer set.
  4. Build cn[class, kmer] from those k-mer sets.
  5. Within-block discriminative filter (carriers in [1, n_classes-1]).
  6. Cross-block dedup of k-mers shared across blocks (same as the
     FASTA-dedup builder).

Why cn_var dedup:
  - SNPs + SVs are first-class in panel; assembly noise is not in cn_var
    so it can't fragment classes.
  - Founders with identical real biology in a block collapse to one class
    even if their FASTAs differ at noise positions.
  - At founder-level (after scatter), k-mers in a multi-founder class
    have carrier_count >= class size — natural multi-carrier structure
    that gives the EM Poisson-SNR a fighting chance.
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
    canonical_kmer_hashes, hash_to_kmer, _BASE_LUT, _open_fa,
)


def process_block_cnvar(args_tuple):
    """Build per-block class × kmer cn using cn_var-derived dedup.

    Returns (block_idx, kmer_hashes_kept (uint64 array),
             cn_kept_csr (n_classes, K_kept), founder_to_class int32 (N_eco,)).

    args_tuple now also accepts variant_positions for variant-anchored
    k-mer masking (N-out positions outside ±(k-1) of any cn_var record).
    """
    args_tuple = list(args_tuple)
    (block_idx, chrom, pos_start, pos_end,
     ecotypes, fastas_dir, k,
     cnvar_slice_path, in_block_record_idx) = args_tuple[:9]
    variant_positions_in_block = args_tuple[9] if len(args_tuple) > 9 else None
    variant_anchored = args_tuple[10] if len(args_tuple) > 10 else False
    allele_specific_cc1 = args_tuple[11] if len(args_tuple) > 11 else False
    n_eco = len(ecotypes)

    # Read the cn_var slice for this block
    # cnvar_slice_path is a pre-computed numpy memmap of (n_eco, n_records_in_block)
    # OR for simplicity we pass record indices and let workers reload cn_var.
    # Here we use the simpler path: pass cn_var slice directly via the tuple.
    # But cn_var is huge; passing per-block via pickle is wasteful. Instead,
    # the parent passes a binary genotype-profile matrix (n_eco × n_var_in_block)
    # already sliced.
    # `cnvar_slice_path` is actually just a small bytes blob we pass directly:
    # the parent passes a numpy array (n_eco, n_var_in_block) uint8.
    cnvar_block = cnvar_slice_path  # parent passes np.ndarray directly
    # cnvar_block is (n_eco, n_var_in_block), dtype uint8, 0/1 entries

    if cnvar_block.shape[1] == 0:
        # Block has no cn_var records → all founders identical at variant level
        # → 1 class containing all 231. No discriminative k-mers possible.
        founder_to_class = np.zeros(n_eco, dtype=np.int32)
        return (block_idx, np.empty(0, dtype=np.uint64),
                csr_matrix((1, 0), dtype=np.uint8), founder_to_class)

    # Dedup founders by their genotype-profile rows; track ALL founders per class
    seen: dict[bytes, int] = {}
    founder_to_class = np.empty(n_eco, dtype=np.int32)
    members_for_class: list[list[int]] = []
    class_genotype_rows: list[np.ndarray] = []  # one cn_var row per class
    for f in range(n_eco):
        key = cnvar_block[f].tobytes()
        cid = seen.get(key)
        if cid is None:
            cid = len(members_for_class)
            seen[key] = cid
            members_for_class.append([f])
            class_genotype_rows.append(cnvar_block[f].copy())
        else:
            members_for_class[cid].append(f)
        founder_to_class[f] = cid
    n_classes = len(members_for_class)

    # Optional Hamming-distance allele grouping: collapse classes whose cn_var
    # genotype profiles differ at <= hamming_threshold variants. Builds a
    # connectivity graph and uses union-find to merge.
    hamming_threshold = args_tuple[12] if len(args_tuple) > 12 else 0
    if hamming_threshold > 0 and n_classes > 1:
        # Pairwise Hamming distances via matrix XOR + sum
        # Each row is uint8 0/1; XOR gives 1 where they differ
        cgr = np.stack(class_genotype_rows)  # (n_classes, n_var)
        # Compute Hamming distances pairwise
        # dist[i,j] = sum(cgr[i] != cgr[j])
        # Use broadcasting; n_classes is typically <100 so OK
        diff = (cgr[:, None, :] != cgr[None, :, :])  # (n_classes, n_classes, n_var) bool
        dist = diff.sum(axis=2)  # (n_classes, n_classes) int
        # Union-find merge
        parent = list(range(n_classes))
        def find(x):
            while parent[x] != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x
        # Iterate upper triangle; merge pairs with dist <= threshold
        iu, ju = np.where(np.triu(dist <= hamming_threshold, k=1))
        for a, b in zip(iu, ju):
            ra, rb = find(int(a)), find(int(b))
            if ra != rb:
                parent[ra] = rb
        # Renumber components 0..M-1
        comp_to_members: dict[int, list[int]] = defaultdict(list)
        for i in range(n_classes):
            comp_to_members[find(i)].append(i)
        comps_sorted = sorted(comp_to_members.values(), key=lambda m: min(m))
        old_to_new = np.empty(n_classes, dtype=np.int32)
        new_members: list[list[int]] = []
        for new_id, mems in enumerate(comps_sorted):
            # Union of original members across merged classes
            union_members: list[int] = []
            for m in mems:
                union_members.extend(members_for_class[m])
                old_to_new[m] = new_id
            new_members.append(union_members)
        members_for_class = new_members
        founder_to_class = old_to_new[founder_to_class].astype(np.int32)
        n_classes = len(members_for_class)

    if n_classes <= 1:
        # All founders in one class → no discrimination possible
        return (block_idx, np.empty(0, dtype=np.uint64),
                csr_matrix((n_classes, 0), dtype=np.uint8), founder_to_class)

    # Pre-compute variant-anchored mask if requested. Positions WITHIN ±(k-1)
    # of any panel-variant record are KEPT; everything else gets N-out so its
    # k-mers fall to canonical_kmer_hashes_np's N-skip filter. This implements
    # PanGenie's "k-mers come from variant-spanning sequence" property at the
    # haplotype-block scale.
    block_len = pos_end - pos_start + 1
    keep_mask: np.ndarray | None = None
    if variant_anchored and variant_positions_in_block is not None and \
            len(variant_positions_in_block) > 0:
        keep_mask = np.zeros(block_len, dtype=bool)
        offsets = np.asarray(variant_positions_in_block, dtype=np.int64) - pos_start
        # Keep positions within ±(k-1) of each variant
        radius = k - 1
        for off in offsets:
            lo = max(int(off) - radius, 0)
            hi = min(int(off) + radius + 1, block_len)
            keep_mask[lo:hi] = True

    # Hash each founder once; per class, take INTERSECTION of class members'
    # k-mer sets. Founders within a class are genotype-identical at all panel
    # variants (by cn_var dedup), so variant-position k-mers must be in ALL
    # members. Per-base FASTA noise (assembly artifacts, IUPAC, etc.) differs
    # across founders → dropped by intersection.
    founder_hashes_cache: dict[int, np.ndarray] = {}
    def _hash_founder(f_idx: int) -> np.ndarray:
        h = founder_hashes_cache.get(f_idx)
        if h is None:
            eco = ecotypes[f_idx]
            fa = _open_fa(os.path.join(fastas_dir, f"{eco}.chr.fa"))
            seq = fa.fetch(chrom, pos_start - 1, pos_end)
            if keep_mask is not None:
                # N-out positions OUTSIDE the variant-anchored mask
                seq_arr = np.frombuffer(seq.encode("ascii"), dtype=np.uint8).copy()
                # Pad keep_mask to seq_arr length if rounding differs
                if seq_arr.size != keep_mask.size:
                    pad = np.zeros(seq_arr.size, dtype=bool)
                    n = min(seq_arr.size, keep_mask.size)
                    pad[:n] = keep_mask[:n]
                    drop = ~pad
                else:
                    drop = ~keep_mask
                seq_arr[drop] = ord("N")
                seq = seq_arr.tobytes().decode("ascii")
            h = canonical_kmer_hashes_np(seq, k=k)
            founder_hashes_cache[f_idx] = h
        return h

    class_hashes: list[np.ndarray] = []
    for c in range(n_classes):
        members = members_for_class[c]
        if len(members) == 1:
            # Singleton class: use the sole founder's k-mer set as-is.
            class_hashes.append(_hash_founder(members[0]))
        else:
            # Intersection of all members' k-mer sets
            inter = _hash_founder(members[0])
            for m in members[1:]:
                if inter.size == 0:
                    break
                inter = np.intersect1d(inter, _hash_founder(m), assume_unique=True)
            class_hashes.append(inter)

    # Build union of k-mer hashes across classes
    if all(h.size == 0 for h in class_hashes):
        return (block_idx, np.empty(0, dtype=np.uint64),
                csr_matrix((n_classes, 0), dtype=np.uint8), founder_to_class)
    union_hashes = np.unique(np.concatenate(class_hashes))

    # Build cn[class, kmer]: sparse rows assembly
    rows: list[np.ndarray] = []
    cols: list[np.ndarray] = []
    for c, h in enumerate(class_hashes):
        if h.size == 0:
            continue
        idx = np.searchsorted(union_hashes, h)
        cols.append(idx)
        rows.append(np.full(idx.size, c, dtype=np.int32))
    rows_arr = np.concatenate(rows) if rows else np.empty(0, dtype=np.int32)
    cols_arr = np.concatenate(cols) if cols else np.empty(0, dtype=np.int64)
    data = np.ones(rows_arr.size, dtype=np.uint8)
    cn = csr_matrix((data, (rows_arr, cols_arr)),
                    shape=(n_classes, union_hashes.size), dtype=np.uint8)

    # Discriminative filter (default: cc in [1, n_classes-1])
    col_sum = np.asarray(cn.sum(axis=0)).flatten()
    if allele_specific_cc1:
        # PanGenie-style: keep only k-mers carried by exactly 1 class
        discrim = col_sum == 1
    else:
        discrim = (col_sum > 0) & (col_sum < n_classes)
    if not discrim.any():
        return (block_idx, np.empty(0, dtype=np.uint64),
                csr_matrix((n_classes, 0), dtype=np.uint8), founder_to_class)
    keep_idx = np.where(discrim)[0]
    kmer_hashes_kept = union_hashes[keep_idx]
    cn_kept = cn[:, keep_idx].tocsr()
    return (block_idx, kmer_hashes_kept, cn_kept, founder_to_class)


def canonical_kmer_hashes_np(seq: str, k: int) -> np.ndarray:
    """Local copy to avoid name shadowing in pickled args."""
    if k > 32:
        raise ValueError(f"k must be <= 32; got {k}")
    if len(seq) < k:
        return np.empty(0, dtype=np.uint64)
    s = np.frombuffer(seq.encode("ascii"), dtype=np.uint8)
    codes = _BASE_LUT[s]
    L = codes.size
    n_kmers = L - k + 1
    if n_kmers <= 0:
        return np.empty(0, dtype=np.uint64)
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
    return np.unique(canonical)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--block-index", required=True,
                    help="hapfire_block_index_<chrom>.npz (for block coords)")
    ap.add_argument("--cn-var", required=True,
                    help="panel cn_var .cn_var.npz (founder × variant)")
    ap.add_argument("--cn-var-meta", required=True,
                    help="panel cn_var .meta.npz (chrom, pos, ref_len, alt_len, founders)")
    ap.add_argument("--fastas-dir", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--k", type=int, default=31)
    ap.add_argument("--chrom-prefix", default="Chr")
    ap.add_argument("--chrom-filter", default=None)
    ap.add_argument("--max-blocks", type=int, default=None)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--variant-anchored", action="store_true",
                    help="N-mask FASTA positions outside ±(k-1) of any cn_var "
                         "record (PanGenie-style: k-mers must span a real variant).")
    ap.add_argument("--allele-specific-cc1", action="store_true",
                    help="filter to k-mers with carrier_count == 1 within block "
                         "(PanGenie-style allele-specific filter).")
    ap.add_argument("--hamming-threshold", type=int, default=0,
                    help="merge classes whose cn_var genotype profiles differ "
                         "at <= this many variants (allele grouping). "
                         "0 = no grouping (default).")
    args = ap.parse_args()

    print(f"Loading block index from {args.block_index}...", flush=True)
    z = np.load(args.block_index, allow_pickle=True)
    ecotypes = np.asarray(z["ecotypes"]).astype(str)
    block_chrom = np.asarray(z["block_chrom"]).astype(str)
    block_pos_start = np.asarray(z["block_pos_start"]).astype(np.int64)
    block_pos_end = np.asarray(z["block_pos_end"]).astype(np.int64)
    n_eco = len(ecotypes)

    print(f"Loading cn_var from {args.cn_var}...", flush=True)
    cn_var = load_npz(args.cn_var)
    cv_meta = np.load(args.cn_var_meta, allow_pickle=True)
    cv_founders = np.asarray(cv_meta["founders"]).astype(str)
    cv_chrom = np.asarray(cv_meta["chrom"]).astype(str)
    cv_pos = np.asarray(cv_meta["pos"]).astype(np.int64)
    print(f"  cn_var shape: {cn_var.shape}, n_records: {cn_var.shape[1]:,}", flush=True)

    # Reorder cn_var rows to match the block-index ecotype order
    if not np.array_equal(cv_founders, ecotypes):
        print("  reordering cn_var rows to match block-index ecotype order...", flush=True)
        rev = {e: i for i, e in enumerate(cv_founders)}
        order = np.array([rev[e] for e in ecotypes])
        cn_var = cn_var[order, :]

    # Build chrom labels (e.g., '1' → 'Chr1')
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
    print(f"Processing {len(keep_idx)} / {n_blocks_total} blocks "
          f"(chrom_filter={args.chrom_filter})", flush=True)

    # Pre-compute per-block cn_var slices for jobs (avoid passing full cn_var
    # to workers which is huge). Slice once in parent, pass small uint8 arrays.
    print("Pre-slicing cn_var per block...", flush=True)
    t = time.time()
    jobs = []
    for b in keep_idx:
        chrom_b = str(block_chrom_full[b])
        chrom_mask = cv_chrom == chrom_b
        if not chrom_mask.any():
            chrom_mask = cv_chrom == chrom_b.replace("Chr", "")
        recs = np.where(chrom_mask & (cv_pos >= int(block_pos_start[b]))
                        & (cv_pos <= int(block_pos_end[b])))[0]
        cnvar_block = cn_var[:, recs].toarray().astype(np.uint8)
        # 1-based positions of variants in this block (for variant-anchored masking)
        var_positions = cv_pos[recs] if recs.size > 0 else np.empty(0, dtype=np.int64)
        jobs.append((
            int(b), chrom_b,
            int(block_pos_start[b]), int(block_pos_end[b]),
            ecotypes, args.fastas_dir, args.k,
            cnvar_block, recs,
            var_positions,
            bool(args.variant_anchored),
            bool(args.allele_specific_cc1),
            int(args.hamming_threshold),
        ))
    print(f"  pre-sliced {len(jobs)} blocks [{time.time()-t:.0f}s]", flush=True)
    sizes = [j[7].shape[1] for j in jobs]
    print(f"  per-block n_var: median={int(np.median(sizes))}, "
          f"min={min(sizes)}, max={max(sizes)}", flush=True)

    # Worker pool
    t0 = time.time()
    results: dict[int, tuple] = {}
    if args.threads > 1:
        with Pool(args.threads) as pool:
            for i, res in enumerate(pool.imap_unordered(
                process_block_cnvar, jobs, chunksize=4)):
                bidx, kmers, cn, f2c = res
                results[bidx] = (kmers, cn, f2c)
                if (i + 1) % 200 == 0:
                    print(f"  {i+1}/{len(jobs)} blocks [{time.time()-t0:.0f}s]",
                          flush=True)
    else:
        for i, j in enumerate(jobs):
            bidx, kmers, cn, f2c = process_block_cnvar(j)
            results[bidx] = (kmers, cn, f2c)
            if (i + 1) % 200 == 0:
                print(f"  {i+1}/{len(jobs)} blocks [{time.time()-t0:.0f}s]",
                      flush=True)
    print(f"per-block hashing done [{time.time()-t0:.0f}s]", flush=True)

    # Cross-block global dedup: drop k-mers in >1 block's discriminative set
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
    n_global = uniq.size if all_hashes else 0
    n_unique_global = int(unique_set_view.size)
    n_dropped = n_global - n_unique_global
    print(f"  global k-mers: {n_global:,}; unique to one block: {n_unique_global:,}; "
          f"dropped: {n_dropped:,} ({n_dropped/max(1,n_global)*100:.1f}%) "
          f"[{time.time()-t:.0f}s]", flush=True)

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
        kept_hashes = kmer_hashes[keep_local_idx] if keep_local_idx.size else \
            np.empty(0, dtype=np.uint64)
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
    print(f"  decoding {kmer_hashes_concat.size:,} k-mer hashes to strings...", flush=True)
    kmer_strings_all = [hash_to_kmer(int(h), args.k) for h in kmer_hashes_concat]

    # Concat per-block sparse cn matrices into one big CSR
    print(f"  concatenating cn rows for {len(cn_local_blocks)} blocks...", flush=True)
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
