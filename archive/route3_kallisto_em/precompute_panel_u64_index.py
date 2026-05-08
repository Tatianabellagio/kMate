"""Precompute uint64 k-mer index for cn_full panels — saves a sorted uint64
array + sort_idx so the kallisto-EM driver doesn't have to load 2 GB of
panel k-mer strings at runtime.

Output (per chrom):
  <prefix>_<chrom>_kmers_u64.npz
    - sorted_u64: (K,) uint64, ACGT 2-bit packed canonical k-mers, sorted ascending
    - sort_idx:   (K,) int64, original column index in cn_full[:, k]
    - n_invalid:  int — count of k-mers that contained N (encoded as 0xFFFF... sentinel; kept in the array)

Usage:
    python precompute_panel_u64_index.py \
        --cn-kmer-prefix poolfreq/data/cn_full_231_v2/cn \
        --chroms Chr1 Chr2 Chr3 Chr4 Chr5
"""
from __future__ import annotations
import argparse
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                 "..", "src"))
from per_sample_kallisto_em import _BASE2BIT


def encode_kmer_strings_to_u64(kmer_strings: np.ndarray, k: int = 31) -> np.ndarray:
    """Vectorized 2-bit encoding of an array of bytes/string k-mers to uint64.

    K-mers containing any non-ACGT base get the sentinel 0xFFFFFFFFFFFFFFFF.
    Memory-conscious: processes in chunks of 1M to bound peak memory.
    """
    N = len(kmer_strings)
    out = np.zeros(N, dtype=np.uint64)
    sentinel = np.uint64(0xFFFFFFFFFFFFFFFF)
    chunk = 1_000_000
    for start in range(0, N, chunk):
        end = min(start + chunk, N)
        sub = kmer_strings[start:end]
        # Convert each string to bytes if needed — vectorized via bytes view
        if N > 0 and isinstance(sub[0], str):
            sub = np.array([s.encode("ascii") for s in sub], dtype=f"|S{k}")
        # Build (chunk_size, k) uint8 array of base codes
        sub_bytes = np.frombuffer(sub.tobytes(), dtype=np.uint8).reshape(-1, k)
        codes = _BASE2BIT[sub_bytes]              # (chunk, k) uint8
        invalid = (codes == 0xFF).any(axis=1)     # (chunk,)
        # Pack 2-bit codes into one uint64 per row
        # Use shift-OR loop over k positions
        vals = np.zeros(end - start, dtype=np.uint64)
        for j in range(k):
            vals = (vals << np.uint64(2)) | codes[:, j].astype(np.uint64)
        vals[invalid] = sentinel
        out[start:end] = vals
    return out


def encode_canonical_u64(u64_fwd: np.ndarray, k: int = 31) -> np.ndarray:
    """Given forward 2-bit-packed k-mer uint64s, return the canonical (min of fwd, revcomp) form."""
    sentinel = np.uint64(0xFFFFFFFFFFFFFFFF)
    valid = u64_fwd != sentinel
    rev = np.zeros_like(u64_fwd)
    # Reverse-complement at uint64 level: complement bits and reverse 2-bit positions
    # Complement: XOR with mask of 0b11 in each 2-bit slot of length 2*k
    full_mask = (np.uint64(1) << np.uint64(2 * k)) - np.uint64(1)
    cmp_mask = np.uint64(0)
    for j in range(k):
        cmp_mask |= np.uint64(0b11) << np.uint64(2 * j)
    cmp_mask &= full_mask
    comp = u64_fwd ^ cmp_mask
    # Reverse 2-bit positions: build out incrementally
    rev = np.zeros_like(u64_fwd)
    for j in range(k):
        bits = (comp >> np.uint64(2 * j)) & np.uint64(0b11)
        rev |= bits << np.uint64(2 * (k - 1 - j))
    canon = np.where(valid, np.minimum(u64_fwd, rev), sentinel)
    return canon


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cn-kmer-prefix", required=True)
    ap.add_argument("--chroms", nargs="+", required=True)
    ap.add_argument("--k", type=int, default=31)
    args = ap.parse_args()

    for c in args.chroms:
        meta_path = args.cn_kmer_prefix + f"_{c}.meta.npz"
        out_path = args.cn_kmer_prefix + f"_{c}_kmers_u64.npz"
        if not os.path.exists(meta_path):
            print(f"  WARN: {meta_path} missing; skip", flush=True)
            continue
        if os.path.exists(out_path):
            print(f"  {out_path} exists; skipping", flush=True)
            continue
        print(f"[{c}] loading kmer_index from {meta_path}...", flush=True)
        t = time.time()
        meta = np.load(meta_path, allow_pickle=True)
        kmer_index = meta["kmer_index"]
        K = len(kmer_index)
        print(f"  K={K:,} k-mers; loaded in {time.time()-t:.0f}s", flush=True)

        print(f"[{c}] encoding to uint64 (forward)...", flush=True)
        t = time.time()
        fwd_u64 = encode_kmer_strings_to_u64(kmer_index, k=args.k)
        print(f"  encoded forward in {time.time()-t:.0f}s; "
              f"{fwd_u64.nbytes/1e9:.2f} GB", flush=True)
        del kmer_index, meta

        # The PanGenie kmer_index is already canonical (alphabetically minimum
        # of fwd/revcomp at index time). We confirm this by checking that
        # canonicalize(fwd_u64) == fwd_u64 for valid (non-N) k-mers.
        # If not all match, we re-canonicalize.
        print(f"[{c}] canonicalizing (sanity check)...", flush=True)
        t = time.time()
        # Sample-based check — encode small subset and ensure canonical
        canon_u64 = encode_canonical_u64(fwd_u64, k=args.k)
        same = canon_u64 == fwd_u64
        sentinel = np.uint64(0xFFFFFFFFFFFFFFFF)
        valid = fwd_u64 != sentinel
        n_canonical = int((same & valid).sum())
        n_valid = int(valid.sum())
        if n_canonical != n_valid:
            print(f"  panel kmers NOT all canonical ({n_canonical}/{n_valid}); "
                  f"using canon_u64 from re-canonicalization", flush=True)
            fwd_u64 = canon_u64
        else:
            print(f"  panel is canonical ({n_canonical}/{n_valid} match) "
                  f"[{time.time()-t:.0f}s]", flush=True)
            del canon_u64

        print(f"[{c}] sorting...", flush=True)
        t = time.time()
        sort_idx = np.argsort(fwd_u64).astype(np.int64)
        sorted_u64 = fwd_u64[sort_idx]
        del fwd_u64
        print(f"  sorted in {time.time()-t:.0f}s", flush=True)

        n_invalid = int((sorted_u64 == sentinel).sum())
        print(f"[{c}] saving to {out_path} (n_invalid={n_invalid})...",
              flush=True)
        np.savez_compressed(out_path,
                            sorted_u64=sorted_u64,
                            sort_idx=sort_idx,
                            n_invalid=np.int64(n_invalid),
                            k=np.int32(args.k))
        print(f"  saved.", flush=True)


if __name__ == "__main__":
    main()
