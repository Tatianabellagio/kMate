"""Derive per-block per-ecotype frequency matrix from hapFIRE outputs.

Inputs:
  --block-index       hapfire_block_index.npz from build_hapfire_block_index.py
  --uniq-haplo-tsv    sample's <name>_unique_haplotype_frequency.txt (per-block per-haplotype freqs)
  --out               output .npz with per-sample h_per_block (n_blocks × n_ecotypes)

For each (block, ecotype) we compute:
    h_per_block[block, eco] = (pool_freq[block, hap_d1[block, eco]] + pool_freq[block, hap_d2[block, eco]]) / 2

For mostly-selfing A. thaliana, hap_d1 == hap_d2 for almost all ecotypes, so this
collapses to h[block, eco] = pool_freq[block, hap_idx_of_eco_in_block].

Sanity check: sum over ecotypes of h_per_block[b] should equal sum over
unique haplotypes of pool_freq[b] × (count_of_eco_carrying_each_uniq_in_block)
divided by 2 × n_eco. For a normalized pool this is close to 1.
"""
from __future__ import annotations
import argparse, time
from pathlib import Path
import numpy as np
import pandas as pd


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--block-index', required=True)
    ap.add_argument('--uniq-haplo-tsv', required=True,
                    help='Sample\'s _unique_haplotype_frequency.txt')
    ap.add_argument('--out', required=True)
    args = ap.parse_args()

    print(f'Loading block index from {args.block_index}', flush=True)
    idx = np.load(args.block_index, allow_pickle=True)
    block_chrom = np.asarray(idx['block_chrom']).astype(str)
    block_pos_start = np.asarray(idx['block_pos_start']).astype(np.int64)
    block_pos_end = np.asarray(idx['block_pos_end']).astype(np.int64)
    block_n_uniq = np.asarray(idx['block_n_uniq']).astype(np.int32)
    hap_idx_d1 = np.asarray(idx['hap_idx_d1']).astype(np.int32)
    hap_idx_d2 = np.asarray(idx['hap_idx_d2']).astype(np.int32)
    ecotypes = np.asarray(idx['ecotypes']).astype(str)
    n_blocks = len(block_chrom)
    n_eco = len(ecotypes)
    print(f'  n_blocks={n_blocks:,}  n_ecotypes={n_eco}', flush=True)

    # Build "block-name → block_idx" lookup keyed by 'chrom@start-end'
    block_key = np.array([f'{c}@{s}-{e}' for c, s, e in
                          zip(block_chrom, block_pos_start, block_pos_end)])
    block_key_to_idx = {k: i for i, k in enumerate(block_key)}

    print(f'Reading {args.uniq_haplo_tsv}', flush=True)
    t = time.time()
    pool_freq_per_block = [None] * n_blocks  # list of arrays per block
    n_rows = 0
    with open(args.uniq_haplo_tsv) as f:
        for line in f:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 2:
                continue
            name, pool_freq_str = parts[0], parts[1]
            # name format: <chrom>@<start>-<end>_<hap_idx>
            try:
                block_part, hap_str = name.rsplit('_', 1)
                hap_i = int(hap_str)
            except ValueError:
                continue
            b_idx = block_key_to_idx.get(block_part)
            if b_idx is None:
                continue
            if pool_freq_per_block[b_idx] is None:
                pool_freq_per_block[b_idx] = np.zeros(block_n_uniq[b_idx], dtype=np.float32)
            try:
                pf = float(pool_freq_str)
            except ValueError:
                pf = 0.0
            if 0 <= hap_i < block_n_uniq[b_idx]:
                pool_freq_per_block[b_idx][hap_i] = pf
            n_rows += 1
    print(f'  parsed {n_rows:,} rows [{time.time()-t:.0f}s]', flush=True)

    # Project: pool_freq[u] is normalized to sum to 1 per block; multiple
    # ecotypes may share unique haplotype u in the block. Distribute pool_freq
    # evenly across the ecotypes carrying u:
    #   h_per_block[b, eco] = 0.5 × ( pf[d1[eco]] / count_d1[d1[eco]]
    #                                 + pf[d2[eco]] / count_d2[d2[eco]] )
    # For selfing (d1==d2, count_d1==count_d2): collapses to
    #   h_per_block[b, eco] = pf[d1[eco]] / count_d1[d1[eco]]
    # so Σ_eco h_per_block[b, eco] = Σ_u pf[u] = 1 per block.
    print(f'Projecting per-block per-eco freqs (with even distribution across shared haplotypes)...', flush=True)
    t = time.time()
    h_per_block = np.zeros((n_blocks, n_eco), dtype=np.float32)
    for b in range(n_blocks):
        pf = pool_freq_per_block[b]
        if pf is None:
            continue
        d1 = hap_idx_d1[b]
        d2 = hap_idx_d2[b]
        # Per-block ecotype-counts per unique haplotype in d1 and d2
        n_uniq = block_n_uniq[b]
        cnt_d1 = np.bincount(d1, minlength=n_uniq).astype(np.float32)
        cnt_d2 = np.bincount(d2, minlength=n_uniq).astype(np.float32)
        # Avoid divide-by-zero for unique haplotypes carried by 0 ecotypes
        cnt_d1_safe = np.where(cnt_d1 > 0, cnt_d1, 1.0)
        cnt_d2_safe = np.where(cnt_d2 > 0, cnt_d2, 1.0)
        contrib_d1 = pf[d1] / cnt_d1_safe[d1]
        contrib_d2 = pf[d2] / cnt_d2_safe[d2]
        h_per_block[b] = 0.5 * (contrib_d1 + contrib_d2)
        if (b + 1) % 5000 == 0:
            print(f'  {b+1}/{n_blocks} blocks [{time.time()-t:.0f}s]', flush=True)

    # Per-block sanity: row sums should be ~1 (pool-normalized)
    sums = h_per_block.sum(axis=1)
    nonzero = sums > 0
    print(f'  per-block sum-of-eco-freqs (target=1.0): '
          f'median={np.median(sums[nonzero]):.4f}, '
          f'mean={sums[nonzero].mean():.4f}, '
          f'std={sums[nonzero].std():.4f}, '
          f'n_zero_blocks={int((~nonzero).sum())}', flush=True)

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    np.savez(
        args.out,
        ecotypes=ecotypes,
        block_chrom=block_chrom,
        block_pos_start=block_pos_start,
        block_pos_end=block_pos_end,
        block_n_uniq=block_n_uniq,
        h_per_block=h_per_block,
    )
    print(f'wrote {args.out}', flush=True)


if __name__ == '__main__':
    main()
