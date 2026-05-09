"""Project per-block hapFIRE ecotype freqs through cn_var → per-record AFs.

For each VCF record at (chrom, pos):
  1. find the LD block whose [pos_start, pos_end] contains pos
  2. h_block = h_per_block[block_idx, :]  (231-vec)
  3. pred_af[r] = h_block · cn_var[:, r]

Writes a TSV with the same schema as cactus_em outputs so the comparison is
direct: chrom, pos, ref_len, alt_len, alt_freq.

Records outside any block are written with alt_freq = NaN.

Note: hapFIRE's blocks use chrom IDs '1','2',...; cn_var uses 'Chr1','Chr2',...
We translate.
"""
from __future__ import annotations
import argparse, time
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.sparse import load_npz


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--perblock-h', required=True,
                    help='per-sample npz from derive_hapfire_perblock_h.py')
    ap.add_argument('--cn-var', default='/carnegie/nobackup/scratch/tbellagio/hapfire_sv/poolfreq/data/cn_var_231_v2.cn_var.npz')
    ap.add_argument('--cn-var-meta', default='/carnegie/nobackup/scratch/tbellagio/hapfire_sv/poolfreq/data/cn_var_231_v2.meta.npz')
    ap.add_argument('--out', required=True)
    args = ap.parse_args()

    print('Loading per-block h...', flush=True)
    pb = np.load(args.perblock_h, allow_pickle=True)
    h = np.asarray(pb['h_per_block']).astype(np.float32)            # (n_blocks, n_eco)
    block_chrom_raw = np.asarray(pb['block_chrom']).astype(str)
    block_pos_start = np.asarray(pb['block_pos_start']).astype(np.int64)
    block_pos_end = np.asarray(pb['block_pos_end']).astype(np.int64)
    block_eco_order = [str(e) for e in pb['ecotypes']]
    n_blocks, n_eco = h.shape

    # Translate hapFIRE chrom IDs ('1') → cn_var IDs ('Chr1')
    block_chrom = np.array([f'Chr{c}' if not c.startswith('Chr') else c
                            for c in block_chrom_raw])

    print('Loading cn_var + meta...', flush=True)
    cn = load_npz(args.cn_var).tocsc()       # 231 × N records
    meta = np.load(args.cn_var_meta, allow_pickle=True)
    cn_eco_order = [str(e) for e in meta['founders']]
    rec_chrom = np.asarray(meta['chrom']).astype(str)
    rec_pos = np.asarray(meta['pos']).astype(np.int64)
    rec_ref = np.asarray(meta['ref_len']).astype(np.int32)
    rec_alt = np.asarray(meta['alt_len']).astype(np.int32)
    N = len(rec_chrom)
    print(f'  cn_var: {cn.shape}, records={N:,}', flush=True)

    # Reorder h's columns to match cn_var founder order
    if cn_eco_order != block_eco_order:
        idx_in_h = {e: i for i, e in enumerate(block_eco_order)}
        col_order = np.array([idx_in_h.get(e, -1) for e in cn_eco_order])
        if (col_order < 0).any():
            missing = [e for e, i in zip(cn_eco_order, col_order) if i < 0]
            raise ValueError(f'cn_var has ecotypes not in hapFIRE panel: {missing[:5]} ...')
        h = h[:, col_order]
        print(f'  reordered h columns to match cn_var founder order', flush=True)

    # Per-record block lookup: for each record, which block index?
    print('Assigning records to blocks...', flush=True)
    t = time.time()
    rec_block = np.full(N, -1, dtype=np.int32)
    for chrom in np.unique(block_chrom):
        bm = block_chrom == chrom
        b_idx = np.flatnonzero(bm)
        # Sort blocks on this chrom by start position (should already be sorted)
        starts = block_pos_start[b_idx]
        ends = block_pos_end[b_idx]
        order = np.argsort(starts)
        b_idx = b_idx[order]
        starts = starts[order]
        ends = ends[order]

        rm = rec_chrom == chrom
        if not rm.any():
            continue
        rec_pos_c = rec_pos[rm]
        # binary-search starts: position of largest start <= pos
        ii = np.searchsorted(starts, rec_pos_c, side='right') - 1
        ii = np.clip(ii, 0, len(starts) - 1)
        ok = (rec_pos_c >= starts[ii]) & (rec_pos_c <= ends[ii])
        block_for_rec = np.where(ok, b_idx[ii], -1)
        rec_block[rm] = block_for_rec
    n_assigned = int((rec_block >= 0).sum())
    print(f'  {n_assigned:,}/{N:,} records assigned to blocks ({100*n_assigned/N:.2f}%) [{time.time()-t:.0f}s]', flush=True)

    # Project: pred_af[r] = h_per_block[rec_block[r], :] @ cn[:, r]
    # Vectorize per-block: for each block, gather its records and do a matvec
    print('Projecting...', flush=True)
    t = time.time()
    pred_af = np.full(N, np.nan, dtype=np.float32)
    unique_blocks = np.unique(rec_block[rec_block >= 0])
    for b in unique_blocks:
        recs_in_b = np.flatnonzero(rec_block == b)
        if len(recs_in_b) == 0:
            continue
        h_b = h[b]       # (n_eco,)
        cv_sub = cn[:, recs_in_b]            # (n_eco, n_recs_in_b) sparse
        af = (h_b @ cv_sub).A1 if hasattr((h_b @ cv_sub), 'A1') \
             else np.asarray(h_b @ cv_sub).flatten()
        pred_af[recs_in_b] = af.astype(np.float32)
    print(f'  projection done [{time.time()-t:.0f}s]', flush=True)

    # Write TSV
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, 'w') as f:
        f.write('chrom\tpos\tref_len\talt_len\talt_freq\n')
        for i in range(N):
            af = pred_af[i]
            af_str = f'{af:.5f}' if np.isfinite(af) else 'NaN'
            f.write(f'{rec_chrom[i]}\t{rec_pos[i]}\t{rec_ref[i]}\t{rec_alt[i]}\t{af_str}\n')
    print(f'wrote {args.out}', flush=True)


if __name__ == '__main__':
    main()
