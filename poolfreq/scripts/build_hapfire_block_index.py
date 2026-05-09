"""Build a per-block ecotype-to-haplotype-index mapping for hapFIRE blocks.

Inputs:
  --vcf      panel VCF (greneNet_final_v1.1.recode.vcf, ~3.2M SNPs, 231 ecotypes)
  --blocks   genome-wide partition file (chrom, start_snp_idx, end_snp_idx)
  --out      output .npz with per-block lookup

For each block, we want to know "which ecotype carries which unique haplotype",
matching hapFIRE's internal indexing (`drop_duplicates(keep='first')` over
[hap_d1; hap_d2]). This is needed to convert hapFIRE's per-block per-haplotype
frequency outputs into per-block per-ecotype frequencies, which we then
project through cn_var to get per-record AFs (including SVs).

Saved npz:
  ecotypes:        N_eco-vec of ecotype IDs (str)
  block_chrom:     N_block-vec of chrom strings ('1','2',...) — matches partition file
  block_snp_start: N_block-vec of starting SNP idx (within-chrom)
  block_snp_end:   N_block-vec of ending SNP idx (within-chrom)
  block_pos_start: N_block-vec of starting bp position
  block_pos_end:   N_block-vec of ending bp position
  block_n_snps:    N_block-vec of SNP count per block
  block_n_uniq:    N_block-vec of unique-haplotype count per block
  hap_idx_d1:      N_block × N_eco int matrix; ecotype's d1 haplotype index in this block
  hap_idx_d2:      N_block × N_eco int matrix; ecotype's d2 haplotype index in this block
"""
from __future__ import annotations
import argparse, sys, time
from pathlib import Path
import numpy as np
import pandas as pd


def parse_vcf(vcf_path):
    """Parse phased SNP VCF → (chroms_per_snp, positions_per_snp, ecotypes,
    hap_d1, hap_d2). hap_d1, hap_d2 are bool matrices [n_eco × n_snp]."""
    print(f'Parsing {vcf_path} ...', flush=True)
    t = time.time()
    eco_names = None
    chroms = []
    positions = []
    d1_chunks = []
    d2_chunks = []
    chunk_d1 = []
    chunk_d2 = []
    n = 0
    chunk_size = 100000
    with open(vcf_path) as f:
        for line in f:
            if line.startswith('##'):
                continue
            if line.startswith('#CHROM'):
                eco_names = line.rstrip('\n').split('\t')[9:]
                continue
            parts = line.rstrip('\n').split('\t')
            chrom, pos = parts[0], int(parts[1])
            gts = parts[9:]
            d1 = np.empty(len(gts), dtype=np.uint8)
            d2 = np.empty(len(gts), dtype=np.uint8)
            for i, g in enumerate(gts):
                # accept '0|0', '0|1', '1|0', '1|1' etc.
                if g[0] == '.':
                    d1[i] = 0
                else:
                    d1[i] = int(g[0])
                if len(g) >= 3 and g[2] != '.':
                    d2[i] = int(g[2])
                else:
                    d2[i] = 0
            chroms.append(chrom)
            positions.append(pos)
            chunk_d1.append(d1)
            chunk_d2.append(d2)
            n += 1
            if n % chunk_size == 0:
                d1_chunks.append(np.stack(chunk_d1, axis=1))
                d2_chunks.append(np.stack(chunk_d2, axis=1))
                chunk_d1 = []
                chunk_d2 = []
                if n % 500000 == 0:
                    print(f'  {n:,} SNPs parsed [{time.time()-t:.0f}s]', flush=True)
    if chunk_d1:
        d1_chunks.append(np.stack(chunk_d1, axis=1))
        d2_chunks.append(np.stack(chunk_d2, axis=1))
    hap_d1 = np.concatenate(d1_chunks, axis=1)  # n_eco × n_snp
    hap_d2 = np.concatenate(d2_chunks, axis=1)
    print(f'  total: n_snps={n:,}, n_ecotypes={len(eco_names)}, '
          f'shape={hap_d1.shape} [{time.time()-t:.0f}s]', flush=True)
    return (np.asarray(chroms), np.asarray(positions, dtype=np.int64),
            np.asarray(eco_names), hap_d1, hap_d2)


def haplotype_index_block(block_d1, block_d2):
    """For a block's [n_eco, block_len] d1 + d2 matrices, return per-ecotype d1
    and d2 unique-haplotype indices following hapFIRE's
    `pd.concat([d1, d2]).drop_duplicates(keep='first')` ordering.

    Matches hapFIRE's '_l' indexing in unique_haplotype_frequency outputs.
    """
    n_eco = block_d1.shape[0]
    # Concatenate haplotype strings: rows 0..n-1 are d1, n..2n-1 are d2.
    haps = np.concatenate([block_d1, block_d2], axis=0)  # 2n × block_len
    # Hash via tobytes() per row for fast uniqueness
    seen = {}
    out_d1 = np.empty(n_eco, dtype=np.int32)
    out_d2 = np.empty(n_eco, dtype=np.int32)
    next_id = 0
    for i in range(haps.shape[0]):
        key = haps[i].tobytes()
        if key not in seen:
            seen[key] = next_id
            next_id += 1
        idx = seen[key]
        if i < n_eco:
            out_d1[i] = idx
        else:
            out_d2[i - n_eco] = idx
    return out_d1, out_d2, next_id


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--vcf', default='/carnegie/nobackup/scratch/xwu/GrENE_net/hapFIRE_frequencies/greneNet_final_v1.1.recode.vcf')
    ap.add_argument('--blocks', default='/carnegie/nobackup/scratch/xwu/GrENE_net/hapFIRE_frequencies/greneNet_final_v1.1_density0.5_genomewide_partition.txt')
    ap.add_argument('--out', default='/carnegie/nobackup/scratch/tbellagio/hapfire_sv/poolfreq/data/hapfire_block_index.npz')
    args = ap.parse_args()

    chroms, positions, ecotypes, hap_d1, hap_d2 = parse_vcf(args.vcf)
    print(f'Reading blocks from {args.blocks}')
    blocks_df = pd.read_csv(args.blocks, sep='\t', header=None,
                            names=['chrom', 'start_idx', 'end_idx'])
    blocks_df['chrom'] = blocks_df['chrom'].astype(str)
    print(f'  n_blocks={len(blocks_df):,}')

    # The partition uses within-chrom SNP indices. Build chrom-offset table.
    snp_chrom_arr = np.asarray(chroms).astype(str)
    chrom_to_first_idx = {}
    for i, c in enumerate(snp_chrom_arr):
        if c not in chrom_to_first_idx:
            chrom_to_first_idx[c] = i

    n_blocks = len(blocks_df)
    n_eco = len(ecotypes)
    block_chrom = blocks_df['chrom'].to_numpy()
    block_start = blocks_df['start_idx'].to_numpy(dtype=np.int64)
    block_end = blocks_df['end_idx'].to_numpy(dtype=np.int64)

    block_pos_start = np.zeros(n_blocks, dtype=np.int64)
    block_pos_end = np.zeros(n_blocks, dtype=np.int64)
    block_n_snps = np.zeros(n_blocks, dtype=np.int32)
    block_n_uniq = np.zeros(n_blocks, dtype=np.int32)
    hap_idx_d1 = np.zeros((n_blocks, n_eco), dtype=np.int32)
    hap_idx_d2 = np.zeros((n_blocks, n_eco), dtype=np.int32)

    print('Computing per-block haplotype indexing...', flush=True)
    t = time.time()
    for b in range(n_blocks):
        c = block_chrom[b]
        # Convert within-chrom SNP indices to global SNP indices
        offset = chrom_to_first_idx.get(c, 0)
        gs = offset + block_start[b]
        ge = offset + block_end[b]
        if ge < gs:
            ge = gs
        block_d1 = hap_d1[:, gs:ge + 1]
        block_d2 = hap_d2[:, gs:ge + 1]
        if block_d1.shape[1] == 0:
            continue
        block_pos_start[b] = positions[gs]
        block_pos_end[b] = positions[ge]
        block_n_snps[b] = block_d1.shape[1]

        d1_idx, d2_idx, n_uniq = haplotype_index_block(block_d1, block_d2)
        hap_idx_d1[b] = d1_idx
        hap_idx_d2[b] = d2_idx
        block_n_uniq[b] = n_uniq
        if (b + 1) % 1000 == 0:
            print(f'  {b+1}/{n_blocks} blocks [{time.time()-t:.0f}s]', flush=True)
    print(f'  done [{time.time()-t:.0f}s]', flush=True)

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    np.savez(
        args.out,
        ecotypes=ecotypes,
        block_chrom=block_chrom,
        block_snp_start=block_start,
        block_snp_end=block_end,
        block_pos_start=block_pos_start,
        block_pos_end=block_pos_end,
        block_n_snps=block_n_snps,
        block_n_uniq=block_n_uniq,
        hap_idx_d1=hap_idx_d1,
        hap_idx_d2=hap_idx_d2,
    )
    print(f'wrote {args.out}', flush=True)
    print(f'  block_n_snps: min={block_n_snps.min()}, median={int(np.median(block_n_snps))}, max={block_n_snps.max()}, mean={block_n_snps.mean():.1f}')
    print(f'  block_n_uniq: min={block_n_uniq.min()}, median={int(np.median(block_n_uniq))}, max={block_n_uniq.max()}, mean={block_n_uniq.mean():.1f}')


if __name__ == '__main__':
    main()
