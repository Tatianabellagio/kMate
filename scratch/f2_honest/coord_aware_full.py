#!/usr/bin/env python3
"""
Coordinate-aware coalescence — applied to ALL Chr1 SNP cn_var records.

For each cn_var record r (the "source"):
  For each offset in [0, len(REF[r])):
    target_pos = pos[r] + offset
    if offset < len(ALT[r]):
      target_base = ALT[r][offset]
      For each SNP cn_var col at (chrom=Chr1, pos=target_pos, ALT[0]=target_base):
        merge cn_var[:, r] carriers into that SNP's extra-carrier mask

After processing all records: for each SNP, AC_new = | (own carriers) OR (extras) |.
This is the principled fix for the bcftools-norm position-shifting bug.
"""
import numpy as np
import pandas as pd
import scipy.sparse as sp
from pathlib import Path
import time

ROOT = Path('/carnegie/nobackup/scratch/tbellagio/hapfire_sv')

print('[1] Load cn_var + meta')
cv = sp.load_npz(ROOT / 'poolfreq/data/cn_var_231_v3qc_v3.cn_var.npz').tocsc()
m  = np.load(ROOT / 'poolfreq/data/cn_var_231_v3qc_v3.meta.npz', allow_pickle=True)
chrom_arr   = m['chrom']; pos_arr = m['pos']
ref_arr     = m['ref'];   alt_arr = m['alt']
ref_len_arr = m['ref_len']; alt_len_arr = m['alt_len']
print(f'    cn_var: {cv.shape}, {cv.nnz:,} nnz')

# Chr1 records (all variant types) — these are the "sources"
chr1_mask = (chrom_arr == 'Chr1')
chr1_cols = np.where(chr1_mask)[0]
print(f'    Chr1 records (all types): {len(chr1_cols):,}')

# Chr1 SNP records — these are the "targets"
snp_mask = chr1_mask & (ref_len_arr == 1) & (alt_len_arr == 1)
snp_cols = np.where(snp_mask)[0]
print(f'    Chr1 SNP records:         {len(snp_cols):,}')

# Build snp index: (pos, alt_base) -> list of snp col indices (usually 1, occasionally 2+ for multi-allelic)
print('[2] Build SNP index (pos, alt_base) -> col_idx')
snp_index = {}
snp_pos_array = pos_arr[snp_cols]
snp_alt_array = np.array([str(a) for a in alt_arr[snp_cols]])
for i, col in enumerate(snp_cols):
    key = (int(snp_pos_array[i]), snp_alt_array[i])
    snp_index.setdefault(key, []).append(i)   # i is local idx into snp_cols
print(f'    unique (pos, alt_base) keys: {len(snp_index):,}')

# Sorted SNP positions for binary search
snp_pos_sorted_uniq = np.array(sorted(set(int(p) for p in snp_pos_array)))
print(f'    unique SNP positions: {len(snp_pos_sorted_uniq):,}')

# Pre-extract record arrays for Chr1
print('[3] Pre-extract Chr1 record arrays')
src_pos = pos_arr[chr1_cols].astype(np.int64)
src_ref_len = ref_len_arr[chr1_cols].astype(np.int32)
src_alt_str = [str(a) for a in alt_arr[chr1_cols]]   # Python list to avoid numpy str alloc
print(f'    max ref_len on Chr1: {src_ref_len.max()}')
print(f'    median ref_len     : {int(np.median(src_ref_len))}')

# Allocate extras mask: (n_snps, 231) int8
n_snps = len(snp_cols)
extras_mask = np.zeros((n_snps, 231), dtype=np.int8)
print(f'    extras_mask shape: {extras_mask.shape}  size={extras_mask.nbytes/1e6:.1f} MB')

# Iterate Chr1 records, find spanning SNPs, merge carriers
print('[4] Iterate Chr1 records, find spanning SNPs, OR-merge carriers')
t0 = time.time()
n_processed = 0
n_carrier_merges = 0
for ix, src_col in enumerate(chr1_cols):
    p = int(src_pos[ix]); rl = int(src_ref_len[ix]); alt = src_alt_str[ix]
    # Find SNP positions in [p, p + rl)
    lo = np.searchsorted(snp_pos_sorted_uniq, p, side='left')
    hi = np.searchsorted(snp_pos_sorted_uniq, p + rl, side='left')
    if hi <= lo:
        n_processed += 1
        continue
    # For each spanning SNP position, compute offset & check ALT-base match
    src_carriers = None  # lazy-load
    for target_pos in snp_pos_sorted_uniq[lo:hi]:
        offset = int(target_pos) - p
        if offset >= len(alt):
            continue
        target_alt_base = alt[offset]
        key = (int(target_pos), target_alt_base)
        snp_locals = snp_index.get(key)
        if snp_locals is None:
            continue
        # Lazy-load source carrier vec
        if src_carriers is None:
            src_carriers = np.asarray(cv[:, src_col].todense()).ravel().astype(np.int8)
        # OR-merge into each matching SNP target's extras_mask
        for sl in snp_locals:
            if sl == ix:
                # source IS the target — skip (own carriers are tracked separately)
                # but actually src_col == snp_cols[sl] only if record is itself a SNP at this pos
                # and target_pos==p, target_alt_base==alt[0]. Don't OR own into own's extras.
                pass
            # Actually we should also skip when src_col == snp_cols[sl]
            if chr1_cols[ix] == snp_cols[sl]:
                continue
            extras_mask[sl] = np.maximum(extras_mask[sl], src_carriers)
            n_carrier_merges += 1
    n_processed += 1
    if (n_processed) % 200000 == 0:
        print(f'    processed {n_processed:,} / {len(chr1_cols):,} '
              f'records ({time.time()-t0:.0f}s)  carrier_merges={n_carrier_merges:,}')
print(f'    DONE in {time.time()-t0:.0f}s. total carrier_merges = {n_carrier_merges:,}')

# Compute new AC = | own carriers OR extras | per SNP
print('[5] Compute new AC per SNP and re-project AF')
own_AC = np.zeros(n_snps, dtype=np.int32)
new_AC = np.zeros(n_snps, dtype=np.int32)
extras_only = np.zeros(n_snps, dtype=np.int32)
for i, snp_col in enumerate(snp_cols):
    own = np.asarray(cv[:, snp_col].todense()).ravel().astype(np.int8)
    own_AC[i] = int(own.sum())
    combined = np.maximum(own, extras_mask[i])
    new_AC[i] = int(combined.sum())
    extras_only[i] = int(((extras_mask[i] == 1) & (own == 0)).sum())

print(f'    n SNPs with extras > 0: {(extras_only>0).sum():,} / {n_snps:,}  '
      f'({100*(extras_only>0).mean():.1f}%)')

# Build summary DataFrame for joining with F1 table
print('[6] Build summary, join with F1 table, evaluate per cell')
snp_meta = pd.DataFrame({
    'chrom': chrom_arr[snp_cols],
    'pos':   pos_arr[snp_cols],
    'ref':   [str(x) for x in ref_arr[snp_cols]],
    'alt':   [str(x) for x in alt_arr[snp_cols]],
    'AC_orig_cnvar': own_AC,
    'AC_coord_aware': new_AC,
    'extras_only':    extras_only,
})
# Read cactus_em alt_freq for these SNPs
cem_tsv = pd.read_csv(ROOT / 'scratch/v3qc_v3_mixedloose_chr1/SEEDMIX_S1.tsv', sep='\t')
snp_meta['alt_freq'] = cem_tsv.iloc[snp_cols]['alt_freq'].values

F1 = pd.read_csv(ROOT / 'scratch/hetmask_diag/f1_test_joined.tsv', sep='\t')
M = F1.merge(snp_meta, on=['chrom','pos','ref','alt'], how='inner', suffixes=('','_x'))
print(f'    joined: {len(M):,} rows')
# Sanity
assert np.allclose(M['alt_freq'], M['alt_freq_x'], atol=1e-5)
M = M.drop(columns=['alt_freq_x'])
M['cell'] = (M['high_fmiss'].astype(int)*2 + M['is_mixed_bubble']).map({
    0:'a',1:'b',2:'c',3:'d'})

# AF: under uniform-h, AF_coord_aware = AC_coord_aware / AN_post
M['recipe_coord_aware'] = M['AC_coord_aware'] / M['AN_post']
# cactus_em variant: AF_cem + extras/AN_post (uniform-h approx)
M['cactus_em_coord_aware'] = (M['alt_freq'] + M['extras_only'] / M['AN_post']).clip(0,1)

print('\n=== Per-cell agreement (hapFIRE − estimator) ===')
rows = []
for cell, sub in M.groupby('cell'):
    r_o = sub['hapfire_af'] - sub['alt_freq']
    r_b = sub['hapfire_af'] - sub['cactus_em_coord_aware']
    rows.append({
        'cell': cell, 'n': len(sub),
        'mean_orig':        r_o.mean(),
        'mean_coordaware':  r_b.mean(),
        'MAE_orig':         r_o.abs().mean(),
        'MAE_coordaware':   r_b.abs().mean(),
        'mae_reduction_%':  (1 - r_b.abs().mean()/r_o.abs().mean())*100,
        'out05_orig':       (r_o.abs()>0.05).mean()*100,
        'out05_coordaware': (r_b.abs()>0.05).mean()*100,
    })
print(pd.DataFrame(rows).round(4).to_string(index=False))

print('\n=== Genome-wide totals (516k joined Chr1 SNPs) ===')
for est, name in [('alt_freq',              'cactus_em (orig)'),
                  ('cactus_em_coord_aware', 'cactus_em + coord-aware'),
                  ('recipe_af',             'recipe (orig)'),
                  ('recipe_coord_aware',    'recipe + coord-aware')]:
    r = M['hapfire_af'] - M[est]
    print(f'  {name:32s}  mean={r.mean():+.5f}  MAE={r.abs().mean():.5f}  '
          f'R²={1 - (r**2).sum()/((M["hapfire_af"]-M["hapfire_af"].mean())**2).sum():.5f}  '
          f'%|d|>0.05={(r.abs()>0.05).mean()*100:.2f}%  %|d|>0.10={(r.abs()>0.10).mean()*100:.2f}%')

# Spot-checks
print('\n=== Spot-check ===')
for p_spot in [13843898, 10421645]:
    sp_d = M[(M['chrom']=='Chr1') & (M['pos']==p_spot)]
    if len(sp_d):
        print(sp_d[['pos','ref','alt','AC_orig_cnvar','AC_coord_aware','extras_only',
                    'AN_post','alt_freq','cactus_em_coord_aware','hapfire_af']].to_string(index=False))

M.to_csv(ROOT / 'scratch/f2_honest/coord_aware_full.tsv', sep='\t', index=False)
print(f'\nWrote coord_aware_full.tsv')
