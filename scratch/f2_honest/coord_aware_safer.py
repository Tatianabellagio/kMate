#!/usr/bin/env python3
"""
SAFER coord-aware coalescence: only allow coalesce with records where
ref_len == alt_len (clean SNP/MNP). Skip indel-containing records because
their per-character ALT-to-coord mapping is ambiguous.

This trades coalesce coverage for correctness:
  - keeps the legitimate cases (MNPs that share a base at the target coord)
  - excludes the over-coalescence (indel-containing records that coincidentally
    have a matching character at the target offset)
"""
import numpy as np
import pandas as pd
import scipy.sparse as sp
from pathlib import Path
import time

ROOT = Path('/carnegie/nobackup/scratch/tbellagio/kmate')

print('[1] Load cn_var + meta')
cv = sp.load_npz(ROOT / 'poolfreq/data/cn_var_231_v3qc_v3.cn_var.npz').tocsc()
m  = np.load(ROOT / 'poolfreq/data/cn_var_231_v3qc_v3.meta.npz', allow_pickle=True)
chrom_arr   = m['chrom']; pos_arr = m['pos']
ref_arr     = m['ref'];   alt_arr = m['alt']
ref_len_arr = m['ref_len']; alt_len_arr = m['alt_len']

# Chr1 + SAFER restriction: only ref_len == alt_len sources
chr1_mask = (chrom_arr == 'Chr1')
clean_mask = chr1_mask & (ref_len_arr == alt_len_arr)
chr1_cols = np.where(chr1_mask)[0]
clean_cols = np.where(clean_mask)[0]
print(f'    Chr1 records (all):              {len(chr1_cols):,}')
print(f'    Chr1 records (ref_len==alt_len): {len(clean_cols):,}')
print(f'    fraction "clean" (MNP/SNP):      {100*len(clean_cols)/len(chr1_cols):.1f}%')

snp_mask = chr1_mask & (ref_len_arr == 1) & (alt_len_arr == 1)
snp_cols = np.where(snp_mask)[0]

# SNP index
snp_pos_array = pos_arr[snp_cols]
snp_alt_array = [str(a) for a in alt_arr[snp_cols]]
snp_index = {}
for i, col in enumerate(snp_cols):
    key = (int(snp_pos_array[i]), snp_alt_array[i])
    snp_index.setdefault(key, []).append(i)
snp_pos_sorted_uniq = np.array(sorted(set(int(p) for p in snp_pos_array)))

# Source record arrays (CLEAN only)
src_pos = pos_arr[clean_cols].astype(np.int64)
src_ref_len = ref_len_arr[clean_cols].astype(np.int32)
src_alt_str = [str(a) for a in alt_arr[clean_cols]]
print(f'    max ref_len among "clean" sources: {src_ref_len.max()}')
print(f'    median ref_len: {int(np.median(src_ref_len))}')

# Extras mask
n_snps = len(snp_cols)
extras_mask = np.zeros((n_snps, 231), dtype=np.int8)

print('[2] Iterate clean records, find spanning SNPs, OR-merge carriers')
t0 = time.time()
n_carrier_merges = 0
for ix, src_col in enumerate(clean_cols):
    p = int(src_pos[ix]); rl = int(src_ref_len[ix]); alt = src_alt_str[ix]
    lo = np.searchsorted(snp_pos_sorted_uniq, p, side='left')
    hi = np.searchsorted(snp_pos_sorted_uniq, p + rl, side='left')
    if hi <= lo:
        continue
    src_carriers = None
    for target_pos in snp_pos_sorted_uniq[lo:hi]:
        offset = int(target_pos) - p
        if offset >= len(alt):
            continue
        target_alt_base = alt[offset]
        key = (int(target_pos), target_alt_base)
        snp_locals = snp_index.get(key)
        if snp_locals is None:
            continue
        if src_carriers is None:
            src_carriers = np.asarray(cv[:, src_col].todense()).ravel().astype(np.int8)
        for sl in snp_locals:
            if snp_cols[sl] == src_col:
                continue
            extras_mask[sl] = np.maximum(extras_mask[sl], src_carriers)
            n_carrier_merges += 1
    if (ix+1) % 200000 == 0:
        print(f'    processed {ix+1:,}/{len(clean_cols):,}  ({time.time()-t0:.0f}s)  merges={n_carrier_merges:,}')
print(f'    DONE in {time.time()-t0:.0f}s. total merges = {n_carrier_merges:,}')

# Compute new AC
print('[3] Compute new AC')
own_AC = np.zeros(n_snps, dtype=np.int32)
new_AC = np.zeros(n_snps, dtype=np.int32)
extras_only = np.zeros(n_snps, dtype=np.int32)
for i, snp_col in enumerate(snp_cols):
    own = np.asarray(cv[:, snp_col].todense()).ravel().astype(np.int8)
    own_AC[i] = int(own.sum())
    combined = np.maximum(own, extras_mask[i])
    new_AC[i] = int(combined.sum())
    extras_only[i] = int(((extras_mask[i] == 1) & (own == 0)).sum())

print(f'    SNPs with extras > 0: {(extras_only>0).sum():,}/{n_snps:,}  '
      f'({100*(extras_only>0).mean():.1f}%)')

# Build summary DF + join + project + compare
print('[4] Project + compare')
snp_meta = pd.DataFrame({
    'chrom': chrom_arr[snp_cols],
    'pos':   pos_arr[snp_cols],
    'ref':   [str(x) for x in ref_arr[snp_cols]],
    'alt':   [str(x) for x in alt_arr[snp_cols]],
    'AC_orig_cnvar': own_AC,
    'AC_safer':      new_AC,
    'extras_only_safer': extras_only,
})
cem_tsv = pd.read_csv(ROOT / 'scratch/v3qc_v3_mixedloose_chr1/SEEDMIX_S1.tsv', sep='\t')
snp_meta['alt_freq'] = cem_tsv.iloc[snp_cols]['alt_freq'].values

F1 = pd.read_csv(ROOT / 'scratch/hetmask_diag/f1_test_joined.tsv', sep='\t')
M = F1.merge(snp_meta, on=['chrom','pos','ref','alt'], how='inner', suffixes=('','_x'))
assert np.allclose(M['alt_freq'], M['alt_freq_x'], atol=1e-5)
M = M.drop(columns=['alt_freq_x'])
M['cell'] = (M['high_fmiss'].astype(int)*2 + M['is_mixed_bubble']).map({0:'a',1:'b',2:'c',3:'d'})

# Also pull in the previous coord_aware (unsafe) result for comparison
prev = pd.read_csv(ROOT / 'scratch/f2_honest/coord_aware_full.tsv', sep='\t',
                  usecols=['chrom','pos','ref','alt','cactus_em_coord_aware','recipe_coord_aware','extras_only'])
prev = prev.rename(columns={'cactus_em_coord_aware':'cem_unsafe',
                            'recipe_coord_aware':'rec_unsafe',
                            'extras_only':'extras_unsafe'})
M = M.merge(prev, on=['chrom','pos','ref','alt'], how='left')

M['cactus_em_safer'] = (M['alt_freq'] + M['extras_only_safer'] / M['AN_post']).clip(0,1)
M['recipe_safer']    = (M['AC_post'] + M['extras_only_safer']) / M['AN_post']

print('\n=== Per-cell agreement (hapFIRE − estimator) ===')
for cell, sub in M.groupby('cell'):
    rows = []
    r_o = sub['hapfire_af'] - sub['alt_freq']
    r_unsafe = sub['hapfire_af'] - sub['cem_unsafe']
    r_safe = sub['hapfire_af'] - sub['cactus_em_safer']
    print(f'cell [{cell}] n={len(sub):>6}  ',
          f'orig mean={r_o.mean():+.4f} MAE={r_o.abs().mean():.4f}  |',
          f'unsafe mean={r_unsafe.mean():+.4f} MAE={r_unsafe.abs().mean():.4f}  |',
          f'safer mean={r_safe.mean():+.4f} MAE={r_safe.abs().mean():.4f}')

print('\n=== Genome-wide totals ===')
for est, name in [('alt_freq',          'cactus_em (orig)'),
                  ('cem_unsafe',        'coord-aware UNSAFE (all records)'),
                  ('cactus_em_safer',   'coord-aware SAFER (ref_len==alt_len only)')]:
    r = M['hapfire_af'] - M[est]
    print(f'  {name:42s}  mean={r.mean():+.5f}  MAE={r.abs().mean():.5f}  '
          f'%|d|>0.05={(r.abs()>0.05).mean()*100:.2f}%  %|d|>0.10={(r.abs()>0.10).mean()*100:.2f}%')

# Check the over-coalescence rate
print('\n=== Over-coalescence check ===')
delta_unsafe = (M['hapfire_af']-M['cem_unsafe']).abs() - (M['hapfire_af']-M['alt_freq']).abs()
delta_safer  = (M['hapfire_af']-M['cactus_em_safer']).abs() - (M['hapfire_af']-M['alt_freq']).abs()
print(f'  UNSAFE: records improved={(delta_unsafe < -0.01).sum():,}  unchanged={(delta_unsafe.abs() <= 0.01).sum():,}  '
      f'WORSENED={(delta_unsafe > 0.01).sum():,}')
print(f'  SAFER : records improved={(delta_safer  < -0.01).sum():,}  unchanged={(delta_safer.abs()  <= 0.01).sum():,}  '
      f'WORSENED={(delta_safer  > 0.01).sum():,}')

# Spot-checks: did the safer version save the 5870018-style over-coalescence?
print('\n=== Spot-checks ===')
for p_spot in [13843898, 10421645, 5870018, 1937813]:
    sp_d = M[M['pos']==p_spot]
    if len(sp_d):
        print(sp_d[['pos','ref','alt','AC_orig_cnvar','AC_safer','extras_only_safer','extras_unsafe',
                    'alt_freq','cactus_em_safer','cem_unsafe','hapfire_af']].to_string(index=False))

M.to_csv(ROOT / 'scratch/f2_honest/coord_aware_safer.tsv', sep='\t', index=False)
