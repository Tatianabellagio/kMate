#!/usr/bin/env python3
"""
Test: rebuild cn_var with base-at-position semantics on Chr1, re-project AF,
compare to hapFIRE.

Algorithm:
  - For each Chr1 SNP cn_var record at (chrom, pos, REF, ALT):
    target_base = ALT[0]
  - Look up all other cn_var records at the same (chrom, pos)
  - For each (REF_i, ALT_i) record: this record's alt-carriers have base ALT_i[0] at pos
  - Sum carrier mass across all records where ALT_i[0] == target_base
  - new cn_var_carrier_vec = element-wise OR of those columns
    (founders are mutually exclusive across biallelics at the same position
     since they took exactly one ALT in the original multi-allelic)
  - new AC = sum(new carriers); AN unchanged.
"""
import numpy as np
import pandas as pd
import scipy.sparse as sp
from pathlib import Path
import time

ROOT = Path('/carnegie/nobackup/scratch/tbellagio/kmate')

# ----- 1. Load cn_var + meta -----
print('[1] Load cn_var + meta')
cv = sp.load_npz(ROOT / 'data/cn_var_231_v3qc_v3.cn_var.npz').tocsc()
m  = np.load(ROOT / 'data/cn_var_231_v3qc_v3.meta.npz', allow_pickle=True)
chrom_arr   = m['chrom']
pos_arr     = m['pos']
ref_arr     = m['ref']
alt_arr     = m['alt']
ref_len_arr = m['ref_len']
alt_len_arr = m['alt_len']
print(f'    cn_var: {cv.shape}, {cv.nnz:,} nnz')

# Restrict to Chr1
chr1_mask = (chrom_arr == 'Chr1')
chr1_cols = np.where(chr1_mask)[0]
print(f'    Chr1 cols: {len(chr1_cols):,}')

# Build position → list of col indices (Chr1)
print('[2] Build (chrom, pos) → list of cn_var column indices')
t0 = time.time()
pos_chr1   = pos_arr[chr1_cols]
order      = np.argsort(pos_chr1, kind='stable')
cols_sorted= chr1_cols[order]
pos_sorted = pos_chr1[order]
boundaries = np.concatenate([[0],
                             np.where(np.diff(pos_sorted) != 0)[0] + 1,
                             [len(pos_sorted)]])
pos_to_cols = {}
for i in range(len(boundaries)-1):
    p = int(pos_sorted[boundaries[i]])
    cols = cols_sorted[boundaries[i]:boundaries[i+1]]
    pos_to_cols[p] = cols
print(f'    distinct Chr1 positions: {len(pos_to_cols):,}  (took {time.time()-t0:.1f}s)')

# ----- 3. For each Chr1 SNP record, compute base-at-pos coalescence -----
print('[3] Apply base-at-pos coalescence to Chr1 SNP records')
snp_mask  = chr1_mask & (ref_len_arr == 1) & (alt_len_arr == 1)
snp_cols  = np.where(snp_mask)[0]
print(f'    Chr1 SNP cn_var cols: {len(snp_cols):,}')

# For each SNP col, find the carrier-mass-extra by coalescing with same base[0] records
extra_carriers_per_col = np.zeros(len(snp_cols), dtype=np.int32)
n_coalesced_recs_per_col = np.zeros(len(snp_cols), dtype=np.int32)
t0 = time.time()

# Precompute first base of each Chr1 record's ALT (cheap)
chr1_alt_first = np.array([str(alt_arr[c])[0] if alt_arr[c] is not None else '.' for c in chr1_cols])
chr1_col_to_idx = {c: i for i, c in enumerate(chr1_cols)}

for ix, snp_col in enumerate(snp_cols):
    pos = int(pos_arr[snp_col])
    target_base = str(alt_arr[snp_col])[0]   # ALT[0] for SNP
    cols_at_pos = pos_to_cols[pos]
    if len(cols_at_pos) == 1:
        continue  # no other records at this pos; nothing to coalesce
    # Find sibling records with ALT[0] == target_base AND different col index
    siblings = []
    for c in cols_at_pos:
        if c == snp_col:
            continue
        ci = chr1_col_to_idx[c]
        if chr1_alt_first[ci] == target_base:
            siblings.append(c)
    if not siblings:
        continue
    n_coalesced_recs_per_col[ix] = len(siblings)
    # Sum carriers across sibling cn_var columns
    # Founders are mutually exclusive across these — but safer to sum then clip to {0,1}
    sub = cv[:, siblings]                                   # (231, n_sib) CSC slice
    sib_or = np.asarray((sub.sum(axis=1) > 0).astype(np.int8)).ravel()  # 231-vec
    # Founders already in snp_col's cn_var=1 should not double-count
    own = np.asarray(cv[:, snp_col].todense()).ravel().astype(np.int8)
    new_carriers = (sib_or & (own == 0))
    extra_carriers_per_col[ix] = int(new_carriers.sum())
    if (ix + 1) % 200000 == 0:
        print(f'    processed {ix+1:,} / {len(snp_cols):,} SNP cols ({time.time()-t0:.0f}s)')

print(f'    done in {time.time()-t0:.0f}s')
print(f'    SNP cols with ≥1 coalesced sibling: {(n_coalesced_recs_per_col>0).sum():,}'
      f' / {len(snp_cols):,}  ({100*(n_coalesced_recs_per_col>0).mean():.1f}%)')
print(f'    extra_carriers_per_col stats:')
arr = extra_carriers_per_col[extra_carriers_per_col > 0]
print(f'      n>0: {len(arr):,}  mean: {arr.mean():.2f}  median: {int(np.median(arr))}  max: {arr.max()}')

# ----- 4. Re-project AF with the base-at-pos coalescence -----
print('[4] Re-project AF and compare to hapFIRE / xwu / cactus_em current')
F1 = pd.read_csv(ROOT / 'scratch/hetmask_diag/f1_test_joined.tsv', sep='\t')

# Build a per-SNP-record DataFrame with extra_carriers info
snp_meta = pd.DataFrame({
    'col':   snp_cols,
    'chrom': chrom_arr[snp_cols],
    'pos':   pos_arr[snp_cols],
    'ref':   [str(x) for x in ref_arr[snp_cols]],
    'alt':   [str(x) for x in alt_arr[snp_cols]],
    'extra_carriers': extra_carriers_per_col,
    'n_coalesced_recs': n_coalesced_recs_per_col,
})
cem_tsv = pd.read_csv(ROOT / 'scratch/v3qc_v3_mixedloose_chr1/SEEDMIX_S1.tsv', sep='\t')
snp_meta['alt_freq_orig'] = cem_tsv.iloc[snp_cols]['alt_freq'].values

M = F1.merge(snp_meta, on=['chrom','pos','ref','alt'], how='inner')
M['cell'] = (M['high_fmiss'].astype(int)*2 + M['is_mixed_bubble']).map({
    0:'a',1:'b',2:'c',3:'d'})
print(f'    joined: {len(M):,} rows')

# RECIPE-space exact base-at-pos AF:
#   AC_new = AC_post + extra_carriers  (extra = same-base carriers we were missing)
#   AN_new = AN_post                    (unchanged)
M['recipe_base'] = (M['AC_post'] + M['extra_carriers']) / M['AN_post']

# cactus_em-space approximation (uniform-h)
M['cactus_em_base'] = (M['alt_freq'] + M['extra_carriers'] / M['AN_post']).clip(0,1)

# ----- 5. Per-cell agreement vs hapFIRE -----
def summarize(df, label='all'):
    rows = []
    for cell, sub in df.groupby('cell'):
        if cell not in ('a','b','c','d'): continue
        r_orig = sub['hapfire_af'] - sub['alt_freq']
        r_base = sub['hapfire_af'] - sub['cactus_em_base']
        rows.append({
            'cell': cell, 'n': len(sub),
            'mean_orig':    r_orig.mean(),
            'MAE_orig':     r_orig.abs().mean(),
            'mean_baseatpos': r_base.mean(),
            'MAE_baseatpos':  r_base.abs().mean(),
            'out05_orig':   (r_orig.abs()>0.05).mean()*100,
            'out05_base':   (r_base.abs()>0.05).mean()*100,
            'reduction_%':  (1 - r_base.abs().mean()/r_orig.abs().mean())*100,
        })
    return pd.DataFrame(rows)

print('\n=== Per-cell agreement: cactus_em (orig) vs cactus_em (base-at-pos) vs hapFIRE ===')
print(summarize(M).round(4).to_string(index=False))

print('\n=== Genome-wide totals (all 516k joined Chr1 SNPs) ===')
for est, name in [('alt_freq',        'cactus_em (orig)'),
                  ('cactus_em_F2honest','cactus_em + F2-honest'),
                  ('cactus_em_base',  'cactus_em + base-at-pos'),
                  ('recipe_af',       'recipe (orig)'),
                  ('recipe_F2honest', 'recipe + F2-honest'),
                  ('recipe_base',     'recipe + base-at-pos')]:
    if est not in M.columns:
        # Need F2-honest cols if missing
        if est == 'cactus_em_F2honest':
            f2h = pd.read_csv(ROOT/'scratch/f2_honest/f2_honest_test.tsv', sep='\t',
                              usecols=['chrom','pos','ref','alt','cactus_em_F2honest','recipe_F2honest'])
            M = M.merge(f2h, on=['chrom','pos','ref','alt'], how='left')
        if est not in M.columns:
            continue
    r = M['hapfire_af'] - M[est]
    print(f'  {name:32s}  mean={r.mean():+.5f}  MAE={r.abs().mean():.5f}  '
          f'%|d|>0.05={(r.abs()>0.05).mean()*100:.2f}%  %|d|>0.10={(r.abs()>0.10).mean()*100:.2f}%')

# Save patched table
M.to_csv(ROOT / 'scratch/f2_honest/base_at_pos_test.tsv', sep='\t', index=False)
print(f'\nWrote {ROOT / "scratch/f2_honest/base_at_pos_test.tsv"}')

# Spot-check: Chr1:13843898 — the trace position
print('\n=== Spot-check: Chr1:13843898 (the trace position) ===')
sp_d = M[(M['chrom']=='Chr1') & (M['pos']==13843898)]
if len(sp_d):
    print(sp_d[['chrom','pos','ref','alt','AC_post','AN_post','extra_carriers','n_coalesced_recs',
                'alt_freq','cactus_em_base','recipe_base','hapfire_af']].to_string(index=False))
