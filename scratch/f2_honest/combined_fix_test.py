#!/usr/bin/env python3
"""
Combine F2-honest (raw VCF lookup) + base-at-position (sibling coalescence).

For each Chr1 SNP record at (chrom, pos, REF, ALT):
  1. current_carriers = cn_var[:, snp_col]  (231-vec, 0/1)
  2. base_at_pos_carriers = OR over sibling cn_var cols at same (chrom, pos)
                            where sibling ALT[0] == this ALT[0]
  3. f2_honest_carriers = raw deconstruct VCF biallelic at (chrom, pos, REF, ALT),
                          cactus founder positions where raw GT == 1
  4. combined_override_carriers = (base_at_pos OR f2_honest) AND NOT current
  5. New AC = current_AC + combined_override_count
  6. New AN = original AN_post (unchanged)
"""
import numpy as np
import pandas as pd
import scipy.sparse as sp
from pathlib import Path
import time

ROOT = Path('/carnegie/nobackup/scratch/tbellagio/hapfire_sv')

# ----- 1. Load cn_var + meta -----
print('[1] Load cn_var + meta')
cv = sp.load_npz(ROOT / 'poolfreq/data/cn_var_231_v3qc_v3.cn_var.npz').tocsc()
m  = np.load(ROOT / 'poolfreq/data/cn_var_231_v3qc_v3.meta.npz', allow_pickle=True)
chrom_arr   = m['chrom']; pos_arr = m['pos']
ref_arr     = m['ref'];   alt_arr = m['alt']
ref_len_arr = m['ref_len']; alt_len_arr = m['alt_len']
founders    = list(m['founders'])

# ----- 2. Load raw VCF GT carrier dict (for F2-honest) -----
print('[2] Load raw VCF carrier dict')
GT_TSV = ROOT / 'scratch/f2_honest/raw_chr1_biallelic_GT.tsv'
# Re-derive Assembly→Accession mapping and raw→cn_var row index
xl = pd.read_excel(ROOT / 'data/ASSEMBLIES_Best_version_of_dataset.xlsx',
                   usecols=['Assembly_ID','Accession_ID'])
xl['Assembly_ID'] = xl['Assembly_ID'].astype(str); xl['Accession_ID'] = xl['Accession_ID'].astype(str)
asm_to_acc = dict(zip(xl['Assembly_ID'], xl['Accession_ID']))
import gzip
with gzip.open(ROOT / 'scratch/f2_honest/raw_chr1.vcf.gz', 'rt') as gz:
    for line in gz:
        if line.startswith('#CHROM'):
            raw_samples = line.rstrip().split('\t')[9:]
            break
founder_idx = {f: i for i, f in enumerate(founders)}
raw_col_to_cnvar_row = []
for asm_id in raw_samples:
    acc_id = asm_to_acc.get(asm_id)
    cn_idx = founder_idx.get(acc_id, -1) if acc_id else -1
    raw_col_to_cnvar_row.append(cn_idx)
raw_col_to_cnvar_row = np.array(raw_col_to_cnvar_row, dtype=np.int32)
valid_raw = raw_col_to_cnvar_row >= 0
valid_raw_idx = np.where(valid_raw)[0]
valid_cnvar_rows = raw_col_to_cnvar_row[valid_raw]
print(f'    raw VCF cactus founders mapped to cn_var rows: {valid_raw.sum()}')

# Build (pos, ref, alt) -> 80-vec int8 carrier mask (for the 80 mapped cactus founders)
print('    streaming GT TSV (SNP-shape rows only) ...')
carrier_dict = {}
t0 = time.time()
n_streamed = 0
with open(GT_TSV) as f:
    for line in f:
        parts = line.rstrip('\n').split('\t')
        if len(parts[2]) != 1 or len(parts[3]) != 1:
            continue
        if parts[2] == '*' or parts[3] == '*':
            continue
        pos = int(parts[1]); ref = parts[2]; alt = parts[3]
        gts = parts[4:]
        # Only mapped (valid) raw columns
        carrier = np.zeros(len(valid_raw_idx), dtype=np.int8)
        for i, ri in enumerate(valid_raw_idx):
            g = gts[ri]
            if g == '1':
                carrier[i] = 1
        key = (pos, ref, alt)
        if key in carrier_dict:
            prev = carrier_dict[key]
            carrier_dict[key] = (prev | carrier).astype(np.int8)
        else:
            carrier_dict[key] = carrier
        n_streamed += 1
print(f'    SNP keys collected: {len(carrier_dict):,}  ({time.time()-t0:.0f}s)')

# ----- 3. Build position groups (Chr1) -----
print('[3] Build (chrom, pos) → cn_var col list (Chr1)')
chr1_mask = (chrom_arr == 'Chr1')
chr1_cols = np.where(chr1_mask)[0]
pos_chr1   = pos_arr[chr1_cols]
order      = np.argsort(pos_chr1, kind='stable')
cols_sorted= chr1_cols[order]
pos_sorted = pos_chr1[order]
boundaries = np.concatenate([[0], np.where(np.diff(pos_sorted) != 0)[0]+1, [len(pos_sorted)]])
pos_to_cols = {}
for i in range(len(boundaries)-1):
    p = int(pos_sorted[boundaries[i]])
    pos_to_cols[p] = cols_sorted[boundaries[i]:boundaries[i+1]]
print(f'    distinct positions: {len(pos_to_cols):,}')

chr1_col_to_idx = {c: i for i, c in enumerate(chr1_cols)}
chr1_alt_first = np.array([str(alt_arr[c])[0] if alt_arr[c] is not None else '.' for c in chr1_cols])

# ----- 4. Compute combined override per Chr1 SNP col -----
print('[4] Combined patch: base-at-pos OR F2-honest per SNP col')
snp_mask = chr1_mask & (ref_len_arr == 1) & (alt_len_arr == 1)
snp_cols = np.where(snp_mask)[0]
n_snp = len(snp_cols)
print(f'    Chr1 SNP cols: {n_snp:,}')

# Storage
n_extra_combined = np.zeros(n_snp, dtype=np.int32)
n_extra_baseatpos = np.zeros(n_snp, dtype=np.int32)
n_extra_f2honest = np.zeros(n_snp, dtype=np.int32)
n_overlap = np.zeros(n_snp, dtype=np.int32)

t0 = time.time()
for ix, snp_col in enumerate(snp_cols):
    pos = int(pos_arr[snp_col])
    ref = str(ref_arr[snp_col]); alt = str(alt_arr[snp_col])
    target_base = alt[0]
    # Current carriers
    own = np.asarray(cv[:, snp_col].todense()).ravel().astype(np.int8)
    # Base-at-pos: OR over sibling cn_var cols
    base_extra_mask = np.zeros(231, dtype=np.int8)
    cols_at_pos = pos_to_cols[pos]
    siblings = []
    for c in cols_at_pos:
        if c == snp_col:
            continue
        if chr1_alt_first[chr1_col_to_idx[c]] == target_base:
            siblings.append(c)
    if siblings:
        sub = cv[:, siblings]
        base_extra_mask = np.asarray((sub.sum(axis=1) > 0).astype(np.int8)).ravel()
    # F2-honest: lookup raw VCF biallelic at (pos, ref, alt); apply to mapped cactus rows
    f2_extra_mask = np.zeros(231, dtype=np.int8)
    raw_vec = carrier_dict.get((pos, ref, alt))
    if raw_vec is not None:
        for i, cn_row in enumerate(valid_cnvar_rows):
            if raw_vec[i] == 1:
                f2_extra_mask[cn_row] = 1
    # Combined override: (base OR f2) AND NOT own
    combined = ((base_extra_mask | f2_extra_mask) & (own == 0)).astype(np.int8)
    # Counts
    base_only = ((base_extra_mask & (own == 0)) & (f2_extra_mask == 0)).sum()
    f2_only   = ((f2_extra_mask & (own == 0)) & (base_extra_mask == 0)).sum()
    both      = ((base_extra_mask & f2_extra_mask) & (own == 0)).sum()
    n_extra_baseatpos[ix] = base_only + both
    n_extra_f2honest[ix]  = f2_only + both
    n_overlap[ix]         = both
    n_extra_combined[ix]  = combined.sum()
    if (ix+1) % 200000 == 0:
        print(f'    processed {ix+1:,} ({time.time()-t0:.0f}s)')
print(f'    done in {time.time()-t0:.0f}s')

# ----- 5. Aggregate & re-project AF -----
print('[5] Re-project AF and compare')
F1 = pd.read_csv(ROOT / 'scratch/hetmask_diag/f1_test_joined.tsv', sep='\t')

cem_tsv = pd.read_csv(ROOT / 'scratch/v3qc_v3_mixedloose_chr1/SEEDMIX_S1.tsv', sep='\t')
snp_meta = pd.DataFrame({
    'chrom': chrom_arr[snp_cols], 'pos': pos_arr[snp_cols],
    'ref':   [str(x) for x in ref_arr[snp_cols]],
    'alt':   [str(x) for x in alt_arr[snp_cols]],
    'alt_freq': cem_tsv.iloc[snp_cols]['alt_freq'].values,
    'extra_combined':  n_extra_combined,
    'extra_baseatpos': n_extra_baseatpos,
    'extra_f2honest':  n_extra_f2honest,
    'extra_overlap':   n_overlap,
})

M = F1.merge(snp_meta, on=['chrom','pos','ref','alt'], how='inner', suffixes=('','_x'))
M['cell'] = (M['high_fmiss'].astype(int)*2 + M['is_mixed_bubble']).map({0:'a',1:'b',2:'c',3:'d'})
print(f'    joined: {len(M):,}')

# AF projections
M['cactus_em_combined'] = (M['alt_freq'] + M['extra_combined'] / M['AN_post']).clip(0,1)
M['cactus_em_baseatpos']= (M['alt_freq'] + M['extra_baseatpos'] / M['AN_post']).clip(0,1)
M['cactus_em_f2honest'] = (M['alt_freq'] + M['extra_f2honest'] / M['AN_post']).clip(0,1)
M['recipe_combined']    = (M['AC_post'] + M['extra_combined']) / M['AN_post']

print('\n=== Per-cell agreement vs hapFIRE ===')
rows = []
for cell, sub in M.groupby('cell'):
    if cell not in ('a','b','c','d'): continue
    r_orig = sub['hapfire_af'] - sub['alt_freq']
    r_base = sub['hapfire_af'] - sub['cactus_em_baseatpos']
    r_f2   = sub['hapfire_af'] - sub['cactus_em_f2honest']
    r_comb = sub['hapfire_af'] - sub['cactus_em_combined']
    rows.append({
        'cell': cell, 'n': len(sub),
        'mean_orig':    r_orig.mean(),    'MAE_orig':    r_orig.abs().mean(),
        'mean_f2hon':   r_f2.mean(),      'MAE_f2hon':   r_f2.abs().mean(),
        'mean_base':    r_base.mean(),    'MAE_base':    r_base.abs().mean(),
        'mean_comb':    r_comb.mean(),    'MAE_comb':    r_comb.abs().mean(),
        'out05_orig':   (r_orig.abs()>0.05).mean()*100,
        'out05_comb':   (r_comb.abs()>0.05).mean()*100,
    })
T = pd.DataFrame(rows)
print(T.round(4).to_string(index=False))

print('\n=== Genome-wide (516k joined Chr1 SNPs) ===')
for est, name in [
    ('alt_freq',           'cactus_em (orig)'),
    ('cactus_em_f2honest', 'cactus_em + F2-honest'),
    ('cactus_em_baseatpos','cactus_em + base-at-pos'),
    ('cactus_em_combined', 'cactus_em + BOTH (combined)'),
]:
    r = M['hapfire_af'] - M[est]
    print(f'  {name:32s}  mean={r.mean():+.5f}  MAE={r.abs().mean():.5f}  '
          f'%|d|>0.05={(r.abs()>0.05).mean()*100:.2f}%  %|d|>0.10={(r.abs()>0.10).mean()*100:.2f}%')

# Overlap diagnostics
print('\n=== Overlap diagnostics (cells where overrides fire) ===')
M['extra_base_pos'] = (M['extra_baseatpos']>0).astype(int)
M['extra_f2_pos']   = (M['extra_f2honest']>0).astype(int)
print(pd.crosstab(M['extra_base_pos'], M['extra_f2_pos'],
                 rownames=['base-at-pos extras > 0'],
                 colnames=['F2-honest extras > 0']).to_string())

print('\n=== Per-cell overlap (cell [d]) ===')
D = M[M['cell']=='d']
print(f'  n cell [d]: {len(D):,}')
print(f'  records with both base AND F2 extras: {((D["extra_baseatpos"]>0) & (D["extra_f2honest"]>0)).sum()}')
print(f'  records with only base-at-pos extras: {((D["extra_baseatpos"]>0) & (D["extra_f2honest"]==0)).sum()}')
print(f'  records with only F2-honest extras  : {((D["extra_baseatpos"]==0) & (D["extra_f2honest"]>0)).sum()}')
print(f'  records with NO extras             : {((D["extra_baseatpos"]==0) & (D["extra_f2honest"]==0)).sum()}')

# Spot-check
print('\n=== Spot-check: Chr1:13843898 ===')
spd = M[(M['chrom']=='Chr1') & (M['pos']==13843898)]
if len(spd):
    print(spd[['chrom','pos','ref','alt','AC_post','AN_post',
               'extra_baseatpos','extra_f2honest','extra_combined','extra_overlap',
               'alt_freq','cactus_em_baseatpos','cactus_em_f2honest','cactus_em_combined',
               'hapfire_af']].to_string(index=False))

M.to_csv(ROOT / 'scratch/f2_honest/combined_fix_test.tsv', sep='\t', index=False)
print('\nWrote combined_fix_test.tsv')
