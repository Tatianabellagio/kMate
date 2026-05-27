#!/usr/bin/env python3
"""
F1 test: mask SV-carrier founders as ./. at colocated SNP records.

Two evaluations:
  (a) recipe-space:    AF_F1_recipe   = AC / (AN - n_sv_carriers)
  (b) cactus_em-space: AF_F1_cactus  ≈ AF_cem * AN / (AN - n_sv_carriers)  (uniform-h approx)

Then re-evaluate vs hapFIRE on cells [a, b, c, d].
"""
import numpy as np
import pandas as pd
import scipy.sparse as sp
from pathlib import Path
from scipy.stats import pearsonr

ROOT = Path('/carnegie/nobackup/scratch/tbellagio/hapfire_sv')

print('[1] Load cn_var + meta')
cv = sp.load_npz(ROOT / 'poolfreq/data/cn_var_231_v3qc_v3.cn_var.npz')   # (231, 6.29M) CSR
m  = np.load(ROOT / 'poolfreq/data/cn_var_231_v3qc_v3.meta.npz', allow_pickle=True)
chrom_arr   = m['chrom']
pos_arr     = m['pos']
ref_len_arr = m['ref_len']
alt_len_arr = m['alt_len']
ref_arr     = m['ref']
alt_arr     = m['alt']
print(f'    cn_var: {cv.shape}, {cv.nnz:,} nnz')

print('[2] Restrict to Chr1, identify SNP vs non-SNP columns')
chr1_mask  = (chrom_arr == 'Chr1')
chr1_cols  = np.where(chr1_mask)[0]
snp_mask   = chr1_mask & (ref_len_arr == 1) & (alt_len_arr == 1)
nonsnp_mask= chr1_mask & ~((ref_len_arr == 1) & (alt_len_arr == 1))
snp_cols   = np.where(snp_mask)[0]
nonsnp_cols= np.where(nonsnp_mask)[0]
print(f'    Chr1 SNP cols    : {len(snp_cols):,}')
print(f'    Chr1 non-SNP cols: {len(nonsnp_cols):,}')

print('[3] Group non-SNP columns by position; compute sv_carrier_at_p (231-bit mask) per mixed-bubble pos')
# Convert cn_var to CSC for column slicing
print('    converting cn_var to CSC for column slicing...')
cv_csc = cv.tocsc()
print('    done.')

# Build dict: pos -> list of non-SNP col indices
nonsnp_pos = pos_arr[nonsnp_cols]
order = np.argsort(nonsnp_pos)
nonsnp_cols_sorted = nonsnp_cols[order]
nonsnp_pos_sorted  = nonsnp_pos[order]

# Group by pos
boundaries = np.concatenate([[0], np.where(np.diff(nonsnp_pos_sorted) != 0)[0]+1, [len(nonsnp_pos_sorted)]])
pos_to_nonsnp = {}
for i in range(len(boundaries)-1):
    p = int(nonsnp_pos_sorted[boundaries[i]])
    cols = nonsnp_cols_sorted[boundaries[i]:boundaries[i+1]]
    pos_to_nonsnp[p] = cols.tolist()
print(f'    distinct positions with non-SNP records on Chr1: {len(pos_to_nonsnp):,}')

# Build dict: pos -> SNP col indices (we only need SNP cols at positions that also have non-SNP records)
snp_pos = pos_arr[snp_cols]
order_s = np.argsort(snp_pos)
snp_cols_sorted = snp_cols[order_s]
snp_pos_sorted  = snp_pos[order_s]
boundaries_s = np.concatenate([[0], np.where(np.diff(snp_pos_sorted) != 0)[0]+1, [len(snp_pos_sorted)]])
pos_to_snp = {}
for i in range(len(boundaries_s)-1):
    p = int(snp_pos_sorted[boundaries_s[i]])
    cols = snp_cols_sorted[boundaries_s[i]:boundaries_s[i+1]]
    pos_to_snp[p] = cols.tolist()
print(f'    distinct positions with SNP records on Chr1: {len(pos_to_snp):,}')

mixed_pos = sorted(set(pos_to_nonsnp.keys()) & set(pos_to_snp.keys()))
print(f'    mixed-bubble positions (have both SNP + non-SNP)  : {len(mixed_pos):,}')

print('[4] Compute n_sv_carriers per mixed-bubble position')
# sv_carrier_at_p[f] = OR over non-SNP cols at p of cn_var[f, col]
n_sv_carriers_at_pos = {}
sv_carriers_at_pos   = {}  # boolean mask per pos for later h-weighted version
for p in mixed_pos:
    cols = pos_to_nonsnp[p]
    # cn_var (CSC) columns -> sum over founders (rows)
    sub = cv_csc[:, cols]                          # 231 × n_cols (CSC)
    # any across columns:
    carrier = (np.asarray(sub.sum(axis=1)).ravel() > 0)
    n_sv_carriers_at_pos[p] = int(carrier.sum())
    sv_carriers_at_pos[p] = carrier               # 231-vector boolean

print('    summary of n_sv_carriers_at_p:')
arr = np.array(list(n_sv_carriers_at_pos.values()))
print(f'      n mixed positions: {len(arr):,}')
print(f'      mean n_sv_carriers : {arr.mean():.2f}')
print(f'      median             : {np.median(arr):.0f}')
print(f'      max                : {arr.max()}')
print(f'      P(n_sv_carriers > 30): {(arr>30).mean()*100:.1f}%')
print(f'      P(n_sv_carriers > 100): {(arr>100).mean()*100:.1f}%')

print('[5] Load existing joined table (cell [a-d] labels), apply F1 to SNPs at mixed positions')
J2_path = ROOT / 'scratch/hetmask_diag/joined_hetmask_chr1.tsv'
J = pd.read_csv(J2_path, sep='\t')

fm = pd.read_csv(
    ROOT / 'scratch/hetmask_diag/panel231_chr1_snp_fmissing.tsv',
    sep='\t', header=None,
    names=['chrom','pos','ref','alt','AC_post','AN_post','F_MISSING_post']
)
for c in ('AC_post','AN_post','F_MISSING_post'):
    fm[c] = pd.to_numeric(fm[c], errors='coerce')
J = J.merge(fm, on=['chrom','pos','ref','alt'], how='inner')

mech2 = pd.read_csv(
    ROOT / 'scratch/v3qc_v3_mixedloose_chr1/mech2_chr1_joined.tsv.gz',
    sep='\t',
    usecols=['chrom','pos','ref','alt','cat','n_recs_at_pos']
)
J = J.merge(mech2, on=['chrom','pos','ref','alt'], how='inner')
J['is_mixed_bubble'] = (J['cat'] == 'C_mixed_bubble').astype(int)
J['high_fmiss']     = (J['F_MISSING_post'] >= 0.10)
print(f'    joined rows (Chr1 biallelic SNPs in all 3 tables): {len(J):,}')

# n_sv_carriers_at_p per record
J['n_sv_carriers_at_p'] = J['pos'].map(n_sv_carriers_at_pos).fillna(0).astype(int)
print(f'\n    cell-[d] records with n_sv_carriers_at_p available: '
      f'{((J["is_mixed_bubble"]==1) & (J["high_fmiss"]) & (J["n_sv_carriers_at_p"]>0)).sum():,}'
      f' / {((J["is_mixed_bubble"]==1) & (J["high_fmiss"])).sum():,}')

# Recipe AF (current) and F1-corrected recipe AF
J['recipe_af']   = np.where(J['AN_post']>0, J['AC_post']/J['AN_post'], np.nan)
J['AN_F1']       = J['AN_post'] - J['n_sv_carriers_at_p']
J['recipe_F1']   = np.where(J['AN_F1']>0, J['AC_post']/J['AN_F1'], np.nan)

# cactus_em F1 approximation
J['cactus_em_F1'] = np.where(J['AN_F1']>0,
                             J['alt_freq'] * J['AN_post'] / J['AN_F1'],
                             np.nan)
# clip to [0, 1]
J['cactus_em_F1'] = J['cactus_em_F1'].clip(0, 1)

# residuals
J['resid_hf']           = J['hapfire_af'] - J['alt_freq']
J['resid_hf_F1cactus']  = J['hapfire_af'] - J['cactus_em_F1']
J['resid_hf_F1recipe']  = J['hapfire_af'] - J['recipe_F1']
J['resid_hf_recipeOrig']= J['hapfire_af'] - J['recipe_af']

print('\n[6] Aggregate metrics, by 2x2 cell')
labels = {
    (False, False): '[a] low_FMISS  pure_SNP',
    (False, True):  '[b] low_FMISS  mixed_bub',
    (True,  False): '[c] high_FMISS pure_SNP',
    (True,  True):  '[d] high_FMISS mixed_bub',
}
rows = []
for (lo,mb), label in labels.items():
    sub = J[(J['high_fmiss']==lo) & (J['is_mixed_bubble']==mb)]
    rows.append({
        'cell': label,
        'n': len(sub),
        # before F1
        'mean_resid_OLD':     sub['resid_hf'].mean(),
        'MAE_OLD':            sub['resid_hf'].abs().mean(),
        'out05_OLD':          (sub['resid_hf'].abs()>0.05).mean()*100,
        # after F1 (cactus_em-space, approx)
        'mean_resid_F1cact':  sub['resid_hf_F1cactus'].mean(),
        'MAE_F1cact':         sub['resid_hf_F1cactus'].abs().mean(),
        'out05_F1cact':       (sub['resid_hf_F1cactus'].abs()>0.05).mean()*100,
        # recipe space (clean)
        'mean_resid_F1rec':   sub['resid_hf_F1recipe'].mean(),
        'MAE_F1rec':          sub['resid_hf_F1recipe'].abs().mean(),
        # mean correction applied
        'mean_corr_factor':   (sub['AN_post']/sub['AN_F1']).mean(),
    })
out = pd.DataFrame(rows).round(5)
print(out.to_string(index=False))

# Save
J.to_csv(ROOT / 'scratch/hetmask_diag/f1_test_joined.tsv', sep='\t', index=False)
print(f'\nWrote {ROOT / "scratch/hetmask_diag/f1_test_joined.tsv"}')

# Genome-wide aggregate
print('\n[7] Genome-wide totals (all 516k joined biallelic Chr1 SNPs):')
for est_col, name in [('alt_freq','cactus_em (orig)'),
                      ('cactus_em_F1','cactus_em + F1 (approx)'),
                      ('recipe_af','recipe (orig)'),
                      ('recipe_F1','recipe + F1')]:
    r = J['hapfire_af'] - J[est_col]
    print(f'  {name:28s}  mean={r.mean():+.5f}  MAE={r.abs().mean():.5f}  '
          f'%|d|>0.05={(r.abs()>0.05).mean()*100:.2f}%')
