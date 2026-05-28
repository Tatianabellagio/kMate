#!/usr/bin/env python3
"""
Diagnostic: does the v3qc_v3 het-mask explain the hapFIRE vs cactus_em residual?

Joins per-record:
  - cactus_em SEEDMIX_S1 AF (v3qc_v3 mixed-loose, Chr1)
  - hapFIRE SEEDMIX_S1 SNP AF (s1_snp_frequency.txt)
  - pre-mask PanGenie AC_Het + AN  (pangenie_153_filled_bi.vcf.gz)

Restricted to biallelic Chr1 SNPs.
Computes residual = hapFIRE_af - cactus_em_af, then bins by pre-mask PG het rate.
Prediction: if het-mask is the primary driver of disagreement, |residual| and the
mean (hapFIRE - cactus_em) scale monotonically with pre-mask het rate.
"""
import numpy as np
import pandas as pd
from pathlib import Path

ROOT   = Path('/carnegie/nobackup/scratch/tbellagio/hapfire_sv')
META   = ROOT / 'poolfreq/data/cn_var_231_v3qc_v3.meta.npz'
CEMTSV = ROOT / 'scratch/v3qc_v3_mixedloose_chr1/SEEDMIX_S1.tsv'
HAPF   = '/global/scratch/projects/fc_moilab/projects/grenenet-phase1/frequency/hapFIRE_frequencies/seed_mix/s1_snp_frequency.txt'
PGPRE  = ROOT / 'scratch/hetmask_diag/pg_premask_chr1_snps.tsv'

print('[1] Loading cn_var meta...')
m = np.load(META, allow_pickle=True)
chrom    = m['chrom']
pos      = m['pos']
ref_len  = m['ref_len']
alt_len  = m['alt_len']
ref_b    = m['ref']
alt_b    = m['alt']
N = chrom.shape[0]
print(f'    cn_var has {N:,} records')

# Index of Chr1 biallelic-SNP records in cn_var. Row index == row index in TSV
# (after header). per-record SEEDMIX_S1.tsv was written by the per-sample driver
# in the same order as cn_var columns.
snp_mask = (chrom == 'Chr1') & (ref_len == 1) & (alt_len == 1)
snp_idx  = np.where(snp_mask)[0]
print(f'    Chr1 SNP records: {len(snp_idx):,}')

# Per-position biallelic count (used downstream to keep positions w/ a unique ALT)
pos_chr1_snp = pos[snp_idx]
# unique-pos -> count
unique_pos, inv, counts = np.unique(pos_chr1_snp, return_inverse=True, return_counts=True)
biallelic_mask = counts[inv] == 1     # row is biallelic at the SNP class
print(f'    Of those, biallelic (1 SNP-ALT per pos): {biallelic_mask.sum():,}')

# Build cactus_em Chr1 SNP table by reading TSV row 0-indexed = snp_idx (after header)
print('[2] Loading cactus_em TSV (Chr1 SNP rows only)...')
# fastest: read all then index
cem_df = pd.read_csv(CEMTSV, sep='\t')
assert len(cem_df) == N, f'TSV row count {len(cem_df)} != cn_var col count {N}'
cem_snp = cem_df.iloc[snp_idx].copy()
cem_snp['ref'] = ref_b[snp_idx]
cem_snp['alt'] = alt_b[snp_idx]
cem_snp['biallelic'] = biallelic_mask
print(f'    cactus_em SNP rows: {len(cem_snp):,}')
print('    head:'); print(cem_snp.head())

print('[3] Loading pre-mask PG Chr1 SNP table...')
pg = pd.read_csv(
    PGPRE, sep='\t', header=None,
    names=['chrom','pos','ref','alt','AC_pre','AN_pre','AC_Het','F_MISSING_pre']
)
# Coerce numeric -- bcftools may write `.` for missing
for c in ('AC_pre','AN_pre','AC_Het','F_MISSING_pre'):
    pg[c] = pd.to_numeric(pg[c], errors='coerce')
print(f'    PG pre-mask SNP rows: {len(pg):,}')

# Note: per-record AN_pre can be < 306 (=2*153) if PanGenie returned ./. for some samples
# het_rate = AC_Het / (number of called samples) = AC_Het / (AN_pre / 2)
pg['n_called'] = pg['AN_pre'] / 2.0
pg['het_rate'] = np.where(pg['n_called'] > 0, pg['AC_Het'] / pg['n_called'], 0.0)
# Pre-mask PG AF (as a sanity check)
pg['pg_af_pre'] = np.where(pg['AN_pre'] > 0, pg['AC_pre'] / pg['AN_pre'], np.nan)

print('[4] Joining cactus_em + PG pre-mask on (chrom,pos,ref,alt)...')
j = cem_snp.merge(pg, on=['chrom','pos','ref','alt'], how='inner')
print(f'    joined rows (any biallelic flag): {len(j):,}')
j = j[j['biallelic']].copy()
print(f'    biallelic Chr1 SNPs joined: {len(j):,}')

print('[5] Loading hapFIRE s1_snp_frequency.txt and joining on (pos)...')
hf = pd.read_csv(HAPF, sep='\t', header=None, names=['chrom_hf','pos','hapfire_af'])
hf = hf[hf['chrom_hf'] == 1][['pos','hapfire_af']].copy()
print(f'    hapFIRE Chr1 SNP rows: {len(hf):,}')

j = j.merge(hf, on='pos', how='inner')
print(f'    final 3-way joined rows (biallelic Chr1 SNPs): {len(j):,}')

print('[6] Computing residuals + summary...')
j['residual'] = j['hapfire_af'] - j['alt_freq']
j['abs_res']  = j['residual'].abs()

# Bin by het_rate
bins   = [-1e-9, 0.0, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0]
labels = ['=0', '(0,0.005]', '(0.005,0.01]', '(0.01,0.02]', '(0.02,0.05]',
          '(0.05,0.1]', '(0.1,0.2]', '(0.2,0.5]', '(0.5,1.0]']
j['het_bin'] = pd.cut(j['het_rate'], bins=bins, labels=labels)

summary = j.groupby('het_bin', observed=True).agg(
    n=('residual','size'),
    mean_residual=('residual','mean'),
    median_residual=('residual','median'),
    mean_abs_res=('abs_res','mean'),
    mean_hapfire=('hapfire_af','mean'),
    mean_cactus_em=('alt_freq','mean'),
    mean_pg_af_pre=('pg_af_pre','mean'),
).round(5)
print('\n=== Residual (hapFIRE - cactus_em) vs pre-mask PG het rate ===')
print(summary.to_string())

# Pearson + Spearman of residual vs het_rate
from scipy.stats import pearsonr, spearmanr
r_p = pearsonr(j['het_rate'], j['residual'])
r_s = spearmanr(j['het_rate'], j['residual'])
print(f'\nPearson  residual vs het_rate: r={r_p.statistic:+.4f}  p={r_p.pvalue:.2e}')
print(f'Spearman residual vs het_rate: r={r_s.statistic:+.4f}  p={r_s.pvalue:.2e}')

# Also: |residual| vs het_rate
r_pa = pearsonr(j['het_rate'], j['abs_res'])
r_sa = spearmanr(j['het_rate'], j['abs_res'])
print(f'\nPearson  |residual| vs het_rate: r={r_pa.statistic:+.4f}  p={r_pa.pvalue:.2e}')
print(f'Spearman |residual| vs het_rate: r={r_sa.statistic:+.4f}  p={r_sa.pvalue:.2e}')

# Simple linear regression of residual on het_rate
x = j['het_rate'].values
y = j['residual'].values
A = np.vstack([x, np.ones_like(x)]).T
slope, intercept = np.linalg.lstsq(A, y, rcond=None)[0]
print(f'\nresidual = {slope:+.4f} * het_rate + {intercept:+.5f}')

# Sanity: distribution of het_rate
print('\nDistribution of pre-mask het_rate at biallelic Chr1 SNPs:')
print(pd.Series(j['het_rate']).describe(percentiles=[0.5,0.9,0.95,0.99]).to_string())

# Save the joined table for re-use
out = ROOT / 'scratch/hetmask_diag/joined_hetmask_chr1.tsv'
j[['chrom','pos','ref','alt','alt_freq','hapfire_af','residual',
   'AC_Het','AN_pre','n_called','het_rate','pg_af_pre','F_MISSING_pre']].to_csv(out, sep='\t', index=False)
print(f'\nWrote {out}')

# Also: tail outliers
print('\nTop 10 |residual| rows (with their het_rate):')
print(j.nlargest(10, 'abs_res')[['chrom','pos','ref','alt','alt_freq','hapfire_af','residual',
                                 'AC_Het','AN_pre','n_called','het_rate','pg_af_pre','F_MISSING_pre']].to_string())
