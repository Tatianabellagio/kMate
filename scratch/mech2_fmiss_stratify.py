"""
Follow-up to mech2_bubble_diag.py: control for F_MISSING.

Predict (Mechanism 2 distinct from missingness): SNPs at positions shared with
non-SNP records show hapfire > cactus_em bias EVEN at low F_MISSING, because the
mechanism puts founders carrying the SV-path at GT=0 (called REF) on the SNP
record — not ./. — so the AC/AN denominator correction can't fix it.

Stream Chr1 records of the merged panel haploid VCF to extract F_MISSING per row,
then re-stratify the joined cactus_em / hapfire table.
"""
import sys, time, gzip
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path('/carnegie/nobackup/scratch/tbellagio/kmate')

t0 = time.time()
print('=== reading F_MISSING + LV for Chr1 from merged panel VCF ===', flush=True)
VCF = ROOT/'pangenie_genotyping/data/v3qc_v3/founders_231_v3qc_v3.haploid.vcf.gz'

# We need F_MISSING + LV for every Chr1 record, in the SAME order as cn_var meta
# (build_cn_var iterates the VCF top-to-bottom, no filter). So row index of Chr1
# records in the VCF == row index of Chr1 rows in meta. Verify by chrom/pos check.

import subprocess, io
BCF = '/home/tbellagio/miniforge3/envs/sequencing_pipeline/bin/bcftools'
# bcftools query — gives us just what we need (TSV)
cmd = [BCF, 'query', '-r', 'Chr1',
       '-f', '%CHROM\t%POS\t%REF\t%ALT\t%INFO/LV\t%INFO/F_MISSING\n',
       str(VCF)]
print(' '.join(cmd))
proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True)

rows = []
for line in proc.stdout:
    f = line.rstrip('\n').split('\t')
    # F_MISSING may be '.' if fill-tags didn't run on a record — coerce
    fm = f[5]
    try:
        fm = float(fm)
    except ValueError:
        fm = np.nan
    lv = f[4]
    try:
        lv = int(lv)
    except ValueError:
        lv = -1
    rows.append((int(f[1]), f[2], f[3], lv, fm))
proc.wait()
print(f'  Chr1 VCF rows: {len(rows):,}')
vcf_df = pd.DataFrame(rows, columns=['pos','ref','alt','LV','F_MISSING'])

# Filter VCF rows to SNP-shaped + ACGT (must align with the SNP filter in mech2 script)
vcf_df['ref_len'] = vcf_df['ref'].str.len()
vcf_df['alt_len'] = vcf_df['alt'].str.len()
ACGT = {'A','C','G','T'}
mask = (vcf_df['ref_len']==1) & (vcf_df['alt_len']==1) & \
       vcf_df['ref'].isin(ACGT) & vcf_df['alt'].isin(ACGT)
snp_vcf = vcf_df[mask][['pos','ref','alt','LV','F_MISSING']].copy()
print(f'  Chr1 SNP-shaped VCF rows: {len(snp_vcf):,}')

# Load joined mech2 table
m = pd.read_csv(ROOT/'scratch/v3qc_v3_mixedloose_chr1/mech2_chr1_joined.tsv.gz', sep='\t')
print(f'  mech2 joined rows: {len(m):,}')
m = m.merge(snp_vcf, on=['pos','ref','alt'], how='left')
print(f'  after merging LV/F_MISSING: {len(m):,}, F_MISSING NaN: {m["F_MISSING"].isna().sum():,}')

# Bin by F_MISSING
m['fmiss_bin'] = pd.cut(m['F_MISSING'],
                       bins=[-0.001, 0.01, 0.05, 0.1, 0.2, 0.5, 1.0],
                       labels=['<0.01','0.01-0.05','0.05-0.1','0.1-0.2','0.2-0.5','>=0.5'])

def stratify(label, sub):
    if len(sub) == 0:
        return None
    return dict(
        label=label, n=len(sub),
        MAE=sub['resid'].abs().mean(),
        mean=sub['resid'].mean(),
        pct_hfg=(sub['resid']>0.05).mean()*100,
        pct_cmg=(sub['resid']<-0.05).mean()*100,
    )

print('\n=== Cross-tab: F_MISSING bin x bubble category (mean residual hf-cem) ===')
piv_mean = m.pivot_table(values='resid', index='fmiss_bin', columns='cat', aggfunc='mean')
piv_n    = m.pivot_table(values='resid', index='fmiss_bin', columns='cat', aggfunc='count')
print('mean residual (hf - cem):')
print(piv_mean.to_string(float_format=lambda x: f'{x:+.4f}'))
print('\nn records:')
print(piv_n.to_string(float_format=lambda x: f'{x:.0f}'))

print('\n=== Cross-tab: F_MISSING bin x bubble category (% records with hf > cem + 0.05) ===')
m['_hfgt'] = (m['resid']>0.05).astype(int)
piv_hfg = m.pivot_table(values='_hfgt', index='fmiss_bin', columns='cat', aggfunc='mean')*100
print(piv_hfg.to_string(float_format=lambda x: f'{x:.2f}%'))

# Key test: hold F_MISSING < 0.01 (essentially no missing); does mixed-bubble still bias?
print('\n=== Key test: F_MISSING < 0.01 (effectively all founders called) ===')
clean = m[m['F_MISSING'] < 0.01]
print(f'  n (F_MISSING<0.01): {len(clean):,}')
for c in ['A_pure_snp','B_snp_only_bubble','C_mixed_bubble']:
    sub = clean[clean['cat']==c]
    if len(sub)==0:
        print(f'  {c}: n=0'); continue
    print(f'  {c}: n={len(sub):,} | mean(hf-cem)={sub["resid"].mean():+.4f} | MAE={sub["resid"].abs().mean():.4f} | %hf>cem+0.05={(sub["resid"]>0.05).mean()*100:.2f}% | %cem>hf+0.05={(sub["resid"]<-0.05).mean()*100:.2f}%')

# Same for nested LV >= 1
print('\n=== Stratify by LV (snarl-tree level) at the SNP record ===')
for lv in sorted(m['LV'].dropna().unique()):
    sub = m[m['LV']==lv]
    print(f'  LV={int(lv)}: n={len(sub):,} | mean(hf-cem)={sub["resid"].mean():+.4f} | MAE={sub["resid"].abs().mean():.4f} | %hf>cem+0.05={(sub["resid"]>0.05).mean()*100:.2f}%')

# n_recs_at_pos x F_MISSING — finer than category
print('\n=== n_recs_at_pos x F_MISSING<0.01 (mean residual, n) ===')
sub = m[m['F_MISSING']<0.01]
g = sub.groupby('n_recs_at_pos')['resid'].agg(['mean','count','std'])
g['mae'] = sub.groupby('n_recs_at_pos')['resid'].apply(lambda x: x.abs().mean())
print(g.to_string())

print(f'\ndone in {time.time()-t0:.0f}s')
