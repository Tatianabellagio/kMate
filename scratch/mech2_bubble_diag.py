"""
Mechanism 2 diagnostic: are high (hapFIRE - cactus_em) residuals at SNPs
enriched at positions where the cactus pangenome bubble contains non-SNP alleles?

Inputs:
- cactus_em AF: scratch/v3qc_v3_mixedloose_chr1/SEEDMIX_S1.tsv (1:1 with cn_var meta rows)
- cn_var meta: poolfreq/data/cn_var_231_v3qc_v3.meta.npz (chrom, pos, ref, alt, ref_len, alt_len)
- hapfire SNP AF: /carnegie/.../GrENE_net/hapFIRE_updatedVCF/s1_1_density0.6_snp_frequency.txt (Chr1)
- greneNet VCF: gives REF/ALT for hapfire SNPs
- merged panel VCF: founders_231_v3qc_v3.haploid.vcf.gz -> for LV + F_MISSING per record

Plan:
  1. Load cactus_em + meta. Restrict to Chr1, ACGT->ACGT SNPs.
  2. Compute multi-allelicity signature: n_records_at_pos, max_alt_len_at_pos,
     n_nonSNP_records_at_pos, has_large_alt_at_pos (>50bp).
  3. Load greneNet VCF (Chr1, REF/ALT) + hapfire SNP AF -> hapfire SNP AF keyed by (pos, REF, ALT).
  4. Inner join cactus_em SNPs ⨝ hapfire SNPs on (pos, REF, ALT).
  5. Compute residual = hapfire_AF - cactus_em_AF.
  6. Stratify by bubble context. Confound check: do the same stratified by F_MISSING.

Outputs printed to stdout; nothing written to disk except a small summary.
"""
import sys, time
from pathlib import Path
import numpy as np
import pandas as pd
from collections import defaultdict

ROOT = Path('/carnegie/nobackup/scratch/tbellagio/hapfire_sv')

t0 = time.time()
print('=== loading cactus_em + meta ===', flush=True)
meta = np.load(ROOT/'poolfreq/data/cn_var_231_v3qc_v3.meta.npz', allow_pickle=True)
chrom_arr  = meta['chrom']
pos_arr    = meta['pos']
ref_arr    = meta['ref']
alt_arr    = meta['alt']
ref_len    = meta['ref_len']
alt_len    = meta['alt_len']
N_meta = len(pos_arr)
print(f'  meta records: {N_meta:,}')

cem = pd.read_csv(ROOT/'scratch/v3qc_v3_mixedloose_chr1/SEEDMIX_S1.tsv', sep='\t')
print(f'  cactus_em records: {len(cem):,}')
assert len(cem) == N_meta, 'meta vs cactus_em row count mismatch'
# Sanity check alignment by chrom/pos/ref_len/alt_len
assert (cem['chrom'].values == chrom_arr).all(), 'chrom mismatch'
assert (cem['pos'].values == pos_arr).all(), 'pos mismatch'
assert (cem['ref_len'].values == ref_len).all(), 'ref_len mismatch'
assert (cem['alt_len'].values == alt_len).all(), 'alt_len mismatch'
print('  alignment OK')

# Restrict to Chr1 first to avoid blowing memory on object-dtype REF/ALT strings
chr1_mask = (chrom_arr == 'Chr1')
print(f'  Chr1 records (raw): {chr1_mask.sum():,}')
df1 = pd.DataFrame({
    'pos':      pos_arr[chr1_mask],
    'ref_len':  ref_len[chr1_mask],
    'alt_len':  alt_len[chr1_mask],
    'cem_af':   cem.loc[chr1_mask, 'alt_freq'].values,
})
# Only materialize REF/ALT for SNP-shaped rows (cheap; the SV strings live untouched in the npz)
chr1_idx = np.where(chr1_mask)[0]
df1['ref_idx'] = chr1_idx  # keep pointer for later
df1['chrom'] = 'Chr1'

# Position-level bubble signatures (over the FULL Chr1 record set, before SNP filter)
print('\n=== computing bubble-context signatures (Chr1) ===', flush=True)
pos_group = df1.groupby('pos')
n_recs_at_pos       = pos_group['pos'].transform('count').values
max_alt_at_pos      = pos_group['alt_len'].transform('max').values
max_ref_at_pos      = pos_group['ref_len'].transform('max').values
# count of records at this pos that are NOT SNP-shaped (ref_len>1 or alt_len>1)
df1['_is_nonsnp_local'] = ((df1['ref_len'] > 1) | (df1['alt_len'] > 1)).astype(int)
n_nonsnp_at_pos = df1.groupby('pos')['_is_nonsnp_local'].transform('sum').values
df1.drop(columns='_is_nonsnp_local', inplace=True)

df1['n_recs_at_pos']    = n_recs_at_pos
df1['max_alt_at_pos']   = max_alt_at_pos
df1['max_ref_at_pos']   = max_ref_at_pos
df1['n_nonsnp_at_pos']  = n_nonsnp_at_pos
df1['has_sv_at_pos']    = (np.maximum(max_alt_at_pos, max_ref_at_pos) >= 50).astype(int)

# Filter to SNP-shaped, then materialize REF/ALT just for those rows (cheap, length-1)
is_snp = (df1['ref_len']==1) & (df1['alt_len']==1)
snp = df1[is_snp].copy().reset_index(drop=True)
snp_idx = snp['ref_idx'].values
# Pull bases as plain str (1 char each)
snp_ref = np.array([str(ref_arr[i]) for i in snp_idx], dtype='<U1')
snp_alt = np.array([str(alt_arr[i]) for i in snp_idx], dtype='<U1')
ACGT = np.array(['A','C','G','T'])
clean = np.isin(snp_ref, ACGT) & np.isin(snp_alt, ACGT)
snp = snp[clean].copy().reset_index(drop=True)
snp['ref'] = snp_ref[clean]
snp['alt'] = snp_alt[clean]
print(f'  SNP-shaped (ACGT) Chr1 records: {len(snp):,}')

# Distribution of bubble context among Chr1 SNPs
print('\n--- Chr1 SNP records by bubble context ---')
print(f'  n_recs_at_pos==1 (pos has ONLY this SNP): {(snp["n_recs_at_pos"]==1).sum():,}')
print(f'  n_recs_at_pos>=2: {(snp["n_recs_at_pos"]>=2).sum():,}')
print(f'    of which has_sv_at_pos: {((snp["n_recs_at_pos"]>=2)&(snp["has_sv_at_pos"]==1)).sum():,}')
print(f'  has_sv_at_pos (>=50bp REF or ALT at same pos): {(snp["has_sv_at_pos"]==1).sum():,}')
print(f'  n_nonsnp_at_pos>=1: {(snp["n_nonsnp_at_pos"]>=1).sum():,}')

# Load greneNet VCF (Chr1) for REF/ALT lookup
print('\n=== loading greneNet VCF (Chr1) for REF/ALT keying ===', flush=True)
gv = []
with open('/global/scratch/projects/fc_moilab/projects/grenenet-phase1/vcf/greneNet_final_v1.1.recode.vcf') as fh:
    for line in fh:
        if line.startswith('#'): continue
        f = line.split('\t', 6)
        if f[0] != '1': continue
        gv.append((int(f[1]), f[3], f[4]))
gv_df = pd.DataFrame(gv, columns=['pos','ref','alt'])
print(f'  greneNet chr1 SNPs: {len(gv_df):,}')

hf = pd.read_csv(
    # TODO: hapFIRE updated VCF run output, not present on moilab mirror — re-generate or skip if no longer needed
    '/carnegie/nobackup/scratch/xwu/GrENE_net/hapFIRE_updatedVCF/s1_1_density0.6_snp_frequency.txt',
    sep='\t', header=None, names=['chr','pos','hf_af'])
hf = hf[hf['chr']==1].copy()
print(f'  hapfire chr1 records: {len(hf):,}')
hf = hf.merge(gv_df, on='pos', how='inner')
print(f'  after merging w/ greneNet REF/ALT: {len(hf):,}')

# Join cactus_em SNPs ⨝ hapfire on (pos, ref, alt)
print('\n=== joining cactus_em SNPs to hapfire ===', flush=True)
merged = snp.merge(hf[['pos','ref','alt','hf_af']], on=['pos','ref','alt'], how='inner')
print(f'  joined records: {len(merged):,}')
merged['resid'] = merged['hf_af'] - merged['cem_af']
merged['abs_resid'] = merged['resid'].abs()

print('\n--- overall residual stats (hapfire - cactus_em) ---')
print(merged['resid'].describe())
print(f'MAE = {merged["abs_resid"].mean():.5f}')
print(f'mean residual (hf - cem) = {merged["resid"].mean():.5f}')
print(f'|resid| > 0.05: {(merged["abs_resid"]>0.05).sum():,} ({100*(merged["abs_resid"]>0.05).mean():.3f}%)')
print(f'resid >  0.05 (hf > cem): {(merged["resid"]> 0.05).sum():,}')
print(f'resid < -0.05 (cem > hf): {(merged["resid"]<-0.05).sum():,}')

# === Stratify by bubble context ===
def show_strat(label, mask):
    sub = merged[mask]
    n = len(sub)
    if n == 0:
        print(f'  {label}: n=0')
        return
    print(f'  {label}: n={n:,} | MAE={sub["abs_resid"].mean():.4f} | '
          f'mean(hf-cem)={sub["resid"].mean():+.4f} | '
          f'%|res|>0.05={100*(sub["abs_resid"]>0.05).mean():.2f}% | '
          f'%hf>cem+0.05={100*(sub["resid"]>0.05).mean():.2f}% | '
          f'%cem>hf+0.05={100*(sub["resid"]<-0.05).mean():.2f}%')

print('\n--- A. Stratify by n_recs_at_pos (was the SNP decomposed from a multi-allelic bubble?) ---')
show_strat('n_recs=1 (pure SNP, no multi-allelic decomp)', merged['n_recs_at_pos']==1)
show_strat('n_recs=2',                                     merged['n_recs_at_pos']==2)
show_strat('n_recs=3-5',                                  (merged['n_recs_at_pos']>=3)&(merged['n_recs_at_pos']<=5))
show_strat('n_recs>=6',                                    merged['n_recs_at_pos']>=6)

print('\n--- B. Stratify by has_sv_at_pos (>=50bp REF or ALT exists at same position) ---')
show_strat('has_sv_at_pos=0', merged['has_sv_at_pos']==0)
show_strat('has_sv_at_pos=1', merged['has_sv_at_pos']==1)

print('\n--- C. Stratify by n_nonsnp_at_pos ---')
show_strat('n_nonsnp=0 (pos contains ONLY SNP records)', merged['n_nonsnp_at_pos']==0)
show_strat('n_nonsnp>=1',                                merged['n_nonsnp_at_pos']>=1)

# Composite categories
print('\n--- D. Composite (no other record at pos / shares pos with SNP only / shares pos with non-SNP) ---')
cat = np.where(merged['n_recs_at_pos']==1, 'A_pure_snp',
       np.where(merged['n_nonsnp_at_pos']==0, 'B_snp_only_bubble', 'C_mixed_bubble'))
merged['cat'] = cat
for c in ['A_pure_snp','B_snp_only_bubble','C_mixed_bubble']:
    show_strat(c, merged['cat']==c)

# Enrichment among outliers
out = merged[merged['abs_resid']>0.05]
all_n = len(merged); out_n = len(out)
print(f'\n--- E. Enrichment of bubble context among |residual|>0.05 outliers (n_out={out_n:,}) ---')
for c in ['A_pure_snp','B_snp_only_bubble','C_mixed_bubble']:
    f_all = (merged['cat']==c).mean()
    f_out = (out['cat']==c).mean() if out_n else 0
    print(f'  {c}: %all={100*f_all:.2f}%, %outliers={100*f_out:.2f}%, '
          f'enrichment={ (f_out/f_all if f_all>0 else float("nan")):.2f}x')

# Sign-specific enrichment: hapfire > cactus_em (the predicted direction)
out_hf = merged[merged['resid']>0.05]
print(f'\n--- F. Same, restricted to hf > cem outliers (direction predicted by Mechanism 2; n={len(out_hf):,}) ---')
for c in ['A_pure_snp','B_snp_only_bubble','C_mixed_bubble']:
    f_all = (merged['cat']==c).mean()
    f_out = (out_hf['cat']==c).mean() if len(out_hf) else 0
    print(f'  {c}: %all={100*f_all:.2f}%, %hf>cem outliers={100*f_out:.2f}%, '
          f'enrichment={ (f_out/f_all if f_all>0 else float("nan")):.2f}x')

# Sign-specific enrichment: cactus_em > hapfire
out_cm = merged[merged['resid']<-0.05]
print(f'\n--- G. Same, restricted to cem > hf outliers (opposite direction; n={len(out_cm):,}) ---')
for c in ['A_pure_snp','B_snp_only_bubble','C_mixed_bubble']:
    f_all = (merged['cat']==c).mean()
    f_out = (out_cm['cat']==c).mean() if len(out_cm) else 0
    print(f'  {c}: %all={100*f_all:.2f}%, %cem>hf outliers={100*f_out:.2f}%, '
          f'enrichment={ (f_out/f_all if f_all>0 else float("nan")):.2f}x')

# Save a per-record dump for downstream notebook use
out_tsv = ROOT/'scratch/v3qc_v3_mixedloose_chr1/mech2_chr1_joined.tsv.gz'
keep = merged[['chrom','pos','ref','alt','cem_af','hf_af','resid',
               'n_recs_at_pos','max_alt_at_pos','max_ref_at_pos','n_nonsnp_at_pos','has_sv_at_pos','cat']]
keep.to_csv(out_tsv, sep='\t', index=False, compression='gzip')
print(f'\nwrote {out_tsv} ({len(keep):,} rows)')

print(f'\ndone in {time.time()-t0:.0f}s')
