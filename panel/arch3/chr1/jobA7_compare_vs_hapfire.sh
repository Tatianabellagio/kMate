#!/bin/bash
#SBATCH --job-name=chr1_a7_cmp
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=logs/A7_cmp_%j.out
#SBATCH --error=logs/A7_cmp_%j.err
mkdir -p logs
set -euo pipefail

# Phase 2 Job A7: compare per-SNP AF on Chr1
#   NEW = SEEDMIX_S1_arch3_chr1.tsv (A6 output: h_v3 @ var_pa_arch3)
#   OLD = /scratch/v3qc_v3_mixedloose_chr1/SEEDMIX_S1.tsv (already-computed v3 panel projection)
#   HAPFIRE = xwu's production hapFIRE on the 231-panel (greneNet_final_v1.1.recode.vcf).
#         This is what SEEDMIX_S1_v3_vs_hapfire.ipynb consumes. The contamination_test folder
#         is a SEPARATE 1141-panel experiment for non-231 contamination detection — NOT the
#         AF reference.
#   HAPFIRE_VCF = the panel VCF hapFIRE was run on, used to pull (REF, ALT) so we can
#         join on the 4-tuple (chrom, pos, ref, alt) instead of just (chrom, pos).
#         CRITICAL: multi-allelic split records in NEW share a `pos` with other rows;
#         joining on pos only manufactures off-diagonal scatter (per notebook L268).
# Per memory feedback_use_mae_not_r2.md: lead with MAE.

cd /global/scratch/users/tbellg/kmate/panel/arch3/chr1
PY=/global/home/users/tbellg/miniforge3/envs/hapfm/bin/python

NEW=SEEDMIX_S1_arch3_chr1.tsv
OLD=/global/scratch/users/tbellg/kmate/scratch/v3qc_v3_mixedloose_chr1/SEEDMIX_S1.tsv
HAPFIRE=/global/scratch/projects/fc_moilab/projects/grenenet-phase1/frequency/hapFIRE_frequencies/seed_mix/s1_snp_frequency.txt
HAPFIRE_VCF=/global/scratch/projects/fc_moilab/projects/grenenet-phase1/vcf/greneNet_final_v1.1.recode.vcf

[ -s "$NEW" ] || { echo "ERROR: missing $NEW (A6 not done)"; exit 1; }
[ -s "$OLD" ] || { echo "ERROR: missing $OLD"; exit 1; }
[ -s "$HAPFIRE" ] || { echo "ERROR: missing $HAPFIRE"; exit 1; }
[ -s "$HAPFIRE_VCF" ] || { echo "ERROR: missing $HAPFIRE_VCF"; exit 1; }

OUT_JOIN=SEEDMIX_S1_chr1_AF_compare.tsv
OUT_SUMMARY=SEEDMIX_S1_chr1_AF_compare_summary.txt
HAPFIRE_REFALT=hapfire_chr1_refalt.tsv

# Extract (chrom, pos, ref, alt) from greneNet_final_v1.1 chr1 once, for 4-tuple join.
echo "[$(date)] Extracting REF/ALT from hapFIRE panel VCF (chr1 only) ..."
awk -F'\t' 'BEGIN{OFS="\t"; print "chrom","pos","ref","alt"} !/^#/ && $1==1 {print $1,$2,$4,$5}' "$HAPFIRE_VCF" > "$HAPFIRE_REFALT"
wc -l "$HAPFIRE_REFALT"

$PY -u <<PYEOF | tee $OUT_SUMMARY
import pandas as pd
import numpy as np

# Per SEEDMIX_S1_v3_vs_hapfire.ipynb L268-285, L620-622:
#   Joining NEW vs hapFIRE on (chrom,pos) only joins each multi-allelic-split row
#   against the same hapFIRE row → manufactured off-diagonal scatter. Fix: 4-tuple join.

print('=== Load NEW Arch 3 (Chr1, SNPs only) ===')
new = pd.read_csv('$NEW', sep='\t')
# strip "Chr" prefix if present so chrom matches hapFIRE (which uses "1" not "Chr1")
new['chrom_num'] = new['chrom'].astype(str).str.replace('Chr','',regex=False).astype(int)
print(f'  total NEW records: {len(new):,}')
new_snp = new[(new.ref_len==1) & (new.alt_len==1)].copy()
print(f'  SNP-only: {len(new_snp):,}')

print()
print('=== Load OLD v3qc_v3 (Chr1, SNPs only) ===')
old = pd.read_csv('$OLD', sep='\t')
old['chrom_num'] = old['chrom'].astype(str).str.replace('Chr','',regex=False).astype(int)
print(f'  total OLD records: {len(old):,}')
old_snp = old[(old.ref_len==1) & (old.alt_len==1)].copy()
print(f'  SNP-only: {len(old_snp):,}')

# NEW + OLD: do we have ref/alt columns? var_pa meta has them, but A6 only writes ref_len/alt_len.
# Need to re-load from meta.npz to attach actual REF/ALT bases.
# CRITICAL: meta['ref']/meta['alt'] are object arrays — DO NOT cast to fixed-width strings
# via np.array([str(x) for x in ref]), some cactus SV REF alleles are 94kb+ → 925 GiB allocation.
print()
print('=== Re-attach REF/ALT bases from meta.npz to NEW ===')
meta = np.load('var_pa_231_arch3_chr1.meta.npz', allow_pickle=True)
# Subset to SNPs FIRST using ref_len/alt_len (which are int arrays), then convert just the
# subset to a list of strings — keeps the 925 GiB monster off the heap.
meta_ref_len = meta['ref_len']; meta_alt_len = meta['alt_len']
snp_mask = (meta_ref_len == 1) & (meta_alt_len == 1)
print(f'  meta records: {len(meta_ref_len):,}, SNPs (ref_len=alt_len=1): {snp_mask.sum():,}')
# A6 writes records in the same order as meta, so we can attach by index.
new_full = pd.read_csv('$NEW', sep='\t')
assert len(new_full) == len(meta['pos']), (
    f"row count mismatch: NEW TSV {len(new_full):,} != meta {len(meta['pos']):,}"
)
# Slice to SNPs only — small enough to materialize ref/alt strings safely.
new_snp = new_full.iloc[snp_mask].copy()
new_snp['ref'] = [str(x) for x in np.asarray(meta['ref'])[snp_mask]]
new_snp['alt'] = [str(x) for x in np.asarray(meta['alt'])[snp_mask]]
new_snp['chrom_num'] = new_snp['chrom'].astype(str).str.replace('Chr','',regex=False).astype(int)
print(f'  NEW SNP records w/ REF/ALT: {len(new_snp):,}')

# OLD: meta is at /data/var_pa_231_v3qc_v3.meta.npz
print('=== Re-attach REF/ALT bases to OLD ===')
old_meta_path = '/global/scratch/users/tbellg/kmate/data/var_pa_231_v3qc_v3.meta.npz'
try:
    om = np.load(old_meta_path, allow_pickle=True)
    om_ref_len = om['ref_len']; om_alt_len = om['alt_len']
    om_snp_mask = (om_ref_len == 1) & (om_alt_len == 1)
    if len(old) == len(om['pos']):
        # Re-slice old TSV to SNPs (was already done above but redoing on the full old to align with om order)
        old_snp = old.iloc[om_snp_mask].copy()
        old_snp['ref'] = [str(x) for x in np.asarray(om['ref'])[om_snp_mask]]
        old_snp['alt'] = [str(x) for x in np.asarray(om['alt'])[om_snp_mask]]
        old_snp['chrom_num'] = old_snp['chrom'].astype(str).str.replace('Chr','',regex=False).astype(int)
        print(f'  OLD SNP records w/ REF/ALT: {len(old_snp):,}')
    else:
        print(f'  WARN: OLD TSV {len(old):,} != OLD meta {len(om["pos"]):,} — skipping ref/alt attach for OLD')
        old_snp = old[(old.ref_len==1) & (old.alt_len==1)].copy()
        old_snp['ref'] = None; old_snp['alt'] = None
except FileNotFoundError:
    print(f'  WARN: {old_meta_path} not found — OLD comparison will fall back to (pos)-only join')
    old_snp = old[(old.ref_len==1) & (old.alt_len==1)].copy()
    old_snp['ref'] = None; old_snp['alt'] = None

print()
print('=== Load hapFIRE AF + REF/ALT from greneNet_final_v1.1 ===')
hf = pd.read_csv('$HAPFIRE', sep='\t', header=None, names=['chrom_num','pos','af_hapfire'])
hf = hf[hf.chrom_num == 1].copy()
ra = pd.read_csv('$HAPFIRE_REFALT', sep='\t')
ra = ra.rename(columns={'chrom':'chrom_num'})
ra['chrom_num'] = ra['chrom_num'].astype(int)
print(f'  hapFIRE chr1 SNPs: {len(hf):,}')
print(f'  refalt records:    {len(ra):,}')
hf_full = hf.merge(ra, on=['chrom_num','pos'], how='left')
n_missing = hf_full['ref'].isna().sum()
print(f'  hapFIRE rows missing REF/ALT after join: {n_missing:,}  (should be 0 or near-0)')
hf_full = hf_full.dropna(subset=['ref','alt']).copy()
print(f'  hapFIRE w/ REF/ALT: {len(hf_full):,}')

print()
print('=== Join NEW ∩ hapFIRE on (chrom, pos, ref, alt) ===')
m_new = new_snp.merge(
    hf_full[['chrom_num','pos','ref','alt','af_hapfire']],
    on=['chrom_num','pos','ref','alt'], how='inner'
)
print(f'  4-tuple intersect: {len(m_new):,}  (of {len(hf_full):,} hapFIRE SNPs = {100*len(m_new)/len(hf_full):.1f}%)')
# Also report pos-only intersect for comparison with notebook
m_new_posonly = new_snp[['chrom_num','pos']].drop_duplicates().merge(
    hf_full[['chrom_num','pos']].drop_duplicates(), on=['chrom_num','pos'], how='inner'
)
print(f'  pos-only intersect (positions present in both, dedup): {len(m_new_posonly):,}')

dif_new = m_new['alt_freq'] - m_new['af_hapfire']
mae_new = np.abs(dif_new).mean()
rmse_new = np.sqrt((dif_new**2).mean())
bias_new = dif_new.mean()
print(f'  NEW vs hapFIRE: MAE={mae_new:.4f}  RMSE={rmse_new:.4f}  bias={bias_new:+.4f}')
for thr in [0.01, 0.05, 0.10]:
    n_out = (np.abs(dif_new) > thr).sum()
    print(f'    |Δ|>{thr}: {n_out:,} ({100*n_out/len(m_new):.2f}%)')

print()
print('=== Join OLD ∩ hapFIRE on (chrom, pos, ref, alt) if possible ===')
if 'ref' in old_snp.columns and old_snp['ref'].notna().any():
    m_old = old_snp.merge(
        hf_full[['chrom_num','pos','ref','alt','af_hapfire']],
        on=['chrom_num','pos','ref','alt'], how='inner'
    )
    print(f'  4-tuple intersect: {len(m_old):,}  ({100*len(m_old)/len(hf_full):.1f}% of hapFIRE)')
    dif_old = m_old['alt_freq'] - m_old['af_hapfire']
    mae_old = np.abs(dif_old).mean()
    rmse_old = np.sqrt((dif_old**2).mean())
    bias_old = dif_old.mean()
    print(f'  OLD vs hapFIRE: MAE={mae_old:.4f}  RMSE={rmse_old:.4f}  bias={bias_old:+.4f}')
    for thr in [0.01, 0.05, 0.10]:
        n_out = (np.abs(dif_old) > thr).sum()
        print(f'    |Δ|>{thr}: {n_out:,} ({100*n_out/len(m_old):.2f}%)')
else:
    print('  skipped (no OLD ref/alt)')
    m_old = None

print()
print('=== NEW ∩ OLD ∩ hapFIRE (3-way, 4-tuple) ===')
if m_old is not None and len(m_old) > 0:
    a = m_new[['chrom_num','pos','ref','alt','alt_freq','af_hapfire']].rename(columns={'alt_freq':'af_new'})
    b = m_old[['chrom_num','pos','ref','alt','alt_freq']].rename(columns={'alt_freq':'af_old'})
    three = a.merge(b, on=['chrom_num','pos','ref','alt'], how='inner')
    print(f'  3-way intersect rows: {len(three):,}')
    if len(three) > 0:
        d_new = three['af_new'] - three['af_hapfire']
        d_old = three['af_old'] - three['af_hapfire']
        d_no = three['af_new'] - three['af_old']
        print(f'  NEW vs hapFIRE (3way): MAE={np.abs(d_new).mean():.4f}  bias={d_new.mean():+.4f}  |Δ|>0.10: {(np.abs(d_new)>0.10).sum():,}')
        print(f'  OLD vs hapFIRE (3way): MAE={np.abs(d_old).mean():.4f}  bias={d_old.mean():+.4f}  |Δ|>0.10: {(np.abs(d_old)>0.10).sum():,}')
        print(f'  NEW vs OLD     (3way): MAE={np.abs(d_no).mean():.4f}  bias={d_no.mean():+.4f}  |Δ|>0.10: {(np.abs(d_no)>0.10).sum():,}')
        delta = np.abs(d_new).mean() - np.abs(d_old).mean()
        sign = 'BETTER' if delta < 0 else 'WORSE'
        print(f'  → Arch 3 is {abs(100*delta/np.abs(d_old).mean()):.1f}% {sign} than v3qc_v3 (MAE delta vs hapFIRE: {delta:+.4f})')
        three.to_csv('$OUT_JOIN', sep='\t', index=False)
        print(f'  wrote $OUT_JOIN ({len(three):,} rows)')

print()
print('=== Spot-checks at 3 canonical positions (NEW vs hapFIRE) ===')
for tp, tr, ta in [(5870018,'T','A'), (10421645,'T','C'), (13843898,'C','T')]:
    sub = m_new[(m_new.pos==tp) & (m_new.ref==tr) & (m_new.alt==ta)]
    if len(sub) == 0:
        # try without ref/alt restriction
        sub2 = m_new[m_new.pos==tp]
        if len(sub2) == 0:
            print(f'  Chr1:{tp} {tr}>{ta}: not in NEW∩hapFIRE 4-tuple, and pos absent')
        else:
            for _, r in sub2.iterrows():
                print(f'  Chr1:{tp} {tr}>{ta}: only different ALT present — {r.ref}>{r.alt} NEW={r.alt_freq:.4f} hf={r.af_hapfire:.4f}')
    else:
        for _, r in sub.iterrows():
            print(f'  Chr1:{tp} {tr}>{ta}: NEW={r.alt_freq:.4f}  hapFIRE={r.af_hapfire:.4f}  Δ={r.alt_freq-r.af_hapfire:+.4f}')
PYEOF

echo
ls -lh $OUT_JOIN $OUT_SUMMARY
echo "[$(date)] DONE A7"
