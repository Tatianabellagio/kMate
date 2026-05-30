#!/bin/bash
#SBATCH --job-name=chr1_atomcmp
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH --time=00:30:00
#SBATCH --output=logs/D2_atomcmp_%j.out
#SBATCH --error=logs/D2_atomcmp_%j.err
mkdir -p logs
set -euo pipefail

# Re-project SEEDMIX_S1 h_v3 through the ATOMIZED var_pa, then compare to hapFIRE
# on the 4-tuple (chrom, pos, ref, alt) join. Quantify encoding-disagreement reduction.

cd /global/scratch/users/tbellg/kmate/panel/arch3/chr1
PY=/global/home/users/tbellg/miniforge3/envs/hapfm/bin/python

H_PATH=/global/scratch/users/tbellg/kmate/scratch/v3qc_v3_mixedloose_chr1/SEEDMIX_S1.h_per_chrom.npz
CN_VAR=var_pa_231_arch3_chr1_atomized.var_pa.npz
VAR_CALLED=var_pa_231_arch3_chr1_atomized.var_called.npz
CN_VAR_META=var_pa_231_arch3_chr1_atomized.meta.npz
HAPFIRE=/global/scratch/projects/fc_moilab/projects/grenenet-phase1/frequency/hapFIRE_frequencies/seed_mix/s1_snp_frequency.txt
HAPFIRE_VCF=/global/scratch/projects/fc_moilab/projects/grenenet-phase1/vcf/greneNet_final_v1.1.recode.vcf

for f in $H_PATH $CN_VAR $VAR_CALLED $CN_VAR_META $HAPFIRE $HAPFIRE_VCF; do
  [ -s "$f" ] || { echo "ERROR: missing $f"; exit 1; }
done

OUT_TSV=SEEDMIX_S1_arch3_chr1_atomized.tsv
OUT_CMP=SEEDMIX_S1_chr1_atomized_AF_compare.tsv
OUT_SUMMARY=SEEDMIX_S1_chr1_atomized_AF_compare_summary.txt
HAPFIRE_REFALT=hapfire_chr1_refalt.tsv

if [ ! -s "$HAPFIRE_REFALT" ]; then
  echo "[$(date)] Extracting REF/ALT from hapFIRE panel VCF (chr1 only) ..."
  awk -F'\t' 'BEGIN{OFS="\t"; print "chrom","pos","ref","alt"} !/^#/ && $1==1 {print $1,$2,$4,$5}' "$HAPFIRE_VCF" > "$HAPFIRE_REFALT"
  wc -l "$HAPFIRE_REFALT"
fi

$PY -u <<PYEOF | tee $OUT_SUMMARY
import numpy as np
import pandas as pd
from scipy.sparse import load_npz
import time

t0 = time.time()

print('=== Load h vector (v3qc_v3) ===')
h_data = np.load('$H_PATH', allow_pickle=True)
h_raw = h_data['Chr1']
h_founders = h_data['founders']
print(f'  h shape: {h_raw.shape}, sum: {h_raw.sum():.4f}')

print('=== Load atomized var_pa + meta ===')
var_pa = load_npz('$CN_VAR').tocsr()
var_called = load_npz('$VAR_CALLED').tocsr()
meta = np.load('$CN_VAR_META', allow_pickle=True)
chrom = meta['chrom']; pos = meta['pos']; ref = meta['ref']; alt = meta['alt']
cn_founders = meta['founders']
F = var_pa.shape[0]
N = var_pa.shape[1]
print(f'  var_pa: {var_pa.shape}, {var_pa.nnz:,} nnz')
print(f'  var_called: {var_called.shape}, {var_called.nnz:,} nnz')
print(f'  meta records: {N:,}')

# Verify founder ordering matches
h_founders_s = [str(x) for x in h_founders]
cn_founders_s = [str(x) for x in cn_founders]
assert h_founders_s == cn_founders_s, 'Founder ordering differs between h and atomized var_pa'
h = h_raw.astype(float)
print('  founder ordering identical, no reindex needed')

print(f'\n=== Project AF = h @ var_pa / h @ var_called ===')
numer = h @ var_pa
denom = h @ var_called
af = np.where(denom > 0, numer / denom, np.nan)
print(f'  projected in {time.time()-t0:.1f}s; AF range [{np.nanmin(af):.3f}, {np.nanmax(af):.3f}], NaN={np.isnan(af).sum():,}')

# Write per-record TSV
print('=== Write atomized per-record TSV ===')
with open('$OUT_TSV', 'w') as fh:
    fh.write('chrom\tpos\tref\talt\talt_freq\n')
    for i in range(N):
        fh.write(f'{chrom[i]}\t{pos[i]}\t{ref[i]}\t{alt[i]}\t{af[i]:.5f}\n')
print(f'  wrote $OUT_TSV  ({N:,} rows)')

# === Compare to hapFIRE (4-tuple join) ===
print(f'\n=== Load hapFIRE AF + REF/ALT ===')
hf = pd.read_csv('$HAPFIRE', sep='\t', header=None, names=['chrom_num','pos','af_hapfire'])
hf = hf[hf.chrom_num == 1].copy()
ra = pd.read_csv('$HAPFIRE_REFALT', sep='\t')
ra = ra.rename(columns={'chrom':'chrom_num'})
ra['chrom_num'] = ra['chrom_num'].astype(int)
hf_full = hf.merge(ra, on=['chrom_num','pos'], how='left').dropna(subset=['ref','alt'])
print(f'  hapFIRE w/ REF/ALT: {len(hf_full):,}')

print('=== Build atomized DataFrame ===')
# chrom in meta is like 'Chr1' — strip and convert
chrom_num = np.array([int(str(c).replace('Chr','')) for c in chrom])
ref_str = np.array([str(r) for r in ref])
alt_str = np.array([str(a) for a in alt])
atom_df = pd.DataFrame({
    'chrom_num': chrom_num,
    'pos': pos,
    'ref': ref_str,
    'alt': alt_str,
    'af_atomized': af,
})

print('=== 4-tuple join: atomized ∩ hapFIRE ===')
m = atom_df.merge(hf_full[['chrom_num','pos','ref','alt','af_hapfire']],
                  on=['chrom_num','pos','ref','alt'], how='inner')
print(f'  4-tuple intersect: {len(m):,}  (of {len(hf_full):,} hapFIRE SNPs = {100*len(m)/len(hf_full):.1f}%)')

m = m.dropna(subset=['af_atomized'])
print(f'  finite-AF after drop NaN: {len(m):,}')

d = m['af_atomized'] - m['af_hapfire']
mae = np.abs(d).mean()
rmse = np.sqrt((d**2).mean())
bias = d.mean()
r = np.corrcoef(m['af_atomized'], m['af_hapfire'])[0,1]
print(f'\n=== ATOMIZED vs hapFIRE ===')
print(f'  MAE={mae:.4f}  RMSE={rmse:.4f}  bias={bias:+.4f}  r={r:.4f}')
for thr in [0.01, 0.05, 0.10]:
    n_out = (np.abs(d) > thr).sum()
    print(f'    |Δ|>{thr}: {n_out:,} ({100*n_out/len(m):.2f}%)')

# Compare to RAW (path-aware) on the same hapFIRE positions
print(f'\n=== Compare to RAW (path-aware) Arch 3 from prior A7 run ===')
raw_cmp = pd.read_csv('SEEDMIX_S1_chr1_AF_compare.tsv', sep='\t')
print(f'  raw 4-tuple rows: {len(raw_cmp):,}')

# Join atomized result with raw on (chrom_num, pos, ref, alt) to compare side-by-side
side = m[['chrom_num','pos','ref','alt','af_atomized','af_hapfire']].merge(
    raw_cmp[['chrom_num','pos','ref','alt','af_new','af_old']].rename(
        columns={'af_new':'af_arch3_raw','af_old':'af_v3qc_raw'}),
    on=['chrom_num','pos','ref','alt'], how='inner'
)
print(f'  3-way join (atomized ∩ arch3_raw ∩ hapFIRE): {len(side):,}')

d_raw = side['af_arch3_raw'] - side['af_hapfire']
d_atom = side['af_atomized'] - side['af_hapfire']
d_v3 = side['af_v3qc_raw'] - side['af_hapfire']

def stats(name, d):
    print(f'  {name:<20}  MAE={np.abs(d).mean():.4f}  RMSE={np.sqrt((d**2).mean()):.4f}  bias={d.mean():+.4f}  |Δ|>0.10: {(np.abs(d)>0.10).sum():,} ({100*(np.abs(d)>0.10).sum()/len(side):.2f}%)')

print(f'\n=== Head-to-head on {len(side):,} shared 4-tuple SNPs ===')
stats('v3qc_v3 raw vs hF',  d_v3)
stats('arch3 raw vs hF',    d_raw)
stats('arch3 ATOMIZED vs hF', d_atom)

# Which raw outliers got fixed by atomization?
raw_outl = (np.abs(d_raw) > 0.10)
atom_outl = (np.abs(d_atom) > 0.10)
fixed = raw_outl & ~atom_outl
introduced = ~raw_outl & atom_outl
print(f'\n=== Outlier transitions (|Δ|>0.10) ===')
print(f'  raw outliers fixed by atomization:        {fixed.sum():,}  ({100*fixed.sum()/raw_outl.sum():.1f}% of raw outliers)')
print(f'  outliers introduced by atomization:       {introduced.sum():,}')
print(f'  raw outliers still outlier post-atomize:  {(raw_outl & atom_outl).sum():,}')

# Spot-check: positions where raw was outlier and atomized fixed it
print(f'\n=== Sample of FIXED outliers (top 10 by |Δ_raw|) ===')
fix_idx = np.where(fixed)[0]
fix_top = side.iloc[fix_idx].assign(d_raw=d_raw.iloc[fix_idx], d_atom=d_atom.iloc[fix_idx]).copy()
fix_top['abs_d_raw'] = np.abs(fix_top['d_raw'])
top = fix_top.sort_values('abs_d_raw', ascending=False).head(10)
print(top[['chrom_num','pos','ref','alt','af_hapfire','af_arch3_raw','af_atomized','d_raw','d_atom']].to_string(index=False))

# Save the 3-way comparison
side.to_csv('$OUT_CMP', sep='\t', index=False)
print(f'\nWrote $OUT_CMP ({len(side):,} rows)')
print(f'Total time: {time.time()-t0:.1f}s')
PYEOF

echo
ls -lh $OUT_TSV $OUT_CMP $OUT_SUMMARY
echo "[$(date)] DONE D2"
