#!/bin/bash
#SBATCH --job-name=chr1_filt2cmp
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH --time=00:30:00
#SBATCH --output=logs/F1_filt2cmp_%j.out
#SBATCH --error=logs/F1_filt2cmp_%j.err
mkdir -p logs
set -euo pipefail

# Project SEEDMIX_S1 h_filt2 through arch3 atomized var_pa, compare to hapFIRE,
# and do a 3-way head-to-head: mixedloose vs filt2 vs hapFIRE.

cd /global/scratch/users/tbellg/kmate/panel/arch3/chr1
PY=/global/home/users/tbellg/miniforge3/envs/hapfm/bin/python

H_FILT2=/global/scratch/users/tbellg/kmate/scratch/v3qc_v3_filt2_chr1/SEEDMIX_S1.h_per_chrom.npz
H_OLD=/global/scratch/users/tbellg/kmate/scratch/v3qc_v3_mixedloose_chr1/SEEDMIX_S1.h_per_chrom.npz
CN_VAR=var_pa_231_arch3_chr1_atomized.var_pa.npz
VAR_CALLED=var_pa_231_arch3_chr1_atomized.var_called.npz
CN_VAR_META=var_pa_231_arch3_chr1_atomized.meta.npz
HAPFIRE=/global/scratch/projects/fc_moilab/projects/grenenet-phase1/frequency/hapFIRE_frequencies/seed_mix/s1_snp_frequency.txt
HAPFIRE_REFALT=hapfire_chr1_refalt.tsv   # built by D2

for f in $H_FILT2 $H_OLD $CN_VAR $VAR_CALLED $CN_VAR_META $HAPFIRE $HAPFIRE_REFALT; do
  [ -s "$f" ] || { echo "ERROR: missing $f"; exit 1; }
done

$PY -u <<'PYEOF'
import numpy as np
import pandas as pd
from scipy.sparse import load_npz
import time

t0 = time.time()

H_FILT2  = '/global/scratch/users/tbellg/kmate/scratch/v3qc_v3_filt2_chr1/SEEDMIX_S1.h_per_chrom.npz'
H_OLD    = '/global/scratch/users/tbellg/kmate/scratch/v3qc_v3_mixedloose_chr1/SEEDMIX_S1.h_per_chrom.npz'
HAPFIRE  = '/global/scratch/projects/fc_moilab/projects/grenenet-phase1/frequency/hapFIRE_frequencies/seed_mix/s1_snp_frequency.txt'
HAPFIRE_REFALT = 'hapfire_chr1_refalt.tsv'

print('=== Load var_pa (atomized) ===')
var_pa        = load_npz('var_pa_231_arch3_chr1_atomized.var_pa.npz').tocsr()
var_called = load_npz('var_pa_231_arch3_chr1_atomized.var_called.npz').tocsr()
meta          = np.load('var_pa_231_arch3_chr1_atomized.meta.npz', allow_pickle=True)
chrom = meta['chrom']; pos = meta['pos']; ref = meta['ref']; alt = meta['alt']
cn_founders = [str(x) for x in meta['founders']]
N = var_pa.shape[1]
print(f'  var_pa: {var_pa.shape}, {var_pa.nnz:,} nnz')

def load_h(path, label):
    d = np.load(path, allow_pickle=True)
    h = d['Chr1'].astype(np.float64)
    founders = [str(x) for x in d['founders']]
    assert founders == cn_founders, f'{label}: founder mismatch vs var_pa'
    return h

h_filt2 = load_h(H_FILT2, 'filt2')
h_old   = load_h(H_OLD,   'mixedloose')

def project(h, var_pa, var_called):
    numer = h @ var_pa
    denom = h @ var_called
    return np.where(denom > 0, numer / denom, np.nan)

print('=== Project both h vectors ===')
af_filt2 = project(h_filt2, var_pa, var_called)
af_old   = project(h_old,   var_pa, var_called)
print(f'  filt2:      AF range [{np.nanmin(af_filt2):.3f}, {np.nanmax(af_filt2):.3f}], NaN={np.isnan(af_filt2).sum():,}')
print(f'  mixedloose: AF range [{np.nanmin(af_old):.3f}, {np.nanmax(af_old):.3f}], NaN={np.isnan(af_old).sum():,}')
print(f'  elapsed: {time.time()-t0:.1f}s')

# Build DataFrames
chrom_num = np.array([int(str(c).replace('Chr','')) for c in chrom])
ref_str   = np.array([str(r) for r in ref])
alt_str   = np.array([str(a) for a in alt])

df = pd.DataFrame({
    'chrom_num': chrom_num,
    'pos':       pos,
    'ref':       ref_str,
    'alt':       alt_str,
    'af_filt2':  af_filt2,
    'af_old':    af_old,
})

print('\n=== Load hapFIRE + REF/ALT ===')
hf = pd.read_csv(HAPFIRE, sep='\t', header=None, names=['chrom_num','pos','af_hapfire'])
hf = hf[hf.chrom_num == 1].copy()
ra = pd.read_csv(HAPFIRE_REFALT, sep='\t').rename(columns={'chrom':'chrom_num'})
ra['chrom_num'] = ra['chrom_num'].astype(int)
hf_full = hf.merge(ra, on=['chrom_num','pos'], how='left').dropna(subset=['ref','alt'])
print(f'  hapFIRE w/ REF/ALT: {len(hf_full):,}')

print('\n=== 4-tuple join ===')
m = df.merge(hf_full[['chrom_num','pos','ref','alt','af_hapfire']],
             on=['chrom_num','pos','ref','alt'], how='inner')
m = m.dropna(subset=['af_filt2','af_old'])
print(f'  joined: {len(m):,}  ({100*len(m)/len(hf_full):.1f}% of hapFIRE SNPs)')

def stats(name, pred, truth):
    d = pred - truth
    mae  = np.abs(d).mean()
    rmse = np.sqrt((d**2).mean())
    bias = d.mean()
    n10  = (np.abs(d) > 0.10).sum()
    n05  = (np.abs(d) > 0.05).sum()
    print(f'  {name:<24}  MAE={mae:.4f}  RMSE={rmse:.4f}  bias={bias:+.4f}  |Δ|>0.05:{n05:,}  |Δ|>0.10:{n10:,} ({100*n10/len(d):.2f}%)')
    return d

print(f'\n=== Head-to-head on {len(m):,} shared SNPs ===')
d_old   = stats('mixedloose vs hapFIRE', m['af_old'],   m['af_hapfire'])
d_filt2 = stats('filt2      vs hapFIRE', m['af_filt2'], m['af_hapfire'])

# Outlier transitions
old_outl   = np.abs(d_old)   > 0.10
filt2_outl = np.abs(d_filt2) > 0.10
fixed    = old_outl & ~filt2_outl
introduced = ~old_outl & filt2_outl
print(f'\n=== Outlier transitions (|Δ|>0.10, mixedloose→filt2) ===')
print(f'  fixed by filt2:                {fixed.sum():,}  ({100*fixed.sum()/old_outl.sum():.1f}% of old outliers)')
print(f'  introduced by filt2:           {introduced.sum():,}')
print(f'  still outlier in both:         {(old_outl & filt2_outl).sum():,}')

# h-bias comparison
print(f'\n=== h-vector comparison ===')
print(f'  mixedloose CV: {h_old.std()/h_old.mean():.4f}')
print(f'  filt2      CV: {h_filt2.std()/h_filt2.mean():.4f}')
# cactus vs PG h-bias
import json
with open('/global/scratch/users/tbellg/kmate/data/founder_split_cactus_pg.json') as f:
    split = json.load(f)
is_cactus = np.array([str(s) in set(map(str, split['cactus'])) for s in cn_founders])
for label, h in [('mixedloose', h_old), ('filt2', h_filt2)]:
    F = len(h)
    uniform = 1.0 / F
    h_bias = (h[is_cactus].mean() / h[~is_cactus].mean())
    print(f'  {label:<12}  cactus_mean_h/PG_mean_h = {h_bias:.4f}  (uniform=1.000)  cactus_sum={h[is_cactus].sum():.4f} (unif={is_cactus.sum()/F:.4f})')

# Save 3-way comparison TSV
out = m[['chrom_num','pos','ref','alt','af_hapfire','af_old','af_filt2']].copy()
out['d_old']   = d_old.values
out['d_filt2'] = d_filt2.values
out.to_csv('SEEDMIX_S1_chr1_filt2_vs_mixedloose_vs_hapfire.tsv', sep='\t', index=False)
print(f'\nWrote SEEDMIX_S1_chr1_filt2_vs_mixedloose_vs_hapfire.tsv  ({len(out):,} rows)')
print(f'Total time: {time.time()-t0:.1f}s')
PYEOF

echo
echo "[$(date)] DONE F1"
