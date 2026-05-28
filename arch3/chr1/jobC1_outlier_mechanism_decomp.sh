#!/bin/bash
#SBATCH --job-name=outl_decomp
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH --time=00:45:00
#SBATCH --output=logs/C1_decomp_%j.out
#SBATCH --error=logs/C1_decomp_%j.err
mkdir -p logs
set -euo pipefail

# Decompose the 3,344 outliers (|Δ|>0.10) from A7 into mechanism classes:
#   1. Encoding disagreement (multi-ALT at pos; UNION AF matches hapFIRE)
#   2. Genuine GT disagreement at single-ALT positions
#   3. Genuine GT disagreement at multi-ALT positions
# Then for each class, split by direction (under vs over) and check whether
# cactus and PG founders internally agree (panel-vs-1001G) or disagree (PG-specific).

cd /global/scratch/users/tbellg/hapfire_sv/arch3/chr1
PY=/global/home/users/tbellg/miniforge3/envs/hapfm/bin/python

$PY -u <<'PYEOF'
import pandas as pd
import numpy as np
from scipy.sparse import load_npz
from collections import defaultdict
import time

t0 = time.time()
print('=== Load A7 comparison TSV ===')
df = pd.read_csv('SEEDMIX_S1_chr1_AF_compare.tsv', sep='\t')
df['delta'] = df['af_new'] - df['af_hapfire']
df['absd']  = np.abs(df['delta'])
outl = df[df['absd'] > 0.10].copy().reset_index(drop=True)
print(f'  total rows: {len(df):,}, outliers (|Δ|>0.10): {len(outl):,}')

print('=== Load Arch 3 cn_var + meta ===')
cn = load_npz('cn_var_231_arch3_chr1.cn_var.npz').tocsc()
cn_called = load_npz('cn_var_231_arch3_chr1.cn_var_called.npz').tocsc()
meta = np.load('cn_var_231_arch3_chr1.meta.npz', allow_pickle=True)
m_pos = meta['pos']
m_ref = meta['ref']
m_alt = meta['alt']
founders = meta['founders']
print(f'  cn_var: {cn.shape}, {cn.nnz:,} nnz')

cactus_ids = set(open('/global/scratch/users/tbellg/hapfire_sv/pangenie_genotyping/data/merged/cactus_overlap_80.txt').read().split())
f_is_cactus = np.array([str(f) in cactus_ids for f in founders])
print(f'  cactus founders: {f_is_cactus.sum()},  PG founders: {(~f_is_cactus).sum()}')

print('=== Indexing positions (full chr1) ===')
pos_to_indices = defaultdict(list)
for i in range(len(m_pos)):
    pos_to_indices[int(m_pos[i])].append(i)
print(f'  {len(pos_to_indices):,} unique positions, {time.time()-t0:.1f}s')

print('=== Computing UNION AF per outlier position ===')
unique_pos = outl['pos'].unique()
union_af_for_pos = {}
for j, pos_i in enumerate(unique_pos):
    rows = pos_to_indices.get(int(pos_i), [])
    if not rows:
        continue
    union_carrier = np.zeros(cn.shape[0], dtype=bool)
    union_called  = np.zeros(cn.shape[0], dtype=bool)
    for ri in rows:
        c = cn[:, ri].toarray().ravel()
        cc = cn_called[:, ri].toarray().ravel()
        union_carrier |= (c > 0)
        union_called  |= (cc > 0)
    union_af_for_pos[int(pos_i)] = (int(union_carrier.sum()), int(union_called.sum()), len(rows),
                                      int(union_carrier[f_is_cactus].sum()), int(union_called[f_is_cactus].sum()),
                                      int(union_carrier[~f_is_cactus].sum()), int(union_called[~f_is_cactus].sum()))
    if (j+1) % 500 == 0:
        print(f'  {j+1}/{len(unique_pos)} positions, {time.time()-t0:.1f}s')

outl['n_alt_rows'] = outl['pos'].map(lambda p: union_af_for_pos.get(int(p), (0,0,0,0,0,0,0))[2])
outl['union_ac']   = outl['pos'].map(lambda p: union_af_for_pos.get(int(p), (0,0,0,0,0,0,0))[0])
outl['union_an']   = outl['pos'].map(lambda p: union_af_for_pos.get(int(p), (0,0,0,0,0,0,0))[1])
outl['union_af']   = outl['union_ac'] / outl['union_an'].clip(lower=1)
outl['union_absd'] = np.abs(outl['union_af'] - outl['af_hapfire'])
outl['union_cact_ac'] = outl['pos'].map(lambda p: union_af_for_pos.get(int(p), (0,0,0,0,0,0,0))[3])
outl['union_cact_an'] = outl['pos'].map(lambda p: union_af_for_pos.get(int(p), (0,0,0,0,0,0,0))[4])
outl['union_pg_ac']   = outl['pos'].map(lambda p: union_af_for_pos.get(int(p), (0,0,0,0,0,0,0))[5])
outl['union_pg_an']   = outl['pos'].map(lambda p: union_af_for_pos.get(int(p), (0,0,0,0,0,0,0))[6])
outl['cact_af'] = outl['union_cact_ac'] / outl['union_cact_an'].clip(lower=1)
outl['pg_af']   = outl['union_pg_ac']   / outl['union_pg_an'].clip(lower=1)

print()
print('=' * 70)
print(f'=== MECHANISM DECOMPOSITION of {len(outl):,} outliers (|Δ|>0.10) ===')
print('=' * 70)

multi_alt = outl[outl['n_alt_rows'] >= 2]
single_alt = outl[outl['n_alt_rows'] == 1]
encoding = outl[(outl['n_alt_rows'] >= 2) & (outl['union_absd'] <= 0.05)]
genuine_multi = outl[(outl['n_alt_rows'] >= 2) & (outl['union_absd'] > 0.05)]
genuine_single = outl[outl['n_alt_rows'] == 1]
genuine = pd.concat([genuine_multi, genuine_single])

print(f'\nAt the position level (across all ALTs):')
print(f'  multi-ALT outliers (n_alt_rows ≥ 2): {len(multi_alt):>6,}  ({100*len(multi_alt)/len(outl):.1f}%)')
print(f'    of which ENCODING-disagreement (union_AF matches hapFIRE within 0.05):')
print(f'      {len(encoding):>6,}  ({100*len(encoding)/len(outl):.1f}%)  ← NOT REAL DISAGREEMENT')
print(f'    of which GENUINE GT disagreement at multi-ALT position:')
print(f'      {len(genuine_multi):>6,}  ({100*len(genuine_multi)/len(outl):.1f}%)')
print(f'  single-ALT outliers (n_alt_rows = 1):  {len(single_alt):>6,}  ({100*len(single_alt)/len(outl):.1f}%)  ← pure GT disagreement')

print()
print('=== Bottom-line classification ===')
total = len(outl)
print(f'  ENCODING-disagreement (cactus-MNP-vs-1001G-SNP splits): {len(encoding):>6,} ({100*len(encoding)/total:.1f}%)')
print(f'  GENUINE per-founder GT disagreement                    : {len(genuine):>6,} ({100*len(genuine)/total:.1f}%)')

print()
print('=== GENUINE outliers — internal agreement between cactus & PG ===')
under = genuine[genuine['delta'] < -0.10]
over  = genuine[genuine['delta'] >  0.10]

# Under-prediction: our panel says LOW, hapFIRE says HIGH
print(f'\n  UNDER-prediction (NEW < hapFIRE): {len(under):,}')
both_low = under[(under['cact_af'] < 0.1) & (under['pg_af'] < 0.1)]
cact_high_pg_low = under[(under['cact_af'] > 0.5) & (under['pg_af'] < 0.1)]
pg_high_cact_low = under[(under['pg_af'] > 0.5) & (under['cact_af'] < 0.1)]
print(f'    cactus + PG both LOW (both AF < 0.1):  {len(both_low):>5,}  ({100*len(both_low)/max(1,len(under)):.0f}%) ← our panel agrees internally; 1001G is the outlier')
print(f'    cactus HIGH, PG LOW (cactus-specific): {len(cact_high_pg_low):>5,}  ({100*len(cact_high_pg_low)/max(1,len(under)):.0f}%)')
print(f'    PG HIGH, cactus LOW (PG-specific):     {len(pg_high_cact_low):>5,}  ({100*len(pg_high_cact_low)/max(1,len(under)):.0f}%)')

# Over-prediction: our panel says HIGH, hapFIRE says LOW
print(f'\n  OVER-prediction (NEW > hapFIRE): {len(over):,}')
both_high = over[(over['cact_af'] > 0.5) & (over['pg_af'] > 0.5)]
cact_high_pg_low2 = over[(over['cact_af'] > 0.5) & (over['pg_af'] < 0.1)]
pg_high_cact_low2 = over[(over['pg_af'] > 0.5) & (over['cact_af'] < 0.1)]
print(f'    cactus + PG both HIGH (both AF > 0.5): {len(both_high):>5,}  ({100*len(both_high)/max(1,len(over)):.0f}%) ← our panel agrees internally; 1001G is the outlier')
print(f'    cactus HIGH, PG LOW (cactus-specific): {len(cact_high_pg_low2):>5,}  ({100*len(cact_high_pg_low2)/max(1,len(over)):.0f}%)')
print(f'    PG HIGH, cactus LOW (PG-specific):     {len(pg_high_cact_low2):>5,}  ({100*len(pg_high_cact_low2)/max(1,len(over)):.0f}%)')

# Write the annotated outliers TSV
out_tsv = 'SEEDMIX_S1_chr1_outliers_annotated.tsv'
outl.to_csv(out_tsv, sep='\t', index=False)
print(f'\nWrote {out_tsv} ({len(outl):,} rows)')
print(f'\nTotal time: {time.time()-t0:.1f}s')
PYEOF

echo
echo "[$(date)] DONE C1"
