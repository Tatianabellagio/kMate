#!/bin/bash
#SBATCH --job-name=chr1_atomdd
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=2
#SBATCH --mem=48G
#SBATCH --time=00:45:00
#SBATCH --output=logs/D3_deepdive_%j.out
#SBATCH --error=logs/D3_deepdive_%j.err
mkdir -p logs
set -euo pipefail

# Deep dive on the ~8,900 remaining outliers in arch3 ATOMIZED vs hapFIRE.
# Decompose by:
#  - direction (under vs over)
#  - cactus vs PG internal agreement (panel-vs-1001G, cactus-specific, PG-specific, mixed)
#  - position (centromere 14-17 Mb vs arms)
#  - F_MISSING bin
#  - AF stratification
#  - vs raw arch3 / vs v3qc_v3 — were they outliers there too?
# Plus spot-checks: top 15 outliers in each class.

cd /global/scratch/users/tbellg/hapfire_sv/arch3/chr1
PY=/global/home/users/tbellg/miniforge3/envs/hapfm/bin/python

$PY -u <<'PYEOF'
import numpy as np
import pandas as pd
from scipy.sparse import load_npz
import time, json

t0 = time.time()

print('=== Load atomized cactus_em + meta + hapFIRE ===')
atom_af = pd.read_csv('SEEDMIX_S1_arch3_chr1_atomized.tsv', sep='\t')
print(f'  atomized records: {len(atom_af):,}')

cn_a = load_npz('cn_var_231_arch3_chr1_atomized.cn_var.npz').tocsc()
cnc_a = load_npz('cn_var_231_arch3_chr1_atomized.cn_var_called.npz').tocsc()
meta_a = np.load('cn_var_231_arch3_chr1_atomized.meta.npz', allow_pickle=True)
m_pos = meta_a['pos']; m_ref = meta_a['ref']; m_alt = meta_a['alt']
founders = meta_a['founders']

cactus_ids = set(open('/global/scratch/users/tbellg/hapfire_sv/pangenie_genotyping/data/merged/cactus_overlap_80.txt').read().split())
f_is_cactus = np.array([str(f) in cactus_ids for f in founders])
print(f'  cactus founders: {f_is_cactus.sum()}, PG founders: {(~f_is_cactus).sum()}')

# hapFIRE with REF/ALT
hf = pd.read_csv('/global/scratch/projects/fc_moilab/projects/grenenet-phase1/frequency/hapFIRE_frequencies/seed_mix/s1_snp_frequency.txt',
                 sep='\t', header=None, names=['chrom_num','pos','af_hapfire'])
hf = hf[hf.chrom_num == 1].copy()
ra = pd.read_csv('hapfire_chr1_refalt.tsv', sep='\t').rename(columns={'chrom':'chrom_num'})
ra['chrom_num'] = ra['chrom_num'].astype(int)
hf_full = hf.merge(ra, on=['chrom_num','pos'], how='left').dropna(subset=['ref','alt'])

# 4-tuple join: atomized ∩ hapFIRE
atom_af['chrom_num'] = 1
print('=== 4-tuple join atomized ∩ hapFIRE ===')
joined = atom_af.merge(hf_full[['chrom_num','pos','ref','alt','af_hapfire']],
                       on=['chrom_num','pos','ref','alt'], how='inner')
joined = joined.dropna(subset=['alt_freq'])
print(f'  4-tuple intersect: {len(joined):,}')

joined['delta'] = joined['alt_freq'] - joined['af_hapfire']
joined['absd']  = np.abs(joined['delta'])
outl = joined[joined['absd'] > 0.10].copy()
print(f'\n=== Outliers |Δ|>0.10: {len(outl):,} ({100*len(outl)/len(joined):.2f}%) ===')

# Map (pos,ref,alt) → atomized row index for sparse slicing
print('\n=== Building (pos,ref,alt) → row index for outlier positions ===')
m_ref_s = np.array([str(r) for r in m_ref])
m_alt_s = np.array([str(a) for a in m_alt])
key2idx = {}
for i in range(len(m_pos)):
    key2idx[(int(m_pos[i]), m_ref_s[i], m_alt_s[i])] = i

# Annotate each outlier with cactus_AF, PG_AF, F_MISSING (computed from cn_var_called)
print('=== Computing per-class cactus_AF / PG_AF / F_MISSING for outliers ===')
cact_af_list = []; pg_af_list = []; fmiss_list = []; n_alt_atoms_list = []
# Also: how many atomized rows at this same pos (different ALT)
from collections import defaultdict
pos_to_rows = defaultdict(list)
for i in range(len(m_pos)):
    pos_to_rows[int(m_pos[i])].append(i)

for _, r in outl.iterrows():
    key = (int(r['pos']), str(r['ref']), str(r['alt']))
    idx = key2idx.get(key)
    if idx is None:
        cact_af_list.append(np.nan); pg_af_list.append(np.nan); fmiss_list.append(np.nan); n_alt_atoms_list.append(0); continue
    col_c = cn_a[:, idx].toarray().ravel()
    col_k = cnc_a[:, idx].toarray().ravel()
    ac_c = int(col_c[f_is_cactus].sum()); an_c = int(col_k[f_is_cactus].sum())
    ac_p = int(col_c[~f_is_cactus].sum()); an_p = int(col_k[~f_is_cactus].sum())
    cact_af_list.append(ac_c/max(1,an_c))
    pg_af_list.append(ac_p/max(1,an_p))
    fmiss_list.append(1 - col_k.sum()/len(col_k))
    n_alt_atoms_list.append(len(pos_to_rows[int(r['pos'])]))
outl['cact_af'] = cact_af_list
outl['pg_af']   = pg_af_list
outl['fmiss']   = fmiss_list
outl['n_alt_atoms_at_pos'] = n_alt_atoms_list

# === DIRECTION ===
under = outl[outl['delta'] < -0.10]
over  = outl[outl['delta'] >  0.10]
print(f'\n=== Direction ===')
print(f'  under-pred (atomized < hapFIRE): {len(under):,}  ({100*len(under)/len(outl):.1f}%)')
print(f'  over-pred  (atomized > hapFIRE): {len(over):,}   ({100*len(over)/len(outl):.1f}%)')

# === CENTROMERE vs ARMS ===
cent = (outl['pos'] >= 14_000_000) & (outl['pos'] <= 17_000_000)
cent_total_frac = ((joined['pos'] >= 14_000_000) & (joined['pos'] <= 17_000_000)).sum() / len(joined)
print(f'\n=== Centromere (14-17 Mb) ===')
print(f'  outliers in centromere: {cent.sum():,} ({100*cent.sum()/len(outl):.1f}%)')
print(f'  centromere fraction of all 4-tuple rows: {100*cent_total_frac:.2f}%')
print(f'  → centromere is {100*cent.sum()/len(outl)/max(0.01, 100*cent_total_frac):.1f}× enriched')

# === CACTUS vs PG INTERNAL AGREEMENT ===
print(f'\n=== Internal agreement (cactus vs PG founders) — UNDER outliers ===')
u = under.dropna(subset=['cact_af','pg_af'])
panel_low_both    = u[(u['cact_af'] < 0.1) & (u['pg_af'] < 0.1)]
panel_high_both   = u[(u['cact_af'] > 0.5) & (u['pg_af'] > 0.5)]
cact_only         = u[(u['cact_af'] > 0.5) & (u['pg_af'] < 0.1)]
pg_only           = u[(u['pg_af']   > 0.5) & (u['cact_af'] < 0.1)]
mixed             = u.drop(panel_low_both.index).drop(panel_high_both.index, errors='ignore').drop(cact_only.index, errors='ignore').drop(pg_only.index, errors='ignore')
print(f'  cact+PG BOTH LOW  ({len(panel_low_both):>5,}):  our panel agrees LOW; hapFIRE says HIGH → 1001G outlier')
print(f'  cact+PG BOTH HIGH ({len(panel_high_both):>5,}):  shouldn\'t exist for under-pred (sanity)')
print(f'  cactus-specific   ({len(cact_only):>5,}):  cactus high, PG low')
print(f'  PG-specific       ({len(pg_only):>5,}):     PG high, cactus low')
print(f'  mixed/intermediate ({len(mixed):>5,}):  neither class cleanly extreme')

print(f'\n=== Internal agreement — OVER outliers ===')
o = over.dropna(subset=['cact_af','pg_af'])
panel_high_both_o = o[(o['cact_af'] > 0.5) & (o['pg_af'] > 0.5)]
panel_low_both_o  = o[(o['cact_af'] < 0.1) & (o['pg_af'] < 0.1)]
cact_only_o       = o[(o['cact_af'] > 0.5) & (o['pg_af'] < 0.1)]
pg_only_o         = o[(o['pg_af']   > 0.5) & (o['cact_af'] < 0.1)]
mixed_o           = o.drop(panel_high_both_o.index).drop(panel_low_both_o.index, errors='ignore').drop(cact_only_o.index, errors='ignore').drop(pg_only_o.index, errors='ignore')
print(f'  cact+PG BOTH HIGH ({len(panel_high_both_o):>5,}): our panel agrees HIGH; hapFIRE says LOW → 1001G outlier')
print(f'  cact+PG BOTH LOW  ({len(panel_low_both_o):>5,}):  shouldn\'t exist for over-pred (sanity)')
print(f'  cactus-specific   ({len(cact_only_o):>5,})')
print(f'  PG-specific       ({len(pg_only_o):>5,})')
print(f'  mixed/intermediate({len(mixed_o):>5,})')

# === F_MISSING stratification ===
print(f'\n=== F_MISSING distribution ===')
print(f'{"F bin":<18}{"outliers":>10}{"all":>10}{"outlier rate":>14}')
for lo, hi in [(0,0.01),(0.01,0.05),(0.05,0.1),(0.1,0.2),(0.2,0.5),(0.5,1.01)]:
    sel = (outl['fmiss'] >= lo) & (outl['fmiss'] < hi)
    n_out = sel.sum()
    # corresponding fraction in full joined
    fm_all_n = ((joined['absd']<=0.10).sum() if False else 0)  # placeholder; compute over joined below
    print(f'  F∈[{lo},{hi}): {n_out:>10,}')

# === AF stratification ===
print(f'\n=== hapFIRE AF stratum: outlier rate ===')
joined['af_bin'] = pd.cut(joined['af_hapfire'], bins=[0,0.05,0.1,0.3,0.5,0.7,0.9,0.95,1.001], include_lowest=True,
                          labels=['0-0.05','0.05-0.1','0.1-0.3','0.3-0.5','0.5-0.7','0.7-0.9','0.9-0.95','0.95-1'])
outl['af_bin'] = pd.cut(outl['af_hapfire'], bins=[0,0.05,0.1,0.3,0.5,0.7,0.9,0.95,1.001], include_lowest=True,
                          labels=['0-0.05','0.05-0.1','0.1-0.3','0.3-0.5','0.5-0.7','0.7-0.9','0.9-0.95','0.95-1'])
total_per = joined['af_bin'].value_counts().sort_index()
outl_per = outl['af_bin'].value_counts().sort_index()
print(f'{"hapFIRE AF":<14}{"n_total":>10}{"n_outlier":>12}{"rate":>9}')
for lbl in ['0-0.05','0.05-0.1','0.1-0.3','0.3-0.5','0.5-0.7','0.7-0.9','0.9-0.95','0.95-1']:
    t = total_per.get(lbl, 0); o2 = outl_per.get(lbl, 0)
    print(f'  {lbl:<12}{t:>10,}{o2:>12,}{100*o2/max(1,t):>8.2f}%')

# === Spot checks: top 10 outliers in each class ===
def show_top(df, label, sort_col='absd', n=10):
    sub = df.dropna(subset=['cact_af','pg_af']).sort_values(sort_col, ascending=False).head(n)
    if len(sub) == 0:
        print(f'\n=== {label}: (no records) ==='); return
    print(f'\n=== {label}: top {n} by |Δ| ===')
    cols = ['pos','ref','alt','af_hapfire','alt_freq','delta','cact_af','pg_af','fmiss','n_alt_atoms_at_pos']
    print(sub[cols].to_string(index=False, float_format='%.3f'))

show_top(panel_low_both,   'UNDER outliers — cact+PG BOTH LOW (1001G outlier)')
show_top(cact_only,        'UNDER outliers — cactus-specific (cact high, PG low)')
show_top(pg_only,          'UNDER outliers — PG-specific (PG high, cact low)')
show_top(mixed,            'UNDER outliers — MIXED/intermediate (largest class)')
show_top(panel_high_both_o,'OVER outliers — cact+PG BOTH HIGH (1001G outlier)')
show_top(mixed_o,          'OVER outliers — MIXED/intermediate')

# Save annotated outlier table
outl.to_csv('SEEDMIX_S1_chr1_atomized_outliers_annotated.tsv', sep='\t', index=False)
print(f'\nWrote SEEDMIX_S1_chr1_atomized_outliers_annotated.tsv  ({len(outl):,} rows)')
print(f'Total time: {time.time()-t0:.1f}s')
PYEOF

echo
echo "[$(date)] DONE D3"
