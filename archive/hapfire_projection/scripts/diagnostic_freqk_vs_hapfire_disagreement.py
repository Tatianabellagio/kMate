"""
Diagnostic: where does freqk disagree with hapFIRE-proj on the SEEDMIX data,
and is that disagreement structural (concentrated in singletons / repeat-prone
SVs) or widespread (also present at mid-AC, non-repeat SVs)?

The answer steers the LD-borrowing design:
  - if disagreement concentrates in repeat / singleton SVs only -> a posterior
    combine + a repeat filter is enough (cheap)
  - if it's widespread -> joint-CVXPY on SNP+SV evidence is worth the work
"""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

ROOT = Path('/carnegie/nobackup/scratch/tbellagio/hapfire_sv')
OUT = ROOT / 'results' / 'diagnostic_disagreement'
OUT.mkdir(exist_ok=True, parents=True)

df = pd.read_csv(ROOT / 'results' / 'seedmix_3way_comparison.tsv.gz', sep='\t')
df = df[df['var_type'].isin(['INS', 'DEL'])].copy()

# proxy 1 for "freqk struggled here": per-SV NaN/zero rate across the 8 samples
per_sv = (df.groupby(['chrom','pos','alt_idx','var_type','sv_size','ac_in_panel'])
            .agg(freqk_nan_n=('freqk_af', lambda x: x.isna().sum()),
                 freqk_zero_n=('freqk_af', lambda x: (x.fillna(-1) == 0).sum()),
                 freqk_med=('freqk_af', 'median'),
                 hap_med=('hapfire_proj_af', 'median'),
                 truth_med=('truth_panel', 'median'))
            .reset_index())
per_sv['freqk_avail_n'] = 8 - per_sv['freqk_nan_n']
per_sv['freqk_repeat_drop'] = (per_sv['freqk_avail_n'] == 0)
per_sv['signed_diff'] = per_sv['freqk_med'] - per_sv['hap_med']
per_sv['abs_diff']    = per_sv['signed_diff'].abs()

# AC binning
def ac_bin(ac):
    if ac == 1:    return '1 (singleton)'
    if ac <= 5:    return '2-5'
    if ac <= 20:   return '6-20'
    if ac <= 50:   return '21-50'
    return '50+'
per_sv['ac_bin'] = per_sv['ac_in_panel'].apply(ac_bin)
order = ['1 (singleton)','2-5','6-20','21-50','50+']

# ---------- TABLE 1: AC-stratified disagreement ----------
print('='*100)
print('TABLE 1 — disagreement summary by AC bin (n_SV, fraction freqk-dropped, |diff| stats)')
print('='*100)
rows = []
for b in order:
    sub = per_sv[per_sv['ac_bin']==b]
    if len(sub)==0: continue
    avail = sub[~sub['freqk_repeat_drop']]
    rows.append({
        'ac_bin': b,
        'n_SV': len(sub),
        'pct_freqk_dropped_8of8': f'{sub["freqk_repeat_drop"].mean():.1%}',
        'pct_with_|diff|>0.05': f'{(avail["abs_diff"]>0.05).mean():.1%}',
        'pct_with_|diff|>0.10': f'{(avail["abs_diff"]>0.10).mean():.1%}',
        'pct_with_|diff|>0.20': f'{(avail["abs_diff"]>0.20).mean():.1%}',
        'median_|diff|':         f'{avail["abs_diff"].median():.4f}',
        'mean_|diff|':           f'{avail["abs_diff"].mean():.4f}',
        'pct_freqk>hap_by_0.10': f'{(avail["signed_diff"]>0.10).mean():.1%}',
        'pct_hap>freqk_by_0.10': f'{(-avail["signed_diff"]>0.10).mean():.1%}',
    })
t1 = pd.DataFrame(rows)
print(t1.to_string(index=False))
t1.to_csv(OUT / 'table1_disagreement_by_AC.tsv', sep='\t', index=False)

# ---------- TABLE 2: agreement metric (per-SV correlation across 8 samples) ----------
# is freqk's residual variance explainable by hapFIRE in non-repeat / non-singleton SVs?
def per_sv_corr(g):
    a, b = g['freqk_af'].values, g['hapfire_proj_af'].values
    mask = ~np.isnan(a) & ~np.isnan(b)
    if mask.sum() < 3 or np.std(a[mask])==0 or np.std(b[mask])==0:
        return np.nan
    return np.corrcoef(a[mask], b[mask])[0,1]
ps = (df.groupby(['chrom','pos','alt_idx','ac_in_panel'])
        .apply(per_sv_corr, include_groups=False).rename('r').reset_index())
ps['ac_bin'] = ps['ac_in_panel'].apply(ac_bin)
print()
print('='*100)
print('TABLE 2 — across-sample Pearson r per SV (freqk vs hapFIRE-proj), by AC')
print('   (a low r means: even *direction* of variation across the 8 samples disagrees)')
print('='*100)
t2 = ps.groupby('ac_bin')['r'].agg(['count','mean','median', lambda x: (x<0).mean()])
t2.columns = ['n_with_corr','mean_r','median_r','pct_negative_r']
t2 = t2.reindex(order)
t2.to_csv(OUT / 'table2_per_sv_correlation_by_AC.tsv', sep='\t')
print(t2.to_string())

# ---------- TABLE 3: same disagreement, but excluding SVs where panel-truth says low AF ----------
# (panel-truth = recipe-based, only counts the 80 panel founders. Singletons in panel
# have truth ~ 0.005 and disagreement is near-trivially explained by panel ascertainment.)
print()
print('='*100)
print('TABLE 3 — disagreement on SVs where panel-truth median > 0.10 (i.e. SVs the panel says are common)')
print('   — this is the "real-info" subset; panel-ascertainment cant explain disagreement here')
print('='*100)
common = per_sv[per_sv['truth_med']>0.10].copy()
print(f'  n_SV with truth_med > 0.10: {len(common)}')
print(f'  median |diff|: {common["abs_diff"].median():.4f}')
print(f'  pct |diff|>0.05: {(common["abs_diff"]>0.05).mean():.1%}')
print(f'  pct |diff|>0.10: {(common["abs_diff"]>0.10).mean():.1%}')
print(f'  pct |diff|>0.20: {(common["abs_diff"]>0.20).mean():.1%}')
print(f'  pct freqk>hap by >0.10: {(common["signed_diff"]>0.10).mean():.1%}')
print(f'  pct hap>freqk by >0.10: {(-common["signed_diff"]>0.10).mean():.1%}')

# ---------- TABLE 4: bin by SV size  ----------
def size_bin(s):
    if s <= 100:   return '<=100bp'
    if s <= 500:   return '101-500bp'
    if s <= 2000:  return '501-2000bp'
    if s <= 10000: return '2-10kb'
    return '>10kb'
per_sv['size_bin'] = per_sv['sv_size'].apply(size_bin)
print()
print('='*100)
print('TABLE 4 — disagreement by SV size  (small SVs have fewer unique flanking k-mers)')
print('='*100)
sb_order = ['<=100bp','101-500bp','501-2000bp','2-10kb','>10kb']
rows=[]
for b in sb_order:
    sub = per_sv[per_sv['size_bin']==b]
    if len(sub)==0: continue
    avail = sub[~sub['freqk_repeat_drop']]
    rows.append({
        'size_bin': b, 'n_SV': len(sub),
        'pct_dropped': f'{sub["freqk_repeat_drop"].mean():.1%}',
        'median_|diff|': f'{avail["abs_diff"].median():.4f}',
        'pct_|diff|>0.10': f'{(avail["abs_diff"]>0.10).mean():.1%}',
    })
t4 = pd.DataFrame(rows)
print(t4.to_string(index=False))
t4.to_csv(OUT / 'table4_disagreement_by_size.tsv', sep='\t', index=False)

# ---------- FIGURE: scatter of freqk vs hapFIRE colored by AC bin ----------
fig, axes = plt.subplots(1, 5, figsize=(20,4), sharey=True)
for ax, b in zip(axes, order):
    sub = per_sv[(per_sv['ac_bin']==b) & (~per_sv['freqk_repeat_drop'])]
    if len(sub)==0:
        ax.set_title(f'{b} (n=0)'); continue
    ax.scatter(sub['hap_med'], sub['freqk_med'], s=2, alpha=0.2, c='steelblue', rasterized=True)
    ax.plot([0,1],[0,1],'k--',lw=0.5)
    ax.set_title(f'{b} (n={len(sub)})')
    ax.set_xlabel('hapFIRE-proj AF (median over 8 samples)')
    ax.set_xlim(-0.02,1.02); ax.set_ylim(-0.02,1.02)
axes[0].set_ylabel('freqk AF (median over 8 samples)')
plt.suptitle('SEEDMIX per-SV: freqk vs hapFIRE-proj, stratified by panel AC')
plt.tight_layout()
plt.savefig(OUT / 'fig_freqk_vs_hapfire_by_AC.png', dpi=140, bbox_inches='tight')
print(f'\nWrote {OUT}/fig_freqk_vs_hapfire_by_AC.png')

# ---------- BOTTOM-LINE summary line ----------
total = len(per_sv)
dropped = per_sv['freqk_repeat_drop'].sum()
avail = per_sv[~per_sv['freqk_repeat_drop']]
sing_drop = per_sv[(per_sv['ac_bin']=='1 (singleton)') & per_sv['freqk_repeat_drop']]
mid = avail[avail['ac_bin'].isin(['6-20','21-50','50+'])]
mid_bad = mid[mid['abs_diff']>0.10]
print()
print('='*100)
print('BOTTOM LINE')
print('='*100)
print(f'  total SVs: {total:,}')
print(f'  freqk dropped in 8/8 samples: {dropped:,} ({dropped/total:.1%}) [repeat-region failure]')
print(f'  of those, singletons (AC=1): {len(sing_drop):,}')
print(f'  freqk available, AC>=6:      {len(mid):,}')
print(f'    of those with |freqk - hap| > 0.10: {len(mid_bad):,} ({len(mid_bad)/max(1,len(mid)):.1%})')
print(f'    => disagreement is {"WIDESPREAD" if len(mid_bad)/max(1,len(mid))>0.05 else "concentrated"} in non-repeat / mid-AC SVs')

per_sv.to_csv(OUT / 'per_sv_disagreement.tsv', sep='\t', index=False)
print(f'\nWrote per-SV table: {OUT}/per_sv_disagreement.tsv')
