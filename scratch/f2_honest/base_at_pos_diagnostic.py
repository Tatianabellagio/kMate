#!/usr/bin/env python3
"""
Focused diagnostic on the base-at-pos fix ALONE (no F2-honest, no raw VCF).

Re-presents the existing base_at_pos_test.tsv results with a clean per-cell
view, and dissects the residual that base-at-pos cannot reach.
"""
import pandas as pd
import numpy as np
from pathlib import Path

ROOT = Path('/carnegie/nobackup/scratch/tbellagio/hapfire_sv')
M = pd.read_csv(ROOT / 'scratch/f2_honest/base_at_pos_test.tsv', sep='\t')
M['cell'] = (M['high_fmiss'].astype(int)*2 + M['is_mixed_bubble']).map({
    0:'[a] low_F pure_SNP',
    1:'[b] low_F mixed_bub',
    2:'[c] high_F pure_SNP',
    3:'[d] high_F mixed_bub',
})
print(f'rows: {len(M):,}')

# ---- Per-cell summary ----
print('\n=== Per-cell agreement (hapFIRE − estimator) ===')
rows = []
for cell, sub in M.groupby('cell'):
    r_orig = sub['hapfire_af'] - sub['alt_freq']
    r_base = sub['hapfire_af'] - sub['cactus_em_base']
    rows.append({
        'cell': cell, 'n': len(sub),
        'mean_orig':    r_orig.mean(),
        'MAE_orig':     r_orig.abs().mean(),
        'mean_baseatpos':r_base.mean(),
        'MAE_baseatpos':  r_base.abs().mean(),
        'mae_reduction_%':(1 - r_base.abs().mean()/r_orig.abs().mean())*100,
        'out05_orig':   (r_orig.abs()>0.05).mean()*100,
        'out05_base':   (r_base.abs()>0.05).mean()*100,
    })
print(pd.DataFrame(rows).round(4).to_string(index=False))

print('\n=== Genome-wide totals (516k joined Chr1 SNPs) ===')
for est, name in [('alt_freq',        'cactus_em (orig)'),
                  ('cactus_em_base',  'cactus_em + base-at-pos'),
                  ('recipe_af',       'recipe (orig)'),
                  ('recipe_base',     'recipe + base-at-pos')]:
    r = M['hapfire_af'] - M[est]
    print(f'  {name:32s}  mean={r.mean():+.5f}  MAE={r.abs().mean():.5f}  '
          f'R²={1 - (r**2).sum()/((M["hapfire_af"]-M["hapfire_af"].mean())**2).sum():.5f}  '
          f'%|d|>0.05={(r.abs()>0.05).mean()*100:.2f}%  %|d|>0.10={(r.abs()>0.10).mean()*100:.2f}%')

# ---- Dissect what's left on cell [d] after base-at-pos ----
print('\n=== What remains on cell [d] after base-at-pos? ===')
D = M[M['cell']=='[d] high_F mixed_bub'].copy()
D['resid_base'] = D['hapfire_af'] - D['cactus_em_base']
print(f'  cell [d] n: {len(D):,}')
print(f'  mean resid post-base-at-pos: {D["resid_base"].mean():+.4f}')
print(f'  MAE post-base-at-pos       : {D["resid_base"].abs().mean():.4f}')

# Stratify by extra_carriers (how much base-at-pos fired)
D['extra_bin'] = pd.cut(D['extra_carriers'],
                       bins=[-1, 0, 1, 10, 50, 230],
                       labels=['0','1','2-10','11-50','>50'])
print('\n  Residual stratified by extra_carriers added (base-at-pos magnitude):')
g = D.groupby('extra_bin', observed=True).apply(lambda s: pd.Series({
    'n': len(s),
    'mean_orig':   (s['hapfire_af']-s['alt_freq']).mean(),
    'mean_base':   (s['hapfire_af']-s['cactus_em_base']).mean(),
    'MAE_orig':    (s['hapfire_af']-s['alt_freq']).abs().mean(),
    'MAE_base':    (s['hapfire_af']-s['cactus_em_base']).abs().mean(),
})).round(4)
print(g.to_string())

# Records where base-at-pos didn't fire (extra=0): these are pure path mismatches
no_fire = D[D['extra_carriers']==0]
fired   = D[D['extra_carriers']>0]
print(f'\n  cell [d] records where base-at-pos DID NOT fire (extra=0): {len(no_fire):,}')
print(f'    mean resid: {(no_fire["hapfire_af"]-no_fire["alt_freq"]).mean():+.4f}'
      f'    MAE: {(no_fire["hapfire_af"]-no_fire["alt_freq"]).abs().mean():.4f}')

print(f'\n  cell [d] records where base-at-pos DID fire: {len(fired):,}')
fired_orig = (fired['hapfire_af']-fired['alt_freq'])
fired_base = (fired['hapfire_af']-fired['cactus_em_base'])
print(f'    mean resid orig: {fired_orig.mean():+.4f}    base-at-pos: {fired_base.mean():+.4f}')
print(f'    MAE orig: {fired_orig.abs().mean():.4f}    base-at-pos: {fired_base.abs().mean():.4f}')

# What % of cell [d] are "fully solved" (resid < 0.05 after base-at-pos)?
fully_solved = (D['resid_base'].abs() < 0.05).sum()
print(f'\n  Cell [d] records fully solved (|resid|<0.05 after base-at-pos): '
      f'{fully_solved:,} of {len(D):,} ({100*fully_solved/len(D):.1f}%)')
print(f'  ... was: {(((D["hapfire_af"]-D["alt_freq"]).abs())<0.05).sum():,} '
      f'({100*((D["hapfire_af"]-D["alt_freq"]).abs()<0.05).mean():.1f}%) before fix.')

# Look at the spot-check
print('\n=== Spot-check Chr1:13843898 ===')
sp = M[(M['chrom']=='Chr1') & (M['pos']==13843898)]
print(sp[['chrom','pos','ref','alt','AC_post','AN_post','extra_carriers','n_coalesced_recs',
          'alt_freq','cactus_em_base','recipe_base','hapfire_af']].to_string(index=False))

# Find a few records where base-at-pos still leaves a big residual
print('\n=== 10 cell [d] records where base-at-pos still leaves the biggest residual ===')
worst = D.nlargest(10, 'resid_base')
print(worst[['pos','ref','alt','AC_post','AN_post','extra_carriers','n_coalesced_recs',
             'alt_freq','cactus_em_base','hapfire_af','resid_base']].round(4).to_string(index=False))
