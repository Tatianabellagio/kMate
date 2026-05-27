#!/usr/bin/env python3
"""
Cell [d] validation v2 — using awk-extracted xwu GT TSV.

For each probe:
  - recipe_xwu = AC/AN computed from xwu's per-founder GTs
  - Compare to recipe_cactus, recipe_F2honest, hapFIRE
  - Per-founder comparison: at each probe, how many founders' calls disagree
    between cactus (cn_var) and xwu?
"""
import pandas as pd
import numpy as np
from pathlib import Path

ROOT = Path('/carnegie/nobackup/scratch/tbellagio/hapfire_sv')

# Load probe table (with cactus/F2-honest/hapFIRE values)
R = pd.read_csv(ROOT / 'scratch/f2_honest/celld_validation.tsv', sep='\t')
print(f'probes: {len(R)}')

# Load xwu GTs
xwu = pd.read_csv(ROOT / 'scratch/f2_honest/xwu_probe_gts.tsv', sep='\t')
sample_cols = list(xwu.columns[3:])
sample_cols = [s for s in sample_cols if s != '']
print(f'xwu samples: {len(sample_cols)}')

def gt_to_carrier(g):
    """Phased diploid GT -> (carrier_count, called_flag)."""
    g = str(g).strip()
    if g in ('.', './.', '.|.', '.|', '|.'):
        return (0, 0, 0)
    # phased "a|b" or unphased "a/b" or haploid "a"
    sep = '|' if '|' in g else ('/' if '/' in g else None)
    if sep is None:
        # haploid
        try:
            v = int(g)
            return (1 if v>0 else 0, 1 if v>0 else 0, 1)   # carrier(any), carrier_strict, called
        except:
            return (0, 0, 0)
    a, b = g.split(sep)
    if a == '.' or b == '.':
        return (0, 0, 0)
    a, b = int(a), int(b)
    # carrier any-ALT: 1 if a>0 or b>0
    carrier_any = 1 if (a>0 or b>0) else 0
    # carrier strict: both alleles ALT (hom-ALT)
    carrier_strict = 1 if (a>0 and b>0) else 0
    return (carrier_any, carrier_strict, 1)

# Compute recipe_xwu per probe
xwu['ref_xwu'] = xwu['ref']
xwu['alt_xwu'] = xwu['alt']

results = []
for _, p in R.iterrows():
    pos = int(p['pos'])
    xwu_row = xwu[xwu['pos']==pos]
    if len(xwu_row)==0:
        continue
    xwu_row = xwu_row.iloc[0]
    ref_xwu = str(xwu_row['ref_xwu']); alt_xwu = str(xwu_row['alt_xwu'])
    # Check REF/ALT alignment with cn_var
    if (p['ref'], p['alt']) == (ref_xwu, alt_xwu):
        swap = False
    elif (p['alt'], p['ref']) == (ref_xwu, alt_xwu):
        swap = True
    else:
        swap = None  # genuinely different alleles
    AC_any = AC_strict = AN = 0
    for s in sample_cols:
        c_any, c_strict, called = gt_to_carrier(xwu_row[s])
        AC_any += c_any
        AC_strict += c_strict
        AN += called
    if swap is None:
        recipe_xwu_any = recipe_xwu_strict = float('nan')
    else:
        af_any = AC_any/AN if AN>0 else float('nan')
        af_strict = AC_strict/AN if AN>0 else float('nan')
        if swap:
            af_any = 1 - af_any
            af_strict = 1 - af_strict
        recipe_xwu_any = af_any
        recipe_xwu_strict = af_strict
    results.append({
        'pos': pos, 'ref': p['ref'], 'alt': p['alt'],
        'ref_xwu': ref_xwu, 'alt_xwu': alt_xwu, 'alleles_swap': swap,
        'cactus_em':        p['cactus_em'],
        'cactus_em_F2hon':  p['cactus_em_F2hon'],
        'recipe_cactus':    p['recipe_cactus'],
        'recipe_F2hon':     p['recipe_F2hon'],
        'recipe_xwu_any':   recipe_xwu_any,
        'recipe_xwu_strict':recipe_xwu_strict,
        'hapfire_af':       p['hapfire_af'],
        'AC_cactus':        p['AC_cactus'],
        'AN_cactus':        p['AN_cactus'],
        'AC_xwu_any':       AC_any,
        'AC_xwu_strict':    AC_strict,
        'AN_xwu':           AN,
        'n_override':       p['n_override'],
        'n_sv_at_p':        p['n_sv_at_p'],
        'resid_orig':       p['resid_orig_vs_hf'],
        'resid_F2h':        p['resid_F2h_vs_hf'],
    })

V = pd.DataFrame(results)
print(f'\nprobes with xwu match: {len(V)}')
print(f'allele swaps: {(V["alleles_swap"]==True).sum()}, '
      f'genuinely-different alleles: {(V["alleles_swap"].isna()).sum()}')

# Summary
ok = V[V['alleles_swap'].notna()].copy()
print(f'\n=== 3-source AF comparison (n={len(ok)} probes with comparable alleles) ===')
print(f'  mean cactus_em            : {ok["cactus_em"].mean():.3f}')
print(f'  mean cactus_em + F2hon    : {ok["cactus_em_F2hon"].mean():.3f}')
print(f'  mean recipe_cactus        : {ok["recipe_cactus"].mean():.3f}')
print(f'  mean recipe_F2hon         : {ok["recipe_F2hon"].mean():.3f}')
print(f'  mean recipe_xwu (any-ALT) : {ok["recipe_xwu_any"].mean():.3f}')
print(f'  mean recipe_xwu (strict)  : {ok["recipe_xwu_strict"].mean():.3f}')
print(f'  mean hapfire_af           : {ok["hapfire_af"].mean():.3f}')

print(f'\nMAE vs hapFIRE:')
for col, name in [('cactus_em', 'cactus_em'),
                  ('cactus_em_F2hon', 'cactus_em + F2hon'),
                  ('recipe_cactus','recipe_cactus'),
                  ('recipe_F2hon', 'recipe_F2honest'),
                  ('recipe_xwu_any', 'recipe_xwu (any-ALT)'),
                  ('recipe_xwu_strict', 'recipe_xwu (strict hom-ALT)')]:
    if col not in ok.columns: continue
    r = ok['hapfire_af'] - ok[col]
    print(f'  {name:30s}  mean_bias={r.mean():+.4f}  MAE={r.abs().mean():.4f}')

print('\n=== Per-probe table (sorted by hapFIRE descending) ===')
ok_sorted = ok.sort_values('hapfire_af', ascending=False)
print(ok_sorted[['pos','ref','alt','cactus_em','cactus_em_F2hon',
                 'recipe_xwu_any','hapfire_af','AC_cactus','AN_cactus',
                 'AC_xwu_any','AN_xwu','n_override','n_sv_at_p']].round(3).to_string(index=False))

V.to_csv(ROOT / 'scratch/f2_honest/celld_3source_validation.tsv', sep='\t', index=False)
print(f'\nWrote celld_3source_validation.tsv')
