#!/usr/bin/env python3
"""
3-source validation at cell [d] positions.

For each cell [d] SNP position, compare per-founder calls across:
  1. CACTUS path-derived  : from raw deconstruct VCF (biallelic, with F2-honest aggregation)
                            Plus the merged 231-panel cn_var (current production)
  2. XWU GrENE-Net SNP catalog : per-founder short-read SNP calls (hapFIRE's source)
  3. (Tie-breaker if needed) direct lookup in founder long-read FASTAs

Compute:
  recipe_cactus  = AC_cactus / AN_cactus on the 231-panel cn_var (current)
  recipe_F2hon   = recipe_cactus + cactus-F2-honest overrides
  recipe_xwu     = AC_xwu / AN_xwu from xwu's catalog at the same positions
  hapfire_af     = hapFIRE's projected AF (from s1_snp_frequency.txt)

If recipe_xwu ≈ hapfire_af AND recipe_xwu ≠ recipe_cactus → it's a calls-source disagreement.
If recipe_xwu ≈ recipe_F2hon AND both ≠ hapfire_af → hapFIRE is adding something on top of xwu.
"""
import numpy as np
import pandas as pd
import subprocess
from pathlib import Path

ROOT = Path('/carnegie/nobackup/scratch/tbellagio/hapfire_sv')
BCF = '/home/tbellagio/miniforge3/envs/sequencing_pipeline/bin/bcftools'
XWU = '/carnegie/nobackup/scratch/xwu/GrENE_net/hapFIRE_updatedVCF/greneNet_final_v1.1_chr1.recode.vcf'

# ---------- 1. Pick 30 cell [d] probes ----------
# Use the F2-honest joined table; pick records with high hapFIRE-cactus_em residual
# and a range of n_override values (some where F2-honest reached, some not).
M = pd.read_csv(ROOT / 'scratch/f2_honest/f2_honest_test.tsv', sep='\t')
M['cell'] = (M['high_fmiss'].astype(int)*2 + M['is_mixed_bubble']).map({
    0:'a', 1:'b', 2:'c', 3:'d'})
D = M[M['cell']=='d'].copy()
D['resid_orig'] = D['hapfire_af'] - D['alt_freq']
D['resid_F2h']  = D['hapfire_af'] - D['cactus_em_F2honest']
print(f'cell [d] candidates: {len(D):,}')

# Sample 30 records: 10 with highest residual after F2-honest, 10 mid, 10 low-residual
D_sorted = D.sort_values('resid_F2h', ascending=False)
probes = pd.concat([
    D_sorted.head(10),
    D_sorted.iloc[len(D_sorted)//2 - 5 : len(D_sorted)//2 + 5],
    D_sorted.tail(10)
]).drop_duplicates(subset=['pos']).head(30)
print(f'probes selected: {len(probes)}')
print(probes[['pos','ref','alt','alt_freq','cactus_em_F2honest','hapfire_af',
              'recipe_af','n_override','n_sv_carriers_at_p','AC_post','AN_post',
              'resid_orig','resid_F2h']].to_string())

# ---------- 2. Look up each probe in xwu's GrENE-Net SNP VCF ----------
print('\n[2] Looking up probes in xwu GrENE-Net SNP catalog')
# Get xwu sample list
out = subprocess.run([BCF, 'view', '-h', XWU], capture_output=True, text=True)
xwu_header = [l for l in out.stdout.split('\n') if l.startswith('#CHROM')][0]
xwu_samples = xwu_header.rstrip().split('\t')[9:]
xwu_samples = [s for s in xwu_samples if s]   # drop blank trailing column
print(f'  xwu n founders: {len(xwu_samples)}')

# Build query: extract per-founder GT at each probe pos
# Use regions argument for efficiency
regions = ','.join([f'1:{int(p)}-{int(p)}' for p in probes['pos']])

# Query: pos, ref, alt, then per-sample GT
fmt = '%CHROM\\t%POS\\t%REF\\t%ALT[\\t%GT]\\n'
out = subprocess.run(
    [BCF, 'query', '-r', regions, '-f', fmt, XWU],
    capture_output=True, text=True
)
xwu_rows = [l for l in out.stdout.split('\n') if l]
print(f'  xwu records matching probe positions: {len(xwu_rows)}')

# Parse into per-founder GT dict
xwu_data = {}
for line in xwu_rows:
    fields = line.split('\t')
    pos = int(fields[1]); ref = fields[2]; alt = fields[3]
    gts = fields[4:]
    xwu_data[(pos, ref, alt)] = gts

# ---------- 3. Compute recipe_xwu and per-position comparison ----------
print('\n[3] Compare per-position AFs across 3 sources')
results = []
for _, p in probes.iterrows():
    pos = int(p['pos']); ref = str(p['ref']); alt = str(p['alt'])
    xwu_key = (pos, ref, alt)
    xwu_gts = xwu_data.get(xwu_key)
    if xwu_gts is None:
        # Maybe REF/ALT swapped or position-only match
        match_pos = [k for k in xwu_data if k[0]==pos]
        if match_pos:
            # Try a variant with any ref/alt
            xwu_key = match_pos[0]
            xwu_gts = xwu_data[xwu_key]
            ref_xwu, alt_xwu = xwu_key[1], xwu_key[2]
            # Flag if REF/ALT mismatch with cn_var
            swap = (ref == alt_xwu) and (alt == ref_xwu)
        else:
            xwu_gts = None
    else:
        ref_xwu, alt_xwu = ref, alt
        swap = False
    if xwu_gts is None:
        ac_xwu = an_xwu = 0
        recipe_xwu = float('nan')
    else:
        # GT can be "0", "1", "0/0", "1/1", "0/1", "./.", "."
        ac = an = 0
        for g in xwu_gts:
            if g in ('.', './.', '.|.'):
                continue
            an += 1
            # carrier if any "1"
            if '1' in g.split('/')[0] or (len(g.split('/'))>1 and '1' in g.split('/')[1]):
                ac += 1
            elif '|' in g:
                if '1' in g.split('|')[0] or '1' in g.split('|')[1]:
                    ac += 1
            elif g == '1':
                ac += 1
        ac_xwu = ac
        an_xwu = an
        recipe_xwu = ac/an if an>0 else float('nan')
        if swap:
            recipe_xwu = 1 - recipe_xwu
    results.append({
        'pos': pos, 'ref': ref, 'alt': alt,
        'cactus_em':       p['alt_freq'],
        'cactus_em_F2hon': p['cactus_em_F2honest'],
        'recipe_cactus':   p['recipe_af'],
        'recipe_F2hon':    p['recipe_F2honest'],
        'recipe_xwu':      recipe_xwu,
        'hapfire_af':      p['hapfire_af'],
        'AC_cactus':       int(p['AC_post']),
        'AN_cactus':       int(p['AN_post']),
        'AC_xwu':          ac_xwu,
        'AN_xwu':          an_xwu,
        'n_override':      int(p['n_override']),
        'n_sv_at_p':       int(p['n_sv_carriers_at_p']),
        'resid_orig_vs_hf': p['resid_orig'],
        'resid_F2h_vs_hf':  p['resid_F2h'],
        'xwu_in_catalog':  xwu_gts is not None,
    })

R = pd.DataFrame(results)
print('\n3-source comparison on 30 cell [d] probes (sorted by hapFIRE):')
print(R.sort_values('hapfire_af', ascending=False)[
    ['pos','ref','alt','cactus_em','cactus_em_F2hon','recipe_cactus','recipe_F2hon',
     'recipe_xwu','hapfire_af','AC_cactus','AN_cactus','AC_xwu','AN_xwu',
     'n_override','n_sv_at_p']
].to_string(index=False))

# Summary stats
print('\n=== Summary: how do the 3 estimators agree on cell [d] probes? ===')
valid = R[R['xwu_in_catalog'] & R['recipe_xwu'].notna()]
print(f'  probes with xwu data: {len(valid)} of {len(R)}')
print(f'  mean recipe_cactus     : {valid["recipe_cactus"].mean():.3f}')
print(f'  mean recipe_F2honest   : {valid["recipe_F2hon"].mean():.3f}')
print(f'  mean recipe_xwu        : {valid["recipe_xwu"].mean():.3f}')
print(f'  mean hapfire_af        : {valid["hapfire_af"].mean():.3f}')
print()
print(f'  |recipe_xwu - hapfire|  mean: {(valid["recipe_xwu"]-valid["hapfire_af"]).abs().mean():.4f}')
print(f'  |recipe_xwu - recipe_cactus| mean: {(valid["recipe_xwu"]-valid["recipe_cactus"]).abs().mean():.4f}')
print(f'  |recipe_xwu - recipe_F2hon| mean: {(valid["recipe_xwu"]-valid["recipe_F2hon"]).abs().mean():.4f}')
print(f'  |recipe_cactus - hapfire|  mean: {(valid["recipe_cactus"]-valid["hapfire_af"]).abs().mean():.4f}')

# Save
R.to_csv(ROOT / 'scratch/f2_honest/celld_validation.tsv', sep='\t', index=False)
print(f'\nWrote {ROOT / "scratch/f2_honest/celld_validation.tsv"}')
