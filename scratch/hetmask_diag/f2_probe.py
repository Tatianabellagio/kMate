#!/usr/bin/env python3
"""
Probe: pick a few high-leverage mixed-bubble SNP records (from cell [d]),
look up the raw VCF at the same (chrom, pos) and inspect:

  - Are there LV>=1 records at the same coords as our SNP record?
  - Are the GTs for cactus founders at the LV>=1 SNP record DIFFERENT from
    their GTs at the LV=0 multi-allelic (which is what cn_var encodes)?
"""
import subprocess
import pandas as pd
from pathlib import Path
import gzip

ROOT = Path('/carnegie/nobackup/scratch/tbellagio/hapfire_sv')
RAW  = '/home/tbellagio/scratch/pang/pang_1001gplus/pang/output/pang_1001gplus_82acc.raw.vcf'

# Load F1 test joined table; pick cell [d] (high F_MISSING + mixed bubble) records
F1 = pd.read_csv(ROOT / 'scratch/hetmask_diag/f1_test_joined.tsv', sep='\t')
D = F1[(F1['high_fmiss']) & (F1['is_mixed_bubble']==1)].copy()
# Sort by n_sv_carriers descending — high-leverage records first
D = D.sort_values('n_sv_carriers_at_p', ascending=False)
print(f'cell [d] records: {len(D):,}')

# Sample 5 high-leverage records on Chr1, spanning n_sv_carriers
probes = D.iloc[:5][['chrom','pos','ref','alt','alt_freq','hapfire_af',
                     'recipe_af','n_sv_carriers_at_p','AC_post','AN_post']]
print('\nProbe records:')
print(probes.to_string(index=False))

# Build sample list from raw VCF header
print('\n--- Raw VCF cactus sample list ---')
with open(RAW) as f:
    for line in f:
        if line.startswith('#CHROM'):
            cactus_samples = line.rstrip().split('\t')[9:]
            break
print(f'  n cactus samples: {len(cactus_samples)}')
print(f'  first 5: {cactus_samples[:5]}')

# For each probe record, scan the raw VCF at the chrom and find any records at that pos
print('\n--- Lookups in raw VCF (linear scan; positions are sorted within chrom) ---')
for _, row in probes.iterrows():
    chrom = row['chrom']; p = int(row['pos']); ref = row['ref']; alt = row['alt']
    print(f'\n== probe Chr1:{p}  REF={ref}  ALT={alt}  '
          f'(cactus_em={row["alt_freq"]:.3f}, hapFIRE={row["hapfire_af"]:.3f}, '
          f'n_sv_carriers_at_p={row["n_sv_carriers_at_p"]})')
    # Linear scan only over Chr1 lines; stop after we pass the pos
    matches = []
    with open(RAW) as f:
        in_chr1 = False
        for line in f:
            if line.startswith('#'):
                continue
            f0 = line.split('\t', 9)
            ch  = f0[0]; po = int(f0[1])
            if ch != chrom:
                if in_chr1:
                    break
                continue
            in_chr1 = True
            if po > p:
                break
            if po == p:
                matches.append(line.rstrip())
    print(f'   raw VCF rows at Chr1:{p}: {len(matches)}')
    for ln in matches:
        f0 = ln.split('\t')
        info = f0[7]
        lv = next((s for s in info.split(';') if s.startswith('LV=')), 'LV=?')
        ps = next((s for s in info.split(';') if s.startswith('PS=')), '')
        ac = next((s for s in info.split(';') if s.startswith('AC=')), '')
        print(f'     pos={f0[1]:>10} REF={f0[3][:10]} ALT={f0[4][:30]} {lv} {ps} {ac}')
