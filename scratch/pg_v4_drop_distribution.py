"""Replicate xwu's exact het-filter metric on the PG-side of v3.

xwu's filter (R script):
  V4 = AC_Het / AC                           # fraction of ALT alleles in het samples
  keep records where V4 < 0.01

Compute V4 per record on the 151 PG founders' GTs from
founders_231_chr.vcf.gz. Show drop counts at xwu's exact threshold + others.
"""
from __future__ import annotations
import json, time
from pathlib import Path
import numpy as np
import pysam

ROOT = Path('/carnegie/nobackup/scratch/tbellagio/hapfire_sv')
SRC = ROOT / 'pangenie_genotyping/data/merged/founders_231_chr.vcf.gz'
SPLIT = json.load(open(ROOT/'data/founder_split_cactus_pg.json'))
pg_set = set(SPLIT['PG'])

vcf = pysam.VariantFile(str(SRC))
samples = list(vcf.header.samples)
pg_idx = [i for i, s in enumerate(samples) if s in pg_set]
print(f'PG founders: {len(pg_idx)}', flush=True)

t = time.time()
ac_all, ac_het_all, ref_lens, alt_lens, poss = [], [], [], [], []
n = 0
for rec in vcf.fetch('Chr1'):
    if rec.alts is None or len(rec.alts) != 1:
        continue
    ac = 0          # total ALT alleles in PG founders
    ac_het = 0      # ALT alleles contributed by het GTs (= # of het samples)
    for i in pg_idx:
        gt = rec.samples[samples[i]]['GT']
        if gt is None or len(gt) == 0 or all(a is None for a in gt):
            continue
        nonnull = [a for a in gt if a is not None]
        n_alt = sum(1 for a in nonnull if a > 0)
        ac += n_alt
        if len(nonnull) >= 2 and 0 in nonnull and any(a > 0 for a in nonnull):
            ac_het += 1   # each het sample contributes 1 alt allele to AC_Het
    ac_all.append(ac); ac_het_all.append(ac_het)
    ref_lens.append(len(rec.ref)); alt_lens.append(len(rec.alts[0]))
    poss.append(rec.pos)
    n += 1
    if n % 200000 == 0:
        print(f'  parsed {n:,} ({time.time()-t:.0f}s)', flush=True)

ac_all = np.array(ac_all, dtype=np.int32)
ac_het_all = np.array(ac_het_all, dtype=np.int32)
ref_lens = np.array(ref_lens, dtype=np.int32)
alt_lens = np.array(alt_lens, dtype=np.int32)
poss = np.array(poss, dtype=np.int64)

# V4 = AC_Het / AC; if AC == 0, set V4 = 0 (no ALT, nothing to filter on)
v4 = np.zeros_like(ac_all, dtype=np.float32)
nz = ac_all > 0
v4[nz] = ac_het_all[nz] / ac_all[nz]

print(f'\nTotal Chr1 biallelic records: {len(v4):,}, {time.time()-t:.0f}s')
print(f'Records with AC > 0: {nz.sum():,}')
print(f'Records with AC = 0: {(~nz).sum():,}')

print(f'\n=== xwu V4 (AC_Het / AC) distribution on records with AC > 0 ===')
v4_nz = v4[nz]
print(f'  mean   : {v4_nz.mean():.5f}')
print(f'  median : {np.median(v4_nz):.5f}')
for p in [25, 50, 75, 90, 95, 99, 99.5, 99.9]:
    print(f'  P{p:5.1f}  : {np.percentile(v4_nz, p):.5f}')

print(f'\n=== Drop counts at xwu-like V4 thresholds ===')
total = len(v4)
print(f'{"keep V4 <":>12s}  {"drop":>10s}  {"drop %":>8s}  {"keep":>10s}')
for thr in [0.001, 0.005, 0.01, 0.02, 0.05, 0.10, 0.20, 0.50, 1.01]:
    drop = ((v4 >= thr) & nz).sum()  # AC=0 records are kept by default
    keep = total - drop
    flag = '  ← xwu' if thr == 0.01 else ''
    print(f'  V4 < {thr:>5.3f}    {drop:>10,}  {100*drop/total:>7.2f}%  {keep:>10,}{flag}')

# Stratify SNPs vs INDELs
is_snp = (ref_lens == 1) & (alt_lens == 1)
print(f'\n=== Stratified by variant type, V4 < 0.01 (xwu threshold) ===')
for label, mask in [('SNPs only', is_snp), ('INDELs/SVs', ~is_snp)]:
    drop = ((v4 >= 0.01) & nz & mask).sum()
    print(f'  {label} ({mask.sum():,} records): drop {drop:>9,} ({100*drop/max(mask.sum(),1):5.2f}%)')

# Save for plot
np.savez(ROOT/'scratch/pg_v4_chr1.npz', v4=v4, ac=ac_all, ac_het=ac_het_all,
         ref_len=ref_lens, alt_len=alt_lens, pos=poss)
print(f'\nsaved scratch/pg_v4_chr1.npz')
