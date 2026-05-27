"""Compute per-record PG-side heterozygosity rate to determine HWE excess-het
filter thresholds for v3. Iterate the merged founders_231_chr.vcf.gz (Chr1
only for tractability), compute het rate among the 151 PG founders per
record, then report drop counts at various thresholds.
"""
from __future__ import annotations
import json, time
from pathlib import Path
import numpy as np
import pysam

ROOT = Path('/carnegie/nobackup/scratch/tbellagio/hapfire_sv')
SRC_VCF = ROOT / 'pangenie_genotyping/data/merged/founders_231_chr.vcf.gz'
SAMPLES_FILE = ROOT / 'data/vcf_samples_231.txt'
SPLIT_FILE = ROOT / 'data/founder_split_cactus_pg.json'

founders_order = [l.strip() for l in open(SAMPLES_FILE)]
F = len(founders_order)
split = json.load(open(SPLIT_FILE))
pg_set = set(split['PG'])
n_pg = len(pg_set)
print(f'PG founders to scan: {n_pg}', flush=True)

vcf = pysam.VariantFile(str(SRC_VCF))
samples = list(vcf.header.samples)
pg_indices_in_vcf = np.array([i for i, s in enumerate(samples) if s in pg_set])
print(f'PG founders in VCF: {len(pg_indices_in_vcf)}', flush=True)

# Iterate Chr1 records and compute PG het rate
het_rates = []
n_called   = []     # how many PG founders had a non-missing GT
n_carriers = []     # how many PG founders had any ALT (het OR hom-alt)
ref_lens   = []
alt_lens   = []
poss       = []

t0 = time.time()
n = 0
for rec in vcf.fetch('Chr1'):
    if rec.alts is None or len(rec.alts) != 1:
        continue   # skip multi-allelic and indels-only records
    n_het = 0
    n_carr = 0
    n_call = 0
    for i in pg_indices_in_vcf:
        gt = rec.samples[samples[i]]['GT']
        if gt is None or len(gt) == 0 or all(a is None for a in gt):
            continue
        n_call += 1
        nonnull = [a for a in gt if a is not None]
        if any(a > 0 for a in nonnull):
            n_carr += 1
        # het: at least one ALT and at least one REF allele
        if len(nonnull) == 2 and nonnull[0] != nonnull[1]:
            n_het += 1
        elif len(nonnull) > 2 and len(set(nonnull)) > 1:
            n_het += 1
    rate = n_het / max(n_call, 1)
    het_rates.append(rate)
    n_called.append(n_call)
    n_carriers.append(n_carr)
    ref_lens.append(len(rec.ref))
    alt_lens.append(len(rec.alts[0]))
    poss.append(rec.pos)
    n += 1
    if n % 200000 == 0:
        print(f'  parsed {n:,} records ({time.time()-t0:.0f}s)', flush=True)

het_rates = np.array(het_rates, dtype=np.float32)
n_called = np.array(n_called, dtype=np.int16)
n_carriers = np.array(n_carriers, dtype=np.int16)
ref_lens = np.array(ref_lens, dtype=np.int32)
alt_lens = np.array(alt_lens, dtype=np.int32)
poss = np.array(poss, dtype=np.int64)

print(f'\nTotal Chr1 biallelic records: {len(het_rates):,}, {time.time()-t0:.0f}s')

# Save for plot/notebook
np.savez(ROOT/'scratch/pg_het_rate_chr1.npz',
         het_rate=het_rates, n_called=n_called, n_carriers=n_carriers,
         ref_len=ref_lens, alt_len=alt_lens, pos=poss)
print(f'saved scratch/pg_het_rate_chr1.npz')

# Quick stats
print(f'\nPG het rate distribution:')
print(f'  mean   : {het_rates.mean():.5f}')
print(f'  median : {np.median(het_rates):.5f}')
for p in [50, 75, 90, 95, 99, 99.5, 99.9]:
    print(f'  P{p:5.1f}  : {np.percentile(het_rates, p):.5f}')

# How many would be dropped at various thresholds?
print('\n=== PG het-rate drop counts at various thresholds (lower = more aggressive) ===')
print(f'{"threshold":>12s}  {"drop":>10s}  {"drop %":>8s}  {"keep":>10s}')
total = len(het_rates)
for t in [0.005, 0.01, 0.02, 0.03, 0.05, 0.10, 0.20, 0.50]:
    drop = (het_rates > t).sum()
    keep = total - drop
    print(f'  > {t:>5.3f}    {drop:>10,}  {100*drop/total:>7.2f}%  {keep:>10,}')

# Stratified by SNP vs indel
print('\n=== Stratified: SNPs vs INDELs ===')
is_snp = (ref_lens == 1) & (alt_lens == 1)
is_indel = ~is_snp
for label, mask in [('SNPs only', is_snp), ('INDELs/SVs', is_indel)]:
    sub = het_rates[mask]
    print(f'\n  {label} ({mask.sum():,} records):')
    print(f'    mean het: {sub.mean():.5f}')
    for t in [0.01, 0.02, 0.05, 0.10]:
        drop = (sub > t).sum()
        print(f'    >{t:.2f}: drop {drop:>9,} ({100*drop/max(mask.sum(),1):5.2f}%)')

# Sanity: how many records have ANY het call at all?
any_het = (het_rates > 0).sum()
print(f'\nRecords with ≥1 PG het call: {any_het:,} ({100*any_het/total:.2f}%)')
print(f'Records with PG het = 0:    {(het_rates == 0).sum():,}')
