#!/usr/bin/env python3
"""Headless version of notebooks/MISSINGNESS_CAUSES_IMBALANCE.ipynb.

Answers "is the cactus/PG imbalance caused by missingness?" from on-disk data:
the imbalance inside fully-called bubbles (where N-on == N-off) is pure biology.
Prints all stats to stdout and saves the dose-response figure + a results npz.
"""
import json, sys, time
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.sparse import load_npz
import pysam

ROOT = Path('/global/scratch/users/tbellg/kmate')
CHROM = sys.argv[1] if len(sys.argv) > 1 else 'Chr1'
CN    = ROOT / f'data/kmer_pa_231_v3qc_v3/cn_{CHROM}.kmer_pa.npz'
META  = ROOT / f'data/kmer_pa_231_v3qc_v3/cn_{CHROM}.meta.npz'
VCF   = ROOT / 'panel/pangenie_genotyping/data/v3qc_v3/founders_231_v3qc_v3.haploid.vcf.gz'
SPLIT = ROOT / 'data/founder_split_cactus_pg.json'
PLOTDIR = ROOT/'notebooks/plots'; PLOTDIR.mkdir(exist_ok=True)
RESULTS = ROOT/'results'; RESULTS.mkdir(exist_ok=True)

t0 = time.time()
def log(*a): print(f'[{time.time()-t0:7.1f}s]', *a, flush=True)

kmer_pa   = load_npz(CN).tocsr()
meta = np.load(META, allow_pickle=True)
founders  = np.asarray(meta['founders']).astype(str)
F = len(founders)
bubble_id = np.asarray(meta['bubble_id']).astype(np.int64)
bub_start = np.asarray(meta['bubble_start']).astype(np.int64)
bub_end   = np.asarray(meta['bubble_end']).astype(np.int64)
n_bubbles = len(bub_start)

with open(SPLIT) as f: split = json.load(f)
cset, pset = set(map(str, split['cactus'])), set(map(str, split['PG']))
is_c = np.array([x in cset for x in founders]); is_p = np.array([x in pset for x in founders])
n_c, n_p = int(is_c.sum()), int(is_p.sum())
log(f'kmer_pa={kmer_pa.shape} nnz={kmer_pa.nnz:,} bubbles={n_bubbles:,} cactus={n_c} PG={n_p}')

ac   = np.asarray(kmer_pa.sum(0)).flatten()
ac_c = np.asarray(kmer_pa[is_c].sum(0)).flatten()
ac_p = np.asarray(kmer_pa[is_p].sum(0)).flatten()
cn_c = kmer_pa[is_c].tocsc()   # pre-slice once; column-subsetting is the hot op
cn_p = kmer_pa[is_p].tocsc()
log('carrier counts done')

# --- per-bubble missingness from the VCF (faithful per-bubble fetch) ---
vcf = pysam.VariantFile(str(VCF))
vsamples = list(vcf.header.samples)
vidx = np.array([vsamples.index(f) for f in founders])  # VCF col -> kmer_pa order
miss_bub = np.zeros((n_bubbles, F), dtype=bool)
for b in range(n_bubbles):
    s, e = int(bub_start[b]), int(bub_end[b])
    try:
        recs = vcf.fetch(CHROM, max(0, s-1), e)
    except Exception:
        continue
    m = np.zeros(len(vsamples), dtype=bool)
    has = False
    for rec in recs:
        has = True
        for j, smp in enumerate(rec.samples.values()):
            gt = smp['GT']
            if gt is None or len(gt) == 0 or gt[0] is None:
                m[j] = True
    if has:
        miss_bub[b] = m[vidx]
    if (b+1) % 20000 == 0:
        log(f'  VCF pass {b+1:,}/{n_bubbles:,}')
log('VCF pass done')

pg_miss_frac  = miss_bub[:, is_p].sum(1) / n_p
any_miss      = miss_bub.any(1)
log(f'fully-called bubbles (0 missing): {(~any_miss).sum():,} / {n_bubbles:,} '
    f'({(~any_miss).mean()*100:.1f}%)')

def priv_ratio(km):
    pmask = km & (ac == 1)
    c = cn_c[:, pmask].sum() / n_c; p = cn_p[:, pmask].sum() / n_p
    return c, p, (c/p if p else np.nan)
def side_ratio(km):
    s = km & (ac >= 2); co = s & (ac_p == 0); po = s & (ac_c == 0)
    c = cn_c[:, co].sum() / n_c; p = cn_p[:, po].sum() / n_p
    return c, p, (c/p if p else np.nan)

allk = np.ones(kmer_pa.shape[1], dtype=bool)
fully = ~any_miss[bubble_id]

print('\n================  HEADLINE: PRIVATE k-mers per founder  ================')
res = {}
for name, km in [('ALL bubbles', allk), ('FULLY-CALLED only', fully)]:
    c, p, r = priv_ratio(km); res[f'priv_{name}'] = (c, p, r)
    print(f'  {name:20s}: cactus/founder={c:8.1f}  PG/founder={p:7.1f}  ratio={r:5.2f}x')
print('================  SIDE-ONLY-SHARED k-mers per founder  ================')
for name, km in [('ALL bubbles', allk), ('FULLY-CALLED only', fully)]:
    c, p, r = side_ratio(km); res[f'side_{name}'] = (c, p, r)
    print(f'  {name:20s}: cactus/founder={c:8.1f}  PG/founder={p:7.1f}  ratio={r:5.2f}x')
print(f'\nac==0 fraction: all={ (ac==0).mean()*100:5.2f}%   fully-called={(ac[fully]==0).mean()*100:5.2f}%')

# --- dose-response ---
edges = [0.0, 1e-9, 0.02, 0.05, 0.10, 0.20, 0.40, 1.01]
labels = ['0%', '0-2%', '2-5%', '5-10%', '10-20%', '20-40%', '>40%']
strat = np.digitize(pg_miss_frac, edges[1:-1], right=False)
priv_r, side_r, nbub = [], [], []
print('\n================  DOSE-RESPONSE (ratio vs PG missingness)  ================')
for s in range(len(labels)):
    km = (strat == s)[bubble_id]; nb = int((strat == s).sum())
    pr = priv_ratio(km)[2]; sr = side_ratio(km)[2]
    priv_r.append(pr); side_r.append(sr); nbub.append(nb)
    print(f'  PG-miss {labels[s]:8s}: {nb:7,} bubbles  priv_ratio={pr:5.2f}x  side_ratio={sr:5.2f}x')

# --- attribution ---
priv = (ac == 1); cac_priv = priv & (ac_c == 1); pg_priv = priv & (ac_p == 1)
pgmiss_km = (pg_miss_frac[bubble_id] > 0)
print('\n================  ATTRIBUTION  ================')
print(f'  cactus private k-mers: {cac_priv.sum():,}; in >0-PG-missing bubbles: '
      f'{(cac_priv & pgmiss_km).sum():,} ({(cac_priv & pgmiss_km).sum()/max(cac_priv.sum(),1)*100:.1f}%)')
print(f'  PG     private k-mers: {pg_priv.sum():,}; in >0-PG-missing bubbles: '
      f'{(pg_priv & pgmiss_km).sum():,} ({(pg_priv & pgmiss_km).sum()/max(pg_priv.sum(),1)*100:.1f}%)')
print(f'  reference: {(pg_miss_frac>0).mean()*100:.1f}% of bubbles have >0 PG missing')

# --- figure ---
fig, ax = plt.subplots(figsize=(10, 5.5))
x = np.arange(len(labels))
ax.plot(x, priv_r, '-o', color='#e31a1c', label='private ratio (cactus/PG per founder)')
ax.plot(x, side_r, '-s', color='#1f78b4', label='side-only-shared ratio')
ax.axhline(1.0, color='k', ls='--', lw=1, label='1x = balanced')
for xi, n in zip(x, nbub): ax.annotate(f'{n:,}', (xi, 0.2), ha='center', fontsize=7, color='gray')
ax.set_xticks(x); ax.set_xticklabels(labels)
ax.set_xlabel('PG missingness in bubble'); ax.set_ylabel('cactus/PG ratio')
ax.set_title(f'{CHROM}: does the imbalance grow with PG missingness?')
ax.legend()
fig.tight_layout()
fig.savefig(PLOTDIR/f'missingness_doseresponse_{CHROM}.png', dpi=130)
log(f'saved {PLOTDIR}/missingness_doseresponse_{CHROM}.png')

np.savez(RESULTS/f'missingness_test_{CHROM}.npz',
         labels=np.array(labels), priv_r=np.array(priv_r), side_r=np.array(side_r),
         nbub=np.array(nbub), pg_miss_frac=pg_miss_frac, any_miss=any_miss)
log('DONE')
