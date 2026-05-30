#!/usr/bin/env python3
"""DIAGNOSTIC: compare the cactus/PG k-mer imbalance between the N-on production
matrix and the N-off diagnostic rebuild, on one chromosome.

This is the airtight confirmation for "is the imbalance caused by missingness?":
if the all-bubble private ratio drops sharply N-on -> N-off, the imbalance is
missingness-driven; if it barely moves, it is real biology. The N-off all-bubble
ratio should also roughly match the fully-called-bubble estimate from
notebooks/MISSINGNESS_CAUSES_IMBALANCE.ipynb.

Usage:
  python scripts/compare_Non_vs_Noff_imbalance.py [Chr1]
Needs: numpy, scipy. Env with the old 'hapfm' gone -- use any env that has them.
"""
import sys, json
from pathlib import Path
import numpy as np
from scipy.sparse import load_npz

ROOT  = Path('/global/scratch/users/tbellg/kmate')
CHROM = sys.argv[1] if len(sys.argv) > 1 else 'Chr1'
NON   = ROOT / f'data/kmer_pa_231_v3qc_v3/cn_{CHROM}.kmer_pa.npz'
NOFF  = ROOT / f'data/kmer_pa_231_v3qc_v3_Noff_diag/cn_{CHROM}.kmer_pa.npz'
META  = ROOT / f'data/kmer_pa_231_v3qc_v3/cn_{CHROM}.meta.npz'   # founder order (same for both)
SPLIT = ROOT / 'data/founder_split_cactus_pg.json'

meta = np.load(META, allow_pickle=True)
founders = np.asarray(meta['founders']).astype(str)
with open(SPLIT) as f: split = json.load(f)
cset, pset = set(map(str, split['cactus'])), set(map(str, split['PG']))
is_c = np.array([x in cset for x in founders]); is_p = np.array([x in pset for x in founders])
n_c, n_p = int(is_c.sum()), int(is_p.sum())

def stats(path, label):
    if not path.exists():
        print(f'[{label}] MISSING: {path}  (run the rebuild first)'); return
    kmer_pa = load_npz(path).tocsr()
    ac   = np.asarray(kmer_pa.sum(0)).flatten()
    ac_c = np.asarray(kmer_pa[is_c].sum(0)).flatten()
    ac_p = np.asarray(kmer_pa[is_p].sum(0)).flatten()
    priv = (ac == 1)
    cac_priv = kmer_pa[is_c][:, priv].sum() / n_c
    pg_priv  = kmer_pa[is_p][:, priv].sum() / n_p
    s = (ac >= 2); co = s & (ac_p == 0); po = s & (ac_c == 0)
    cac_so = kmer_pa[is_c][:, co].sum() / n_c
    pg_so  = kmer_pa[is_p][:, po].sum() / n_p
    print(f'\n===== {label}  ({path.parent.name}) =====')
    print(f'  shape={kmer_pa.shape} nnz={kmer_pa.nnz:,}')
    print(f'  ac==0 fraction:            {(ac==0).mean()*100:6.2f}%')
    print(f'  ac==1 (private) fraction:  {(ac==1).mean()*100:6.2f}%')
    print(f'  private/founder  cactus={cac_priv:8.1f}  PG={pg_priv:7.1f}  ratio={cac_priv/max(pg_priv,1e-9):5.2f}x')
    print(f'  side-only/founder cactus={cac_so:8.1f}  PG={pg_so:7.1f}  ratio={cac_so/max(pg_so,1e-9):5.2f}x')

print(f'chrom={CHROM}  cactus={n_c}  PG={n_p}')
stats(NON,  'N-ON  (production)')
stats(NOFF, 'N-OFF (diagnostic)')
print('\nInterpretation: large drop in the private ratio N-on->N-off => imbalance is '
      'missingness-driven. Little change => real biology.')
