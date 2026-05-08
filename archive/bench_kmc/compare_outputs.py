"""Verify jellyfish and KMC produce the same per-k-mer counts.

Reads:
  - jellyfish output: lines of "<kmer> <count>"
  - KMC dump output:  lines of "<kmer>\t<count>"
  - query.fa:         all 20.8M k-mers (reference set)

Both tools use canonical k-mers; rep may differ (jellyfish picks lex-min of
forward/revcomp; KMC also). For each query k-mer we compare the count from
each tool. Reports: # k-mers seen by each, # disagreements, count distribution.
"""
import argparse, sys
from collections import defaultdict


def revcomp(s):
    comp = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
    try:
        return ''.join(comp[c] for c in s[::-1])
    except KeyError:
        return None


def canon(s):
    s = s.upper()
    rc = revcomp(s)
    if rc is None:
        return None
    return s if s <= rc else rc


def load_jf(path):
    """jellyfish query output: '<kmer> <count>'."""
    counts = {}
    with open(path) as f:
        for line in f:
            parts = line.split()
            if len(parts) >= 2:
                c = canon(parts[0])
                if c is not None:
                    counts[c] = int(parts[1])
    return counts


def load_kmc(path):
    """KMC dump output: '<kmer>\t<count>'."""
    counts = {}
    with open(path) as f:
        for line in f:
            parts = line.rstrip('\n').split('\t')
            if len(parts) >= 2:
                c = canon(parts[0])
                if c is not None:
                    counts[c] = int(parts[1])
    return counts


def load_query_kmers(path):
    kmers = []
    n_skip = 0
    with open(path) as f:
        for line in f:
            if line.startswith('>'): continue
            c = canon(line.strip())
            if c is None or len(c) != 31:
                n_skip += 1
                continue
            kmers.append(c)
    if n_skip:
        print(f'  skipped {n_skip} invalid query k-mers (nan/N/non-ATCG)')
    return kmers


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--jf', required=True)
    ap.add_argument('--kmc', required=True)
    ap.add_argument('--query', required=True)
    args = ap.parse_args()

    print(f'loading jellyfish counts from {args.jf}', flush=True)
    jf = load_jf(args.jf)
    print(f'  {len(jf):,} unique k-mers, sum {sum(jf.values()):,}')
    print(f'loading KMC counts from {args.kmc}', flush=True)
    kmc = load_kmc(args.kmc)
    print(f'  {len(kmc):,} unique k-mers, sum {sum(kmc.values()):,}')
    print(f'loading query set from {args.query}', flush=True)
    query = load_query_kmers(args.query)
    print(f'  {len(query):,} query k-mers')

    n_total = len(query)
    n_match = 0
    n_disagree = 0
    n_jf_only = 0
    n_kmc_only = 0
    n_both_zero = 0
    disagreements = []
    for k in query:
        jc = jf.get(k, 0)
        kc = kmc.get(k, 0)
        if jc == kc:
            n_match += 1
            if jc == 0:
                n_both_zero += 1
        else:
            n_disagree += 1
            if jc > 0 and kc == 0:
                n_jf_only += 1
            elif kc > 0 and jc == 0:
                n_kmc_only += 1
            if len(disagreements) < 10:
                disagreements.append((k, jc, kc))

    print()
    print(f'=== Per-query-kmer comparison (n={n_total:,}) ===')
    print(f'  agree:           {n_match:,}  ({100*n_match/n_total:.4f}%)')
    print(f'    of which both zero: {n_both_zero:,}')
    print(f'    of which both >0:   {n_match - n_both_zero:,}')
    print(f'  disagree:        {n_disagree:,}  ({100*n_disagree/n_total:.4f}%)')
    if n_disagree:
        print(f'    jf>0 kmc=0:    {n_jf_only:,}')
        print(f'    kmc>0 jf=0:    {n_kmc_only:,}')
        print(f'    other:         {n_disagree - n_jf_only - n_kmc_only:,}')
        print('  first 10 disagreements (kmer, jf, kmc):')
        for k, jc, kc in disagreements:
            print(f'    {k}: jf={jc}  kmc={kc}')


if __name__ == '__main__':
    main()
