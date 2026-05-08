"""Dump Chr1's 20.8M cn_full_231_v2 query k-mers to a FASTA for benchmarking."""
import numpy as np
META = '/home/tbellagio/scratch/hapfire_sv/poolfreq/data/cn_full_231_v2/cn_Chr1.meta.npz'
OUT = '/home/tbellagio/scratch/hapfire_sv/poolfreq/bench_kmc/query_chr1.fa'

m = np.load(META, allow_pickle=True)
kmers = list(m['kmer_index'])
print(f'writing {len(kmers):,} k-mers to {OUT}')
with open(OUT, 'w') as f:
    for i, k in enumerate(kmers):
        f.write(f'>k{i}\n{k}\n')
print('done')
