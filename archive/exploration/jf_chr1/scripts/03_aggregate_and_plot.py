#!/usr/bin/env python
"""Aggregate per-anchor pair TSVs into an 82x82 Jaccard matrix and plot heatmap.

Reads <project_root>/jf_chr1/pairs/*.tsv
Writes:
  - jaccard_matrix.tsv  (82x82, k=31 canonical Chr1 Jaccard)
  - jaccard_heatmap.png
"""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import squareform

ROOT = Path(__file__).resolve().parents[1]

# Load sample order
samples = [l.strip() for l in (ROOT/'scripts/asm_ids.txt').read_text().splitlines() if l.strip()]
n = len(samples)
idx = {s: i for i, s in enumerate(samples)}
print(f'{n} samples')

# Read all per-anchor TSVs and build the symmetric matrix
J = np.full((n, n), np.nan, dtype=np.float64)
np.fill_diagonal(J, 1.0)

rows = []
for tsv in sorted((ROOT/'pairs').glob('*.tsv')):
    df = pd.read_csv(tsv, sep='\t')
    if len(df) == 0: continue
    for _, r in df.iterrows():
        i = idx[r.A]; j = idx[r.B]
        J[i, j] = r.jaccard
        J[j, i] = r.jaccard
        rows.append((r.A, r.B, r.nA, r.nB, r.nU, r.nI, r.jaccard))

n_pairs = len(rows)
expected = n*(n-1)//2
print(f'rows read: {n_pairs} (expected {expected})')
missing = np.isnan(J)
np.fill_diagonal(missing, False)
print(f'missing cells: {int(missing.sum())}')

# Save matrix as TSV
out_df = pd.DataFrame(J, index=samples, columns=samples)
out_df.to_csv(ROOT/'jaccard_matrix.tsv', sep='\t', float_format='%.6f')
print(f'wrote {ROOT}/jaccard_matrix.tsv')

# Also save long-form for the notebook
long_df = pd.DataFrame(rows, columns=['A','B','nA','nB','nU','nI','jaccard'])
long_df.to_csv(ROOT/'jaccard_pairs.tsv', sep='\t', index=False, float_format='%.6f')
print(f'wrote {ROOT}/jaccard_pairs.tsv ({len(long_df)} pairs)')

# Hierarchical clustering for sensible ordering
D = 1.0 - np.nan_to_num(J, nan=0.5)  # distance
np.fill_diagonal(D, 0)
condensed = squareform(D, checks=False)
Z = linkage(condensed, method='average')
order = leaves_list(Z)
J_ord = J[order][:, order]
names_ord = np.array(samples)[order]

# Plot
fig, ax = plt.subplots(figsize=(13, 12))
im = ax.imshow(J_ord, cmap='magma', aspect='equal', vmin=0.5, vmax=1.0)
cbar = plt.colorbar(im, ax=ax, fraction=0.04, pad=0.02)
cbar.set_label('Chr1 k-mer Jaccard (k=31, canonical)\nbright = high similarity', fontsize=11)

# Annotate suspect pairs
suspect_pairs = [('101003','100954'), ('100300','100692')]
labels = {'101003':'Set-1 (5772)', '100954':'T980 (6150)',
          '100300':'Ei-2 (6915)',  '100692':'St-0 (8387)'}
for a, b in suspect_pairs:
    ia = np.where(names_ord == a)[0][0]
    ib = np.where(names_ord == b)[0][0]
    for x, y in [(ia, ib), (ib, ia)]:
        ax.add_patch(plt.Rectangle((x-0.5, y-0.5), 1, 1, fill=False, edgecolor='cyan', lw=2.5))

ax.set_xticks(range(n))
ax.set_yticks(range(n))
ax.set_xticklabels(names_ord, rotation=90, fontsize=6)
ax.set_yticklabels(names_ord, fontsize=6)
ax.set_title(f'Pairwise Chr1 k-mer Jaccard across all {n} cactus pangenome assemblies\n'
             f'(jellyfish k=31, computed directly on raw FASTAs from chr_only/)\n'
             f'Two suspect duplicate pairs highlighted in cyan', fontsize=12)
plt.tight_layout()
plt.savefig(ROOT/'jaccard_heatmap.png', dpi=130, bbox_inches='tight')
print(f'wrote {ROOT}/jaccard_heatmap.png')

# Summary stats
ut = J[np.triu_indices_from(J, k=1)]
ut = ut[~np.isnan(ut)]
print('\n=== Pairwise Jaccard summary ({} pairs) ==='.format(len(ut)))
print(f'  min     = {ut.min():.4f}')
print(f'  p1      = {np.percentile(ut, 1):.4f}')
print(f'  median  = {np.median(ut):.4f}')
print(f'  p99     = {np.percentile(ut, 99):.4f}')
print(f'  max     = {ut.max():.4f}')
print()
print('Top 5 pairs by Jaccard (most suspicious):')
top = long_df.nlargest(5, 'jaccard')[['A','B','jaccard']]
print(top.to_string(index=False))
