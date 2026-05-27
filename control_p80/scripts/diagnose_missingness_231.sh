#!/bin/bash
#SBATCH --job-name=miss_231
#SBATCH --partition=bse
#SBATCH --cpus-per-task=2
#SBATCH --mem=96G
#SBATCH --time=1:00:00
#SBATCH --output=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/control_p80/logs/miss_231_%j.out
#SBATCH --error=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/control_p80/logs/miss_231_%j.err

# Compute missingness distribution on the production 231-panel (v3qc_v3)
# stratified by variant class + cactus/PG side. Output a PNG into
# control_p80/results/ (next to the p80 version for side-by-side reference).

set -euo pipefail

/home/tbellagio/miniforge3/envs/hapfm/bin/python <<'PYEOF'
import numpy as np
import pandas as pd
import scipy.sparse as sp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

DATA = Path('/carnegie/nobackup/scratch/tbellagio/hapfire_sv/poolfreq/data')
OUT  = Path('/carnegie/nobackup/scratch/tbellagio/hapfire_sv/control_p80/results')

# Load only the meta arrays we need (skip the heavy 'ref'/'alt' object strings).
print('Loading meta header info ...', flush=True)
with np.load(DATA / 'cn_var_231_v3qc_v3.meta.npz', allow_pickle=True) as m:
    founders = np.asarray(m['founders']).astype(str)
    ref_len  = m['ref_len'].astype(int)
    alt_len  = m['alt_len'].astype(int)
N = len(ref_len)
print(f'  meta: {N:,} records  {len(founders)} founders', flush=True)

print('Loading cn_var_called ...', flush=True)
called = sp.load_npz(DATA / 'cn_var_231_v3qc_v3.cn_var_called.npz')
F = called.shape[0]
assert called.shape[1] == N
print(f'  shape={called.shape}  nnz={called.nnz:,}  density={called.nnz/(F*N)*100:.2f}%', flush=True)

# Identify cactus side via the panel TSV
panel_tsv    = '/carnegie/nobackup/scratch/tbellagio/hapfire_sv/data/sv_panel_to_accession_id.tsv'
exclude_list = '/carnegie/nobackup/scratch/tbellagio/hapfire_sv/data/exclude_list.txt'
with open(exclude_list) as f:
    excluded_asm = {line.strip() for line in f if line.strip()}
sv_panel = pd.read_csv(panel_tsv, sep='\t')
cactus_acc = set(
    sv_panel.loc[~sv_panel['Assembly_ID'].astype(str).isin(excluded_asm),
                 'Accession_ID'].astype(str)
)
is_cactus = np.array([fid in cactus_acc for fid in founders])
F_c = int(is_cactus.sum()); F_p = F - F_c
print(f'  cactus side: {F_c} founders   PG side: {F_p} founders', flush=True)

# Per-record missing counts (overall + per side)
n_called   = np.asarray(called.sum(axis=0)).ravel().astype(int)
n_missing  = F - n_called
n_miss_c   = F_c - np.asarray(called[is_cactus,  :].sum(axis=0)).ravel().astype(int)
n_miss_p   = F_p - np.asarray(called[~is_cactus, :].sum(axis=0)).ravel().astype(int)

# Variant classes
max_len  = np.maximum(ref_len, alt_len)
is_snp   = (ref_len == 1) & (alt_len == 1)
is_sv    = max_len >= 50
is_indel = (~is_snp) & (~is_sv)

# ---------- TEXT REPORT ----------
print()
print(f'{"class":10s} {"n_records":>14s} {"% panel":>8s} {"mean #miss":>12s} {"median":>8s} {"max":>5s} {"% fully called":>16s} {"% w/>50 miss":>14s}')
for label, mask in [('SNP', is_snp), ('indel', is_indel), ('SV', is_sv), ('ALL', np.ones(N, dtype=bool))]:
    n = int(mask.sum())
    if n == 0: continue
    nm = n_missing[mask]
    print(f'{label:10s} {n:>14,} {n/N*100:>7.2f}% '
          f'{nm.mean():>12.2f} {int(np.median(nm)):>8d} {nm.max():>5d} '
          f'{(nm==0).mean()*100:>15.2f}% {(nm>50).mean()*100:>13.2f}%')

print()
print(f'=== Per-side mean missingness across records ===')
print(f'  cactus ({F_c} founders): mean #miss = {n_miss_c.mean():.2f}  '
      f'fully-called rate = {(n_miss_c==0).mean()*100:.2f}%')
print(f'  PG     ({F_p} founders): mean #miss = {n_miss_p.mean():.2f}  '
      f'fully-called rate = {(n_miss_p==0).mean()*100:.2f}%')

print()
print(f'=== Per-side x variant class (mean #missing per side) ===')
print(f'{"class":10s} {"cactus #miss":>15s} {"%cactus called":>16s} {"PG #miss":>12s} {"%PG called":>12s}')
for label, mask in [('SNP', is_snp), ('indel', is_indel), ('SV', is_sv)]:
    c_m = n_miss_c[mask].mean()
    p_m = n_miss_p[mask].mean()
    c_pct = (1 - c_m/F_c) * 100
    p_pct = (1 - p_m/F_p) * 100
    print(f'{label:10s} {c_m:>15.2f} {c_pct:>15.2f}% {p_m:>12.2f} {p_pct:>11.2f}%')

print()
print(f'=== % records retained per class by max-#missing threshold ===')
print(f'{"max_miss":>10s}  {"SNP %kept":>10s}  {"indel %kept":>12s}  {"SV %kept":>10s}  {"ALL %kept":>10s}')
for thr in [0, 1, 5, 10, 25, 50, 100, 150, 200, F-1]:
    snp_k = (n_missing[is_snp]   <= thr).mean() * 100 if is_snp.any()   else 0
    ind_k = (n_missing[is_indel] <= thr).mean() * 100 if is_indel.any() else 0
    sv_k  = (n_missing[is_sv]    <= thr).mean() * 100 if is_sv.any()    else 0
    all_k = (n_missing             <= thr).mean() * 100
    print(f'{thr:>10d}  {snp_k:>9.2f}%  {ind_k:>11.2f}%  {sv_k:>9.2f}%  {all_k:>9.2f}%')

# ---------- PLOTS ----------
colors = {'SNP': 'steelblue', 'indel': 'darkorange', 'SV': 'crimson'}
fig, axes = plt.subplots(2, 2, figsize=(15, 10))

# (1,1) overall counts by class
ax = axes[0, 0]
for label, mask in [('SNP', is_snp), ('indel', is_indel), ('SV', is_sv)]:
    counts = np.bincount(n_missing[mask], minlength=F+1)[:F+1]
    ax.bar(np.arange(F+1), counts, width=1.0, color=colors[label],
           label=f'{label} (n={int(mask.sum()):,})', alpha=0.55, edgecolor='none')
ax.set_yscale('log')
ax.set_xlabel('# missing founders per record (out of 231)')
ax.set_ylabel('# records (log)')
ax.set_title('231-panel: missingness distribution by variant class (counts)')
ax.legend()
ax.set_xlim(-0.5, F+0.5)

# (1,2) within-class proportion
ax = axes[0, 1]
for label, mask in [('SNP', is_snp), ('indel', is_indel), ('SV', is_sv)]:
    counts = np.bincount(n_missing[mask], minlength=F+1)[:F+1]
    frac = counts / max(mask.sum(), 1) * 100
    ax.bar(np.arange(F+1), frac, width=1.0, color=colors[label],
           label=label, alpha=0.55, edgecolor='none')
ax.set_yscale('log')
ax.set_xlabel('# missing founders per record')
ax.set_ylabel('% of class records')
ax.set_title('Within-class % at each missingness level')
ax.legend()
ax.set_xlim(-0.5, F+0.5)

# (2,1) cactus-only missingness by var class
ax = axes[1, 0]
for label, mask in [('SNP', is_snp), ('indel', is_indel), ('SV', is_sv)]:
    counts = np.bincount(n_miss_c[mask], minlength=F_c+1)[:F_c+1]
    frac = counts / max(mask.sum(), 1) * 100
    ax.bar(np.arange(F_c+1), frac, width=1.0, color=colors[label],
           label=label, alpha=0.55, edgecolor='none')
ax.set_yscale('log')
ax.set_xlabel(f'# missing CACTUS founders (out of {F_c})')
ax.set_ylabel('% of class records')
ax.set_title('Cactus-only missingness by variant class')
ax.legend()
ax.set_xlim(-0.5, F_c+0.5)

# (2,2) PG-only missingness by var class
ax = axes[1, 1]
for label, mask in [('SNP', is_snp), ('indel', is_indel), ('SV', is_sv)]:
    counts = np.bincount(n_miss_p[mask], minlength=F_p+1)[:F_p+1]
    frac = counts / max(mask.sum(), 1) * 100
    ax.bar(np.arange(F_p+1), frac, width=1.0, color=colors[label],
           label=label, alpha=0.55, edgecolor='none')
ax.set_yscale('log')
ax.set_xlabel(f'# missing PG founders (out of {F_p})')
ax.set_ylabel('% of class records')
ax.set_title('PG-only missingness by variant class')
ax.legend()
ax.set_xlim(-0.5, F_p+0.5)

plt.suptitle('231-panel (v3qc_v3) missingness — overall + per-side decomposition', y=1.005, fontsize=13)
plt.tight_layout()
plt.savefig(OUT / 'missingness_by_var_class_231panel.png', dpi=130, bbox_inches='tight')
print(f'\nsaved: {OUT / "missingness_by_var_class_231panel.png"}', flush=True)
PYEOF
