"""Plot h values: hapFIRE vs old v3 vs new v3qc cactus_em.
Comparator methods in neutral gray; v3qc colored by cactus/PG to show the bias."""
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

ROOT = Path(__file__).resolve().parent
HF_DIR = Path('/global/scratch/projects/fc_moilab/projects/grenenet-phase1/frequency/hapFIRE_frequencies/seed_mix')
cactus_set = set(json.load(open(ROOT / 'data/founder_split_cactus_pg.json'))['cactus']) - {'5772', '9947'}

def load_cem_h(d, s):
    z = np.load(d / f'SEEDMIX_S{s}.h_per_chrom.npz', allow_pickle=True)
    h_chroms = np.stack([z[c] for c in ['Chr1', 'Chr2', 'Chr3', 'Chr4', 'Chr5']])
    return h_chroms.mean(axis=0), list(np.asarray(z['founders']).astype(str))

def load_hf_h(s):
    df = pd.read_csv(HF_DIR / f's{s}_ecotype_frequency.txt', sep='\t',
                     header=None, names=['ecotype', 'freq'])
    return dict(zip(df['ecotype'].astype(str), df['freq'].astype(float)))

UNIFORM = 1.0 / 231.0
CACT_BLUE = '#1f77b4'
PG_ORANGE = '#ff7f0e'
GRAY_DARK = '#404040'
GRAY_LIGHT = '#909090'

# ===== 2x4 panel: one per replicate =====
fig, axes = plt.subplots(2, 4, figsize=(24, 11), sharey=True)
axes = axes.flatten()

for s in range(1, 9):
    ax = axes[s-1]
    v3_h, v3_f       = load_cem_h(ROOT / 'scratch/seedmix_v3_filt2_dedup', s)
    v3qc_h, v3qc_f   = load_cem_h(ROOT / 'scratch/seedmix_v3qc_filt2_dedup', s)
    hf = load_hf_h(s)
    hf_h = np.array([hf.get(f, np.nan) for f in v3qc_f])
    v3_idx = {f: i for i, f in enumerate(v3_f)}
    v3_h_aligned = np.array([v3_h[v3_idx[f]] if f in v3_idx else np.nan for f in v3qc_f])
    cact_mask = np.array([f in cactus_set for f in v3qc_f])

    order = np.argsort(v3qc_f)
    x = np.arange(len(v3qc_f))

    # Comparator methods in gray (no cactus/PG distinction)
    ax.scatter(x, hf_h[order],        c=GRAY_DARK,  marker='o', s=22, alpha=0.55,
               label='hapFIRE' if s==1 else None, edgecolor='none')
    ax.scatter(x, v3_h_aligned[order], c=GRAY_LIGHT, marker='^', s=22, alpha=0.55,
               label='cactus_em-v3' if s==1 else None, edgecolor='none')

    # v3qc colored by side
    cact_ord = cact_mask[order]
    v3qc_ord = v3qc_h[order]
    ax.scatter(x[cact_ord],  v3qc_ord[cact_ord],  c=CACT_BLUE, marker='x', s=30, lw=1.4,
               alpha=0.9, label='cactus_em-v3qc (cactus)' if s==1 else None)
    ax.scatter(x[~cact_ord], v3qc_ord[~cact_ord], c=PG_ORANGE, marker='x', s=30, lw=1.4,
               alpha=0.9, label='cactus_em-v3qc (PG)' if s==1 else None)

    ax.axhline(UNIFORM, color='red', linestyle='--', lw=1, alpha=0.8,
               label='1/231 uniform' if s==1 else None)
    ax.set_title(f'SEEDMIX_S{s}', fontsize=11)
    ax.set_xlabel('founder (sorted by ID)')
    if s == 1 or s == 5: ax.set_ylabel('h value')
    if s == 1: ax.legend(loc='upper left', fontsize=8)
    ax.set_ylim(-0.001, max(0.035, np.nanmax(v3qc_h) * 1.05))

fig.suptitle('h vectors per founder: hapFIRE & cactus_em-v3 (gray, comparators) vs cactus_em-v3qc (colored by side) — 8 SEEDMIX reps',
             fontsize=13, y=0.995)
plt.tight_layout()
plt.savefig('H_VECTOR_COMPARISON.png', dpi=130, bbox_inches='tight')
plt.close()

# ===== Single-panel S8 detail =====
fig, ax = plt.subplots(figsize=(15, 6))
s = 8
v3_h, v3_f     = load_cem_h(ROOT / 'scratch/seedmix_v3_filt2_dedup', s)
v3qc_h, v3qc_f = load_cem_h(ROOT / 'scratch/seedmix_v3qc_filt2_dedup', s)
hf = load_hf_h(s)
hf_h = np.array([hf.get(f, np.nan) for f in v3qc_f])
v3_idx = {f: i for i, f in enumerate(v3_f)}
v3_h_aligned = np.array([v3_h[v3_idx[f]] if f in v3_idx else np.nan for f in v3qc_f])
cact_mask = np.array([f in cactus_set for f in v3qc_f])
order = np.argsort(v3qc_f)
x = np.arange(len(v3qc_f))

ax.scatter(x, hf_h[order],        c=GRAY_DARK,  marker='o', s=35, alpha=0.6,
           label='hapFIRE (comparator)', edgecolor='none')
ax.scatter(x, v3_h_aligned[order], c=GRAY_LIGHT, marker='^', s=40, alpha=0.65,
           label='cactus_em-v3 (working)', edgecolor='none')

cact_ord = cact_mask[order]
v3qc_ord = v3qc_h[order]
ax.scatter(x[cact_ord],  v3qc_ord[cact_ord],  c=CACT_BLUE, marker='x', s=55, lw=1.8,
           alpha=0.95, label='cactus_em-v3qc — cactus side (78)')
ax.scatter(x[~cact_ord], v3qc_ord[~cact_ord], c=PG_ORANGE, marker='x', s=55, lw=1.8,
           alpha=0.95, label='cactus_em-v3qc — PG side (153)')

ax.axhline(UNIFORM, color='red', linestyle='--', lw=1.5, alpha=0.85,
           label=f'1/231 uniform expectation = {UNIFORM:.5f}')
ax.set_xlabel('founder (sorted by ID)')
ax.set_ylabel('h value (5-chrom mean)')
ax.set_title('SEEDMIX_S8 — h vector per founder per method\n'
             'hapFIRE & v3 cluster at uniform expectation; v3qc puts blue (cactus) ~3× over and orange (PG) ~10× under.')
ax.legend(loc='upper right', fontsize=10)
ax.set_ylim(-0.001, 0.032)
plt.tight_layout()
plt.savefig('H_VECTOR_COMPARISON_S8.png', dpi=130, bbox_inches='tight')
plt.close()
print('Saved H_VECTOR_COMPARISON.png and H_VECTOR_COMPARISON_S8.png')
