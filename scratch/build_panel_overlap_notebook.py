import json
from pathlib import Path
ROOT = Path('/carnegie/nobackup/scratch/tbellagio/kmate')

def code(src):
    return {"cell_type": "code", "execution_count": None, "metadata": {},
            "outputs": [], "source": [l + '\n' for l in src.rstrip().split('\n')]}

def md(src):
    return {"cell_type": "markdown", "metadata": {},
            "source": [l + '\n' for l in src.rstrip().split('\n')]}

cells = []

cells.append(md("""# Panel-level carrier-set agreement: v3, v3 cactus-only, p82, v2 — all vs GrENE-Net

For each panel, compare per-SNP **founder carrier sets** (boolean 0/1 per founder per SNP)
against the GrENE-Net VCF carrier sets at matched (chrom, pos, REF, ALT) on Chr1.

This is **panel level** — the cn_var matrix that AF projection multiplies `h` against.
NOT the EM output. The disagreements measured here are present BEFORE any inference happens.

| panel | source | F overlap with GN |
|---|---|---|
| **v3 (all 231)** | cactus-pangenome (80) + PanGenie short-read (151) | 231 |
| **v3 cactus-only** | restrict v3 to the 80 cactus founders | 80 |
| **p82** | pure cactus pangenome of 82 long-read assemblies (no PanGenie) | 80 (overlapping with GN) |
| **v2** | Beagle-imputed from GN merged with cactus SVs | 231 |

**Two key ratio definitions:**

- **cells O/U** — across all (founder × SNP) cells, count of "panel says carrier, GN says non" vs "panel says non, GN says carrier". Per-cell asymmetry.
- **record O/U** — for each record, classify by NET direction (more cells over → record-over, more cells under → record-under). Per-record direction asymmetry.

These can diverge: a panel can have balanced CELLS (cell O/U ≈ 1) but skewed RECORDS (e.g., over-call broadly across many records, under-call deeply at fewer records).
"""))

cells.append(code("""import json
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

ROOT = Path('/carnegie/nobackup/scratch/tbellagio/kmate')

# Load 4 panel comparisons
results = []
for t in [1, 2, 3, 4]:
    f = ROOT / f'scratch/panel_cmp_task{t}.json'
    if f.exists():
        results.append(json.loads(f.read_text()))

df = pd.DataFrame(results)
df['exact_pct']      = 100 * df['exact'] / df['n_matched']
df['cell_O_U_ratio'] = df['cells_over'] / df['cells_under'].clip(lower=1)
df['rec_O_U_ratio']  = df['rec_over']  / df['rec_under'].clip(lower=1)
df['mean_over_cells_per_over_record']  = df['cells_over']  / df['rec_over'].clip(lower=1)
df['mean_under_cells_per_under_record'] = df['cells_under'] / df['rec_under'].clip(lower=1)
print(df[['label','n_matched','F','exact_pct','cell_O_U_ratio','rec_O_U_ratio',
          'mean_over_cells_per_over_record','mean_under_cells_per_under_record']].to_string(index=False, float_format='%.2f'))"""))

cells.append(md("""## Plot 1 — Composition of matched SNP records per panel

Stacked bars showing: exact match / net over-call / net under-call / net-zero (disagree but count same).
"""))

cells.append(code("""fig, ax = plt.subplots(figsize=(11, 6))
labels = df['label'].tolist()
x = np.arange(len(labels))
exact   = df['exact'].values
over    = df['rec_over'].values
under   = df['rec_under'].values
netzero = df['rec_net_zero'].values

bottom = np.zeros(len(labels))
for height, color, lab in [
    (exact,   '#7fbf7b', 'EXACT match (all founders agree)'),
    (over,    '#fc8d59', 'NET v3-over (panel has MORE carriers)'),
    (under,   '#91bfdb', 'NET v3-under (panel has FEWER carriers)'),
    (netzero, '#bdbdbd', 'NET zero (same count, diff founders)'),
]:
    ax.bar(x, height, bottom=bottom, color=color, label=lab, edgecolor='black', linewidth=0.5)
    bottom += height

ax.set_xticks(x); ax.set_xticklabels(labels, rotation=15, ha='right')
ax.set_ylabel('# matched SNP records')
ax.set_title('Per-record disagreement composition (panel vs GN, Chr1 SNPs)')
ax.legend(loc='upper right', fontsize=9)

# Annotate exact-match %
for xi, (e, total, lab) in enumerate(zip(exact, df['n_matched'].values, labels)):
    ax.text(xi, total + total*0.01, f'n={total:,}\\nexact={100*e/total:.1f}%',
            ha='center', va='bottom', fontsize=9)
ax.set_ylim(0, max(df['n_matched'].values) * 1.13)
plt.tight_layout()
plt.show()"""))

cells.append(md("""## Plot 2 — Same composition, normalized to %

Easier to compare across panels with different total record counts.
"""))

cells.append(code("""fig, ax = plt.subplots(figsize=(11, 6))
exact_pct   = 100 * exact / df['n_matched'].values
over_pct    = 100 * over / df['n_matched'].values
under_pct   = 100 * under / df['n_matched'].values
netzero_pct = 100 * netzero / df['n_matched'].values

bottom = np.zeros(len(labels))
for h, color, lab in [
    (exact_pct,   '#7fbf7b', 'EXACT match'),
    (over_pct,    '#fc8d59', 'NET over'),
    (under_pct,   '#91bfdb', 'NET under'),
    (netzero_pct, '#bdbdbd', 'NET zero'),
]:
    ax.bar(x, h, bottom=bottom, color=color, label=lab, edgecolor='black', linewidth=0.5)
    bottom += h
ax.set_xticks(x); ax.set_xticklabels(labels, rotation=15, ha='right')
ax.set_ylabel('% of matched SNP records')
ax.set_ylim(0, 105)
ax.set_title('Per-record composition (normalized %)')
for xi in range(len(labels)):
    ax.text(xi, 1, f'{exact_pct[xi]:.1f}%', ha='center', va='bottom', fontsize=9, color='black', fontweight='bold')
ax.legend(loc='upper right', fontsize=9)
plt.tight_layout()
plt.show()"""))

cells.append(md("""## Plot 3 — Cell-level totals (over vs under cells across all records)
"""))

cells.append(code("""fig, axes = plt.subplots(1, 2, figsize=(14, 5))

# Left: stacked cell counts
ax = axes[0]
co = df['cells_over'].values
cu = df['cells_under'].values
ax.bar(x - 0.2, co, width=0.4, color='#fc8d59', label='cells panel-over (panel=1, GN=0)', edgecolor='black')
ax.bar(x + 0.2, cu, width=0.4, color='#91bfdb', label='cells panel-under (panel=0, GN=1)', edgecolor='black')
ax.set_xticks(x); ax.set_xticklabels(labels, rotation=15, ha='right')
ax.set_ylabel('# cells (founder × SNP pairs)')
ax.set_title('Cell-level over vs under counts')
for xi in range(len(labels)):
    ratio = co[xi] / max(cu[xi], 1)
    ax.text(xi, max(co[xi], cu[xi]) + max(co.max(), cu.max())*0.02,
            f'O/U = {ratio:.2f}', ha='center', va='bottom', fontsize=10, fontweight='bold')
ax.legend(loc='upper right')

# Right: ratios
ax = axes[1]
ax.bar(x - 0.2, df['cell_O_U_ratio'], width=0.4, color='#9e3d22', label='cell O/U ratio', edgecolor='black')
ax.bar(x + 0.2, df['rec_O_U_ratio'], width=0.4, color='#3a567c', label='record O/U ratio', edgecolor='black')
ax.axhline(1.0, color='black', linestyle='--', lw=0.8, label='unbiased (=1)')
ax.set_xticks(x); ax.set_xticklabels(labels, rotation=15, ha='right')
ax.set_ylabel('over / under ratio')
ax.set_title('Asymmetry ratio: cell-level vs record-level')
for xi, v in enumerate(df['cell_O_U_ratio']):
    ax.text(xi - 0.2, v + 0.03, f'{v:.2f}', ha='center', fontsize=9)
for xi, v in enumerate(df['rec_O_U_ratio']):
    ax.text(xi + 0.2, v + 0.03, f'{v:.2f}', ha='center', fontsize=9)
ax.legend()
plt.tight_layout()
plt.show()"""))

cells.append(md("""## Plot 4 — Mean cell-disagreement per record by direction

Reveals the "broad shallow over-call vs deep concentrated under-call" pattern in p82.
"""))

cells.append(code("""fig, ax = plt.subplots(figsize=(10, 5))
m_over  = df['mean_over_cells_per_over_record']
m_under = df['mean_under_cells_per_under_record']
ax.bar(x - 0.2, m_over, width=0.4, color='#fc8d59',
       label='avg cells over per net-over record', edgecolor='black')
ax.bar(x + 0.2, m_under, width=0.4, color='#91bfdb',
       label='avg cells under per net-under record', edgecolor='black')
for xi, v in enumerate(m_over):
    ax.text(xi - 0.2, v + 0.05, f'{v:.2f}', ha='center', fontsize=9)
for xi, v in enumerate(m_under):
    ax.text(xi + 0.2, v + 0.05, f'{v:.2f}', ha='center', fontsize=9)
ax.set_xticks(x); ax.set_xticklabels(labels, rotation=15, ha='right')
ax.set_ylabel('mean # cells per record')
ax.set_title('How "deep" are over-call and under-call records on average?')
ax.legend()
plt.tight_layout()
plt.show()"""))

cells.append(md("""## Plot 5 — v3 vs GN per-record AF distribution

Distribution of per-record AF diff (v3 carriers / 231 − GN carriers / 231) across 437k matched Chr1 SNPs. Already-precomputed from `scratch/v3_vs_gn_per_record_af.npz`.
"""))

cells.append(code("""npz_path = ROOT / 'scratch/v3_vs_gn_per_record_af.npz'
if npz_path.exists():
    z = np.load(npz_path)
    diff = z['diff']
    fig, axes = plt.subplots(1, 2, figsize=(15, 5.5))
    # Histogram (log y)
    axes[0].hist(diff, bins=200, color='steelblue', edgecolor='none')
    axes[0].axvline(0, color='red', ls='--', lw=1, alpha=0.6)
    axes[0].axvline(diff.mean(), color='green', ls=':', lw=1.5, label=f'mean={diff.mean():+.5f}')
    axes[0].set_yscale('log')
    axes[0].set_xlabel('diff = v3_AF - GN_AF (per record)')
    axes[0].set_ylabel('# SNPs (log)')
    axes[0].set_title(f'v3 cn_var vs GN VCF — per-record AF diff distribution\\n'
                       f'mean={diff.mean():+.5f}  median={np.median(diff):+.5f}\\n'
                       f'big over (>+0.1): {(diff>0.1).sum():,}    big under (<-0.1): {(diff<-0.1).sum():,}    ratio={(diff>0.1).sum()/max((diff<-0.1).sum(),1):.2f}')
    axes[0].legend()

    # Hexbin v3 vs GN (panel level, no EM)
    axes[1].hexbin(z['gn_af'], z['v3_af'], gridsize=80, cmap='viridis', norm=LogNorm(), mincnt=1)
    axes[1].plot([0,1],[0,1],'r--', lw=0.8)
    axes[1].set_xlim(0,1); axes[1].set_ylim(0,1); axes[1].set_aspect('equal')
    axes[1].set_xlabel('GN_AF (carriers/231)'); axes[1].set_ylabel('v3_AF (carriers/231)')
    axes[1].set_title('v3 cn_var carrier-AF vs GN VCF carrier-AF\\n(panel level, BEFORE EM — pure carriers/231)')
    plt.tight_layout()
    plt.show()
else:
    print('per-record AF diff npz not found; run scratch/test_disagreement_directionality.py first')"""))

cells.append(md("""## Bottom line — what each panel says

| panel | exact match % | cell O/U | record O/U | interpretation |
|---|---|---|---|---|
| v2 (Beagle from GN) | ~100% | n/a | n/a | Same panel as GN; sanity baseline |
| **p82 (cactus only, no PG)** | 51.7% | **0.95** | 1.43 | **Symmetric at cell level** — cactus alone calls carriers ~as often as GN. Long-read assembly is NOT the source of v3's over-call asymmetry. |
| v3 cactus-only (in v3 pipeline) | 54.4% | 1.31 | 1.53 | Mild over-call. Comes from v3 build steps (atomization, merge with PG). |
| **v3 (all 231)** | 33.2% | **1.85** | 1.94 | Strongest over-call. The 151 PG founders' addition pushes the ratio from 1.31 → 1.85. |

**The 1.85× cell-level over-call in production v3 is driven by PanGenie genotyping, not by cactus/long-read discovery.** Pure cactus (p82) is balanced at the cell level (0.95). Each step downstream (cactus-only-in-v3 → all-of-v3) progressively introduces over-call asymmetry.

**The record-level O/U is always > cell-level O/U.** This means a panel's over-calls are spread across many records (each adding a few cells), while under-calls are concentrated in fewer records (each losing many founders). p82 shows this most cleanly: cell ratio 0.95 (balanced) but record ratio 1.43 (more records lean over-call).
"""))

nb = {"cells": cells,
      "metadata": {"kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
                   "language_info": {"name": "python", "pygments_lexer": "ipython3"}},
      "nbformat": 4, "nbformat_minor": 5}

out = ROOT / 'panel_overlap_4way.ipynb'
with open(out, 'w') as f: json.dump(nb, f, indent=1)
print(f'wrote {out}')
