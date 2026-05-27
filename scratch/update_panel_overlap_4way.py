"""Surgically update panel_overlap_4way.ipynb:
- Drop v2 from data load + all plots
- Annotate plot 1 (stacked counts) with EXACT numbers per segment
- Add a cell-level match-rate plot (F-comparable)
- Replace plot 5 with 3-panel hexbin (v3, v3 cactus-only, p82)
"""
import json, copy
from pathlib import Path

NB = Path('/carnegie/nobackup/scratch/tbellagio/hapfire_sv/panel_overlap_4way.ipynb')
nb = json.loads(NB.read_text())

def src(s): return s if isinstance(s, list) else s.splitlines(keepends=True)

# Find cells by id (the harness shows ids)
def find_cell(id_prefix):
    for i, c in enumerate(nb['cells']):
        if c.get('id', '').startswith(id_prefix):
            return i
    raise KeyError(id_prefix)

# Cell 0: header markdown — update table to drop v2
i = 0
nb['cells'][i]['source'] = [
    "# Panel-level carrier-set agreement: v3, v3 cactus-only, p82 — all vs GrENE-Net\n",
    "\n",
    "For each panel, compare per-SNP **founder carrier sets** (boolean 0/1 per founder per SNP)\n",
    "against the GrENE-Net VCF carrier sets at matched (chrom, pos, REF, ALT) on Chr1.\n",
    "\n",
    "This is **panel level** — the cn_var matrix that AF projection multiplies `h` against.\n",
    "NOT the EM output. The disagreements measured here are present BEFORE any inference happens.\n",
    "\n",
    "| panel | source | F overlap with GN |\n",
    "|---|---|---|\n",
    "| **v3 (all 231)** | cactus-pangenome (80) + PanGenie short-read (151) | 231 |\n",
    "| **v3 cactus-only** | restrict v3 to the 80 cactus founders | 80 |\n",
    "| **p82** | pure cactus pangenome of 82 long-read assemblies (no PanGenie) | 80 (overlapping with GN) |\n",
    "\n",
    "**Caveat on exact-match %**: comparing exact-match rates across panels with different F is misleading. A panel with F=231 must agree on 231 cells per record to be \"exact\"; F=80 must agree on 80 cells. So 33% (v3, F=231) vs 54% (v3 cactus-only, F=80) is NOT apples-to-apples. The **per-cell match rate** below IS F-comparable.\n",
    "\n",
    "**Two key ratio definitions:**\n",
    "\n",
    "- **cells O/U** — across all (founder × SNP) cells, count of \"panel says carrier, GN says non\" vs \"panel says non, GN says carrier\". Per-cell asymmetry.\n",
    "- **record O/U** — for each record, classify by NET direction (more cells over → record-over, more cells under → record-under). Per-record direction asymmetry.\n",
    "\n",
    "These can diverge: a panel can have balanced CELLS (cell O/U ≈ 1) but skewed RECORDS (e.g., over-call broadly across many records, under-call deeply at fewer records).\n"
]

# Cell 1: load + table — drop task 4 (v2), add cell-level match-rate columns
i = find_cell('2c5f440e')
nb['cells'][i]['source'] = [
    "import json\n",
    "from pathlib import Path\n",
    "import numpy as np\n",
    "import pandas as pd\n",
    "import matplotlib.pyplot as plt\n",
    "from matplotlib.colors import LogNorm\n",
    "\n",
    "ROOT = Path('/carnegie/nobackup/scratch/tbellagio/hapfire_sv')\n",
    "\n",
    "# Load 3 panel comparisons (drop task 4 = v2 — same as GN, not informative)\n",
    "results = []\n",
    "for t in [1, 2, 3]:\n",
    "    f = ROOT / f'scratch/panel_cmp_task{t}.json'\n",
    "    if f.exists():\n",
    "        results.append(json.loads(f.read_text()))\n",
    "\n",
    "df = pd.DataFrame(results)\n",
    "df['exact_pct']      = 100 * df['exact'] / df['n_matched']\n",
    "df['cell_O_U_ratio'] = df['cells_over'] / df['cells_under'].clip(lower=1)\n",
    "df['rec_O_U_ratio']  = df['rec_over']  / df['rec_under'].clip(lower=1)\n",
    "df['mean_over_cells_per_over_record']  = df['cells_over']  / df['rec_over'].clip(lower=1)\n",
    "df['mean_under_cells_per_under_record'] = df['cells_under'] / df['rec_under'].clip(lower=1)\n",
    "# F-comparable per-cell metric\n",
    "df['total_cells']        = df['F'] * df['n_matched']\n",
    "df['disagree_cells']     = df['cells_over'] + df['cells_under']\n",
    "df['cell_match_pct']     = 100 * (1 - df['disagree_cells'] / df['total_cells'])\n",
    "df['cell_disagree_pct']  = 100 * df['disagree_cells'] / df['total_cells']\n",
    "print(df[['label','n_matched','F','total_cells','exact_pct','cell_match_pct',\n",
    "          'cell_O_U_ratio','rec_O_U_ratio',\n",
    "          'mean_over_cells_per_over_record','mean_under_cells_per_under_record']].to_string(index=False, float_format='%.3f'))\n",
]
nb['cells'][i]['outputs'] = []
nb['cells'][i]['execution_count'] = None

# Plot 1 markdown
i_md_plot1 = find_cell('b2d546f7')
nb['cells'][i_md_plot1]['source'] = [
    "## Plot 1 — Composition of matched SNP records per panel (with exact counts)\n",
    "\n",
    "Stacked bars showing: exact match / net over-call / net under-call / net-zero. Each segment annotated with EXACT record count (not just %).\n"
]

# Plot 1 code: annotate each segment with exact count
i = find_cell('57feae5c')
nb['cells'][i]['source'] = [
    "fig, ax = plt.subplots(figsize=(11, 7))\n",
    "labels = df['label'].tolist()\n",
    "x = np.arange(len(labels))\n",
    "exact   = df['exact'].values\n",
    "over    = df['rec_over'].values\n",
    "under   = df['rec_under'].values\n",
    "netzero = df['rec_net_zero'].values\n",
    "\n",
    "bottom = np.zeros(len(labels))\n",
    "for height, color, lab in [\n",
    "    (exact,   '#7fbf7b', 'EXACT match (all founders agree)'),\n",
    "    (over,    '#fc8d59', 'NET panel-over (panel has MORE carriers)'),\n",
    "    (under,   '#91bfdb', 'NET panel-under (panel has FEWER carriers)'),\n",
    "    (netzero, '#bdbdbd', 'NET zero (same count, diff founders)'),\n",
    "]:\n",
    "    bars = ax.bar(x, height, bottom=bottom, color=color, label=lab, edgecolor='black', linewidth=0.5)\n",
    "    # Annotate each segment with its exact count, centered in the segment\n",
    "    for xi, (h, b) in enumerate(zip(height, bottom)):\n",
    "        if h > 0:\n",
    "            ax.text(xi, b + h/2, f'{int(h):,}', ha='center', va='center',\n",
    "                    fontsize=9, color='black', fontweight='bold')\n",
    "    bottom += height\n",
    "\n",
    "ax.set_xticks(x); ax.set_xticklabels(labels, rotation=15, ha='right')\n",
    "ax.set_ylabel('# matched SNP records')\n",
    "ax.set_title('Per-record disagreement composition (panel vs GN, Chr1 SNPs) — exact counts')\n",
    "ax.legend(loc='upper right', fontsize=9)\n",
    "\n",
    "for xi, (e, total) in enumerate(zip(exact, df['n_matched'].values)):\n",
    "    ax.text(xi, total + total*0.015, f'n={total:,}\\nexact={100*e/total:.1f}%',\n",
    "            ha='center', va='bottom', fontsize=9)\n",
    "ax.set_ylim(0, max(df['n_matched'].values) * 1.18)\n",
    "plt.tight_layout()\n",
    "plt.show()\n"
]
nb['cells'][i]['outputs'] = []
nb['cells'][i]['execution_count'] = None

# Plot 2 markdown — keep but acknowledge F-caveat
i_md = find_cell('d3cfeb63')
nb['cells'][i_md]['source'] = [
    "## Plot 2 — Same composition, normalized to %\n",
    "\n",
    "Easier to compare across panels with different total record counts. **Note**: exact-match % is per-record (depends on F); see Plot 3 right panel for F-comparable per-cell match rate.\n"
]

# Plot 2 code: same structure, just operating on the 3-panel df. Annotate counts AND % per segment.
i = find_cell('a96b515a')
nb['cells'][i]['source'] = [
    "fig, ax = plt.subplots(figsize=(11, 7))\n",
    "exact_pct   = 100 * exact / df['n_matched'].values\n",
    "over_pct    = 100 * over / df['n_matched'].values\n",
    "under_pct   = 100 * under / df['n_matched'].values\n",
    "netzero_pct = 100 * netzero / df['n_matched'].values\n",
    "\n",
    "bottom = np.zeros(len(labels))\n",
    "for h_pct, h_count, color, lab in [\n",
    "    (exact_pct,   exact,   '#7fbf7b', 'EXACT match'),\n",
    "    (over_pct,    over,    '#fc8d59', 'NET over'),\n",
    "    (under_pct,   under,   '#91bfdb', 'NET under'),\n",
    "    (netzero_pct, netzero, '#bdbdbd', 'NET zero'),\n",
    "]:\n",
    "    ax.bar(x, h_pct, bottom=bottom, color=color, label=lab, edgecolor='black', linewidth=0.5)\n",
    "    for xi, (hp, hc, b) in enumerate(zip(h_pct, h_count, bottom)):\n",
    "        if hp > 1.5:\n",
    "            ax.text(xi, b + hp/2, f'{hp:.1f}%\\n({int(hc):,})',\n",
    "                    ha='center', va='center', fontsize=8.5, color='black', fontweight='bold')\n",
    "    bottom += h_pct\n",
    "ax.set_xticks(x); ax.set_xticklabels(labels, rotation=15, ha='right')\n",
    "ax.set_ylabel('% of matched SNP records')\n",
    "ax.set_ylim(0, 105)\n",
    "ax.set_title('Per-record composition (normalized %) — count in parens')\n",
    "ax.legend(loc='upper right', fontsize=9)\n",
    "plt.tight_layout()\n",
    "plt.show()\n"
]
nb['cells'][i]['outputs'] = []
nb['cells'][i]['execution_count'] = None

# Plot 3 markdown
i_md = find_cell('db5f8368')
nb['cells'][i_md]['source'] = [
    "## Plot 3 — Cell-level totals + F-comparable per-cell match rate\n",
    "\n",
    "**Right panel is the F-comparable metric**: per-cell match rate divides by total cells (F × n_matched), so a panel with more founders isn't penalized for having more cells to disagree on.\n"
]

# Plot 3 code: cell counts + ratio + ALSO add F-comparable per-cell match rate
i = find_cell('31ffc472')
nb['cells'][i]['source'] = [
    "fig, axes = plt.subplots(1, 3, figsize=(20, 5.5))\n",
    "\n",
    "# Left: stacked cell counts\n",
    "ax = axes[0]\n",
    "co = df['cells_over'].values\n",
    "cu = df['cells_under'].values\n",
    "ax.bar(x - 0.2, co, width=0.4, color='#fc8d59', label='cells panel-over (panel=1, GN=0)', edgecolor='black')\n",
    "ax.bar(x + 0.2, cu, width=0.4, color='#91bfdb', label='cells panel-under (panel=0, GN=1)', edgecolor='black')\n",
    "ax.set_xticks(x); ax.set_xticklabels(labels, rotation=15, ha='right')\n",
    "ax.set_ylabel('# cells (founder × SNP pairs)')\n",
    "ax.set_title('Cell-level over vs under counts')\n",
    "for xi in range(len(labels)):\n",
    "    ax.text(xi - 0.2, co[xi] + max(co.max(), cu.max())*0.01, f'{int(co[xi]):,}', ha='center', fontsize=8.5)\n",
    "    ax.text(xi + 0.2, cu[xi] + max(co.max(), cu.max())*0.01, f'{int(cu[xi]):,}', ha='center', fontsize=8.5)\n",
    "    ratio = co[xi] / max(cu[xi], 1)\n",
    "    ax.text(xi, max(co[xi], cu[xi]) + max(co.max(), cu.max())*0.06,\n",
    "            f'O/U = {ratio:.2f}', ha='center', va='bottom', fontsize=10, fontweight='bold')\n",
    "ax.set_ylim(0, max(co.max(), cu.max()) * 1.18)\n",
    "ax.legend(loc='upper right', fontsize=9)\n",
    "\n",
    "# Middle: ratios\n",
    "ax = axes[1]\n",
    "ax.bar(x - 0.2, df['cell_O_U_ratio'], width=0.4, color='#9e3d22', label='cell O/U ratio', edgecolor='black')\n",
    "ax.bar(x + 0.2, df['rec_O_U_ratio'], width=0.4, color='#3a567c', label='record O/U ratio', edgecolor='black')\n",
    "ax.axhline(1.0, color='black', linestyle='--', lw=0.8, label='unbiased (=1)')\n",
    "ax.set_xticks(x); ax.set_xticklabels(labels, rotation=15, ha='right')\n",
    "ax.set_ylabel('over / under ratio')\n",
    "ax.set_title('Asymmetry ratio: cell-level vs record-level')\n",
    "for xi, v in enumerate(df['cell_O_U_ratio']):\n",
    "    ax.text(xi - 0.2, v + 0.04, f'{v:.2f}', ha='center', fontsize=9)\n",
    "for xi, v in enumerate(df['rec_O_U_ratio']):\n",
    "    ax.text(xi + 0.2, v + 0.04, f'{v:.2f}', ha='center', fontsize=9)\n",
    "ax.legend()\n",
    "\n",
    "# Right: F-COMPARABLE per-cell match rate\n",
    "ax = axes[2]\n",
    "match_pct = df['cell_match_pct'].values\n",
    "disagree_pct = df['cell_disagree_pct'].values\n",
    "bars = ax.bar(x, match_pct, color='#7fbf7b', edgecolor='black', label='cells matching')\n",
    "ax.set_ylim(min(match_pct.min() - 0.5, 98), 100.0)\n",
    "ax.set_xticks(x); ax.set_xticklabels(labels, rotation=15, ha='right')\n",
    "ax.set_ylabel('% cells matching panel↔GN')\n",
    "ax.set_title('F-COMPARABLE: per-cell match rate\\n(disagree_cells / (F × n_matched))')\n",
    "for xi in range(len(labels)):\n",
    "    ax.text(xi, match_pct[xi] - 0.05, f'{match_pct[xi]:.3f}%\\n({int(df[\"disagree_cells\"].values[xi]):,} disagree\\nof {int(df[\"total_cells\"].values[xi]):,} cells)',\n",
    "            ha='center', va='top', fontsize=9, fontweight='bold')\n",
    "plt.tight_layout()\n",
    "plt.show()\n"
]
nb['cells'][i]['outputs'] = []
nb['cells'][i]['execution_count'] = None

# Plot 4 markdown — unchanged but tweak references
i_md = find_cell('9800524d')
nb['cells'][i_md]['source'] = [
    "## Plot 4 — Mean cell-disagreement per record by direction\n",
    "\n",
    "Reveals the \"broad shallow over-call vs deep concentrated under-call\" pattern in p82.\n"
]

# Plot 4 code: same structure, 3 panels
i = find_cell('15c34be7')
nb['cells'][i]['source'] = [
    "fig, ax = plt.subplots(figsize=(10, 5))\n",
    "m_over  = df['mean_over_cells_per_over_record']\n",
    "m_under = df['mean_under_cells_per_under_record']\n",
    "ax.bar(x - 0.2, m_over, width=0.4, color='#fc8d59',\n",
    "       label='avg cells over per net-over record', edgecolor='black')\n",
    "ax.bar(x + 0.2, m_under, width=0.4, color='#91bfdb',\n",
    "       label='avg cells under per net-under record', edgecolor='black')\n",
    "for xi, v in enumerate(m_over):\n",
    "    ax.text(xi - 0.2, v + 0.05, f'{v:.2f}', ha='center', fontsize=9)\n",
    "for xi, v in enumerate(m_under):\n",
    "    ax.text(xi + 0.2, v + 0.05, f'{v:.2f}', ha='center', fontsize=9)\n",
    "ax.set_xticks(x); ax.set_xticklabels(labels, rotation=15, ha='right')\n",
    "ax.set_ylabel('mean # cells per record')\n",
    "ax.set_title('How \"deep\" are over-call and under-call records on average?')\n",
    "ax.legend()\n",
    "plt.tight_layout()\n",
    "plt.show()\n"
]
nb['cells'][i]['outputs'] = []
nb['cells'][i]['execution_count'] = None

# Plot 5 markdown — replace
i_md = find_cell('bd018fa8')
nb['cells'][i_md]['source'] = [
    "## Plot 5 — Per-record carrier-AF: panel vs GN (3 panels, hexbin)\n",
    "\n",
    "For each panel, plot per-record carrier-AF (carriers / F_overlap) against GN's carrier-AF on the matched (chrom, pos, REF, ALT) intersection. Same x/y scale; red dashed = identity. Off-diagonal density = atomization-induced GT disagreement BEFORE EM.\n",
    "\n",
    "Inputs: precomputed at `scratch/per_record_af_{v3,v3_cactus_only,p82}.npz` by `scratch/compute_per_record_af_panels.py` (SLURM array 1-3).\n"
]

# Plot 5 code: 3-panel hexbin
i = find_cell('7d94d7a4')
nb['cells'][i]['source'] = [
    "panel_specs = [\n",
    "    ('v3 (all 231)',         ROOT/'scratch/per_record_af_v3.npz'),\n",
    "    ('v3 cactus-only (80)',  ROOT/'scratch/per_record_af_v3_cactus_only.npz'),\n",
    "    ('p82 (82 cactus, no PG)', ROOT/'scratch/per_record_af_p82.npz'),\n",
    "]\n",
    "missing = [p for _, p in panel_specs if not p.exists()]\n",
    "if missing:\n",
    "    print('Missing per-record AF npz files:')\n",
    "    for p in missing: print(' ', p)\n",
    "    print('Run: sbatch scratch/compute_per_record_af_panels.sbatch  (then re-run this cell)')\n",
    "else:\n",
    "    fig, axes = plt.subplots(1, 3, figsize=(20, 6))\n",
    "    for ax, (lab, p) in zip(axes, panel_specs):\n",
    "        z = np.load(p, allow_pickle=True)\n",
    "        gn_af = z['gn_af']; panel_af = z['panel_af']; diff = z['diff']\n",
    "        F = int(z['F_overlap']) if 'F_overlap' in z.files else None\n",
    "        hb = ax.hexbin(gn_af, panel_af, gridsize=80, cmap='viridis', norm=LogNorm(), mincnt=1)\n",
    "        ax.plot([0,1],[0,1],'r--', lw=0.8)\n",
    "        ax.set_xlim(0,1); ax.set_ylim(0,1); ax.set_aspect('equal')\n",
    "        ax.set_xlabel(f'GN_AF (carriers/{F})')\n",
    "        ax.set_ylabel(f'panel_AF (carriers/{F})')\n",
    "        ax.set_title(f'{lab}\\nn_records={len(diff):,}  '\n",
    "                     f'mean_diff={diff.mean():+.4f}  '\n",
    "                     f'|diff|>0.1: {(np.abs(diff)>0.1).sum():,}  '\n",
    "                     f'(>+0.1 / <-0.1 = {(diff>0.1).sum():,} / {(diff<-0.1).sum():,})')\n",
    "        plt.colorbar(hb, ax=ax, fraction=0.046, pad=0.04, label='# SNPs (log)')\n",
    "    plt.suptitle('Per-record carrier-AF: panel vs GrENE-Net (panel level, BEFORE EM)', y=1.02, fontsize=12)\n",
    "    plt.tight_layout()\n",
    "    plt.show()\n"
]
nb['cells'][i]['outputs'] = []
nb['cells'][i]['execution_count'] = None

# Bottom-line markdown — drop v2 row, mention cell match rate
i_md = find_cell('148e6cea')
nb['cells'][i_md]['source'] = [
    "## Bottom line — what each panel says\n",
    "\n",
    "| panel | F | exact match % (per-record) | per-cell match % (F-comparable) | cell O/U | record O/U | interpretation |\n",
    "|---|---|---|---|---|---|---|\n",
    "| **p82 (cactus only, no PG)** | 80 | 51.7% | (see Plot 3 right) | **0.95** | 1.43 | **Symmetric at cell level** — cactus alone calls carriers ~as often as GN. Long-read assembly is NOT the source of v3's over-call asymmetry. |\n",
    "| v3 cactus-only (in v3 pipeline) | 80 | 54.4% | (see Plot 3 right) | 1.31 | 1.53 | Mild over-call. Comes from v3 build steps (atomization, merge with PG). |\n",
    "| **v3 (all 231)** | 231 | 33.2% | (see Plot 3 right) | **1.85** | 1.94 | Strongest over-call. The 151 PG founders' addition pushes the ratio from 1.31 → 1.85. |\n",
    "\n",
    "**Caveat repeated**: the 33% vs 54% gap on \"exact match per-record\" is partly a denominator artifact — v3 has 231 cells per record vs 80 for the cactus-only restrictions. The **per-cell match %** in Plot 3 (right) is the apples-to-apples metric: it normalizes by F.\n",
    "\n",
    "**The 1.85× cell-level over-call in production v3 is driven by PanGenie genotyping, not by cactus/long-read discovery.** Pure cactus (p82) is balanced at the cell level (0.95). Each step downstream (cactus-only-in-v3 → all-of-v3) progressively introduces over-call asymmetry.\n",
    "\n",
    "**The record-level O/U is always > cell-level O/U.** Over-calls are spread across many records (each adding a few cells), while under-calls are concentrated in fewer records (each losing many founders). p82 shows this most cleanly: cell ratio 0.95 (balanced) but record ratio 1.43 (more records lean over-call).\n"
]

NB.write_text(json.dumps(nb, indent=1))
print('OK — wrote', NB)
