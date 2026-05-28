"""Add Plot 6 (fair comparison: v3 cactus-only vs p82 on the SHARED variant set,
same 80 founders, same GN truth) to panel_overlap_4way.ipynb.
"""
import json
from pathlib import Path

NB = Path('/carnegie/nobackup/scratch/tbellagio/kmate/panel_overlap_4way.ipynb')
nb = json.loads(NB.read_text())

# Find bottom-line markdown to insert before it
def find_cell(id_prefix):
    for i, c in enumerate(nb['cells']):
        if c.get('id', '').startswith(id_prefix): return i
    raise KeyError(id_prefix)

i_bottom = find_cell('148e6cea')

md_cell = {
    'cell_type': 'markdown',
    'id': 'fair6md001',
    'metadata': {},
    'source': [
        "## Plot 6 — FAIR comparison: v3 cactus-only vs p82 on the shared variant set\n",
        "\n",
        "Plot 5 made v3-cactus-only and p82 look very different. That's mostly a **variant-set artifact**: each panel matches a different ~440K subset of GN. Here we restrict to the **357,808 (chrom, pos, ref, alt) keys present in BOTH panels' GN-matched record sets**, then evaluate both panels on the **same 80 cactus founders** (verified earlier: identical sets) using the **same GN truth**.\n",
        "\n",
        "If the pipelines agree on this fair set, the Plot-5 difference is variant-coverage, not pipeline quality. If they still differ, it's a real per-(founder, SNP) call disagreement.\n",
        "\n",
        "Source: `scratch/per_record_af_fair_v3c_vs_p82.npz` from the in-script triple intersection.\n"
    ]
}
code_cell = {
    'cell_type': 'code',
    'id': 'fair6code01',
    'metadata': {},
    'execution_count': None,
    'outputs': [],
    'source': [
        "z = np.load(ROOT/'scratch/per_record_af_fair_v3c_vs_p82.npz', allow_pickle=True)\n",
        "v3c_af = z['v3c_af']; p82_af = z['p82_af']; gn_af = z['gn_af']\n",
        "keep = z['keep'].astype(bool); F = int(z['F'])\n",
        "n_total = int(z['n_total']); n_keep = int(keep.sum())\n",
        "print(f'shared (chrom,pos,ref,alt) keys: {n_total:,}   polymorphic-in-both: {n_keep:,}')\n",
        "\n",
        "edges = (np.arange(F + 2) - 0.5) / F\n",
        "fig, axes = plt.subplots(1, 3, figsize=(20, 6.2))\n",
        "panels = [\n",
        "    ('v3 cactus-only (80) vs GN', v3c_af, gn_af),\n",
        "    ('p82 (80, GN-overlap) vs GN', p82_af, gn_af),\n",
        "    ('p82 vs v3 cactus-only',     p82_af, v3c_af),\n",
        "]\n",
        "last_im = None\n",
        "for ax, (lab, y_af, x_af) in zip(axes, panels):\n",
        "    diff = y_af[keep] - x_af[keep]\n",
        "    H, xe, ye, im = ax.hist2d(x_af[keep], y_af[keep], bins=[edges, edges],\n",
        "                               cmap='viridis', norm=LogNorm(vmin=1))\n",
        "    last_im = im\n",
        "    ax.plot([0,1],[0,1],'r--', lw=0.7)\n",
        "    ax.set_xlim(0,1); ax.set_ylim(0,1); ax.set_aspect('equal')\n",
        "    if 'vs GN' in lab:\n",
        "        ax.set_xlabel(f'GN_AF (carriers/{F})')\n",
        "        ax.set_ylabel(f'panel_AF (carriers/{F})')\n",
        "    else:\n",
        "        ax.set_xlabel(f'v3 cactus-only AF (carriers/{F})')\n",
        "        ax.set_ylabel(f'p82 AF (carriers/{F})')\n",
        "    ax.set_title(\n",
        "        f'{lab}\\n'\n",
        "        f'n={n_keep:,}  mean Δ={diff.mean():+.5f}\\n'\n",
        "        f'|Δ|>0.1: {(np.abs(diff)>0.1).sum():,}  '\n",
        "        f'(>+0.1 / <-0.1 = {(diff>0.1).sum():,} / {(diff<-0.1).sum():,})',\n",
        "        fontsize=10\n",
        "    )\n",
        "fig.subplots_adjust(right=0.92, wspace=0.28)\n",
        "cax = fig.add_axes([0.94, 0.18, 0.012, 0.65])\n",
        "fig.colorbar(last_im, cax=cax, label='# SNPs (log)')\n",
        "fig.suptitle('FAIR comparison: same 80 founders, same 358K shared SNPs, same GN truth',\n",
        "             y=1.02, fontsize=12)\n",
        "plt.show()\n"
    ]
}

# Insert before the bottom-line markdown
nb['cells'] = nb['cells'][:i_bottom] + [md_cell, code_cell] + nb['cells'][i_bottom:]

# Also append a paragraph to the bottom-line cell
i_bottom_new = find_cell('148e6cea')
existing = ''.join(nb['cells'][i_bottom_new]['source']) if isinstance(nb['cells'][i_bottom_new]['source'], list) else nb['cells'][i_bottom_new]['source']
addendum = (
    "\n\n---\n\n"
    "### Plot 6 takeaway (fair comparison)\n\n"
    "On 357,808 SNPs present in BOTH v3-cactus-only and p82 GN-matched sets, evaluated on the same 80 founders against the same GN truth:\n\n"
    "- v3 cactus-only: mean Δ = +0.00207, |Δ|>0.1 in 611 records (434 over / 177 under)\n"
    "- p82:            mean Δ = +0.00203, |Δ|>0.1 in 627 records (426 over / 201 under)\n\n"
    "**The two pipelines agree at essentially the same level on shared variants.** The dramatic difference in Plot 5 was driven by the ~95K records each panel has that the other doesn't (p82 has more cactus-discovered singletons; v3 has more variants from the PG merge surviving). The PIPELINE quality (cactus build → atomization → per-founder GT call) is comparable; only the variant-set COVERAGE differs.\n"
)
nb['cells'][i_bottom_new]['source'] = existing + addendum

NB.write_text(json.dumps(nb, indent=1))
print('OK — added Plot 6')
