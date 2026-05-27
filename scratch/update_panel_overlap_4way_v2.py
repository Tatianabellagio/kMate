"""Patch panel_overlap_4way.ipynb plot 5:
- Filter to polymorphic-in-panel records only (0 < carrier_count < F)
- Replace hexbin with hist2d aligned to discrete carrier-AF values (kills banding)
- Tighter titles + single shared colorbar (no overlap)
"""
import json
from pathlib import Path

NB = Path('/carnegie/nobackup/scratch/tbellagio/hapfire_sv/panel_overlap_4way.ipynb')
nb = json.loads(NB.read_text())

def find_cell(id_prefix):
    for i, c in enumerate(nb['cells']):
        if c.get('id', '').startswith(id_prefix): return i
    raise KeyError(id_prefix)

# Plot 5 markdown
i_md = find_cell('bd018fa8')
nb['cells'][i_md]['source'] = [
    "## Plot 5 — Per-record carrier-AF: panel vs GN, POLYMORPHIC-in-panel only\n",
    "\n",
    "For each panel, plot per-record carrier-AF (carriers / F_overlap) against GN's carrier-AF on the matched (chrom, pos, REF, ALT) intersection. **Filtered to records where the panel itself is polymorphic** (0 < panel_carriers < F) — i.e., we drop records where the variant is monomorphic in the panel subset, since those are uninformative for sub-panel-vs-GN comparisons.\n",
    "\n",
    "**Discrete-grid aware**: both axes are integer multiples of 1/F (F=231 or 80). Plot uses `hist2d` with bin edges centered on `k/F` so each cell counts exact (panel_carriers, gn_carriers) pairs. Eliminates the hex-grid aliasing bands seen with `hexbin`.\n",
    "\n",
    "**Why filter**: a SNP carried only by PG founders has 0 carriers in the cactus-only subset → carrier_AF=0/80=0. Without filtering, all such PG-only SNPs pile up at y=0 and look like systematic under-calling (which they are, but only because we removed the founders that would have carried them). Filtering isolates the apples-to-apples question: where the panel has variation, does it agree with GN on which founders carry it?\n"
]

# Plot 5 code: 3-panel, polymorphic filter, hist2d with discrete-aligned bins, single colorbar
i = find_cell('7d94d7a4')
nb['cells'][i]['source'] = [
    "panel_specs = [\n",
    "    ('v3 (all 231)',           ROOT/'scratch/per_record_af_v3.npz'),\n",
    "    ('v3 cactus-only (80)',    ROOT/'scratch/per_record_af_v3_cactus_only.npz'),\n",
    "    ('p82 (82 cactus, no PG)', ROOT/'scratch/per_record_af_p82.npz'),\n",
    "]\n",
    "missing = [p for _, p in panel_specs if not p.exists()]\n",
    "if missing:\n",
    "    print('Missing:'); [print(' ', p) for p in missing]\n",
    "else:\n",
    "    fig, axes = plt.subplots(1, 3, figsize=(20, 6.2))\n",
    "    last_im = None\n",
    "    for ax, (lab, p) in zip(axes, panel_specs):\n",
    "        z = np.load(p, allow_pickle=True)\n",
    "        gn_af = z['gn_af']; panel_af = z['panel_af']\n",
    "        F = int(z['F_overlap'])\n",
    "        # Recover integer carrier counts (panel_af = carriers/F is exact)\n",
    "        panel_carriers = np.round(panel_af * F).astype(int)\n",
    "        gn_carriers    = np.round(gn_af    * F).astype(int)\n",
    "        # Polymorphic-in-panel filter\n",
    "        keep = (panel_carriers > 0) & (panel_carriers < F)\n",
    "        n_total = len(panel_carriers); n_keep = int(keep.sum())\n",
    "        gn_af_f = gn_af[keep]; panel_af_f = panel_af[keep]\n",
    "        diff_f  = panel_af_f - gn_af_f\n",
    "        # Discrete-aligned bin edges: F+1 cells, each centered on k/F\n",
    "        edges = (np.arange(F + 2) - 0.5) / F\n",
    "        H, xe, ye, im = ax.hist2d(gn_af_f, panel_af_f, bins=[edges, edges],\n",
    "                                   cmap='viridis', norm=LogNorm(vmin=1))\n",
    "        last_im = im\n",
    "        ax.plot([0,1],[0,1], 'r--', lw=0.7)\n",
    "        ax.set_xlim(0,1); ax.set_ylim(0,1); ax.set_aspect('equal')\n",
    "        ax.set_xlabel(f'GN_AF (carriers/{F})')\n",
    "        ax.set_ylabel(f'panel_AF (carriers/{F})')\n",
    "        ax.set_title(\n",
    "            f'{lab}\\n'\n",
    "            f'kept {n_keep:,}/{n_total:,} polymorphic   '\n",
    "            f'mean Δ={diff_f.mean():+.4f}\\n'\n",
    "            f'|Δ|>0.1 : {(np.abs(diff_f)>0.1).sum():,}   '\n",
    "            f'(>+0.1 / <-0.1 = {(diff_f>0.1).sum():,} / {(diff_f<-0.1).sum():,})',\n",
    "            fontsize=10\n",
    "        )\n",
    "    # single colorbar on the right\n",
    "    fig.subplots_adjust(right=0.92, wspace=0.28)\n",
    "    cax = fig.add_axes([0.94, 0.18, 0.012, 0.65])\n",
    "    fig.colorbar(last_im, cax=cax, label='# SNPs (log)')\n",
    "    fig.suptitle('Per-record carrier-AF: panel vs GrENE-Net (polymorphic-in-panel, BEFORE EM)',\n",
    "                 y=1.02, fontsize=12)\n",
    "    plt.show()\n"
]
nb['cells'][i]['outputs'] = []
nb['cells'][i]['execution_count'] = None

NB.write_text(json.dumps(nb, indent=1))
print('OK — patched plot 5 in', NB)
