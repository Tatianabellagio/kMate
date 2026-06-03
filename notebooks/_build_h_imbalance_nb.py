#!/usr/bin/env python3
"""Builder for notebooks/H_IMBALANCE_p231.ipynb.

Defines the notebook cells (markdown + code), writes the .ipynb by hand (no
nbformat dependency), then EXECUTES the code cells so the committed notebook is
exactly the code that produced the PNG. Run with the `basic` env (matplotlib):

    /global/home/users/tbellg/miniforge3/envs/basic/bin/python notebooks/_build_h_imbalance_nb.py
"""
import json, os

ROOT = "/global/scratch/users/tbellg/kmate"
NB_PATH = f"{ROOT}/notebooks/H_IMBALANCE_p231.ipynb"

MD_INTRO = r"""# Effect of k-mer filtering + bubble weighting on the founder vector $h$ (p231 panel)

The 231-founder panel mixes **80 long-read** (cactus-assembly) founders with **151
short-read** (PanGenie-genotyped) founders. Long-read assemblies realize many more
private / discriminating k-mers per founder, so an unweighted k-mer EM **over-credits
the long-read founders** — the recovered $h$ puts too much mass on them. We correct
this in two stages, the rows of the grid below:

| row | approach | k-mer matrix | EM weight |
|---|---|---|---|
| **raw** | no filter, no correction | unfiltered `kmer_pa_p231` | uniform ($\omega_k\equiv1$) |
| **private k-mer filter** | drop singleton/invariant k-mers | `kmer_pa_p231_filt2` | uniform ($\omega_k\equiv1$) |
| **private k-mer filter + bubble weighting** | + per-bubble de-replication (production) | `kmer_pa_p231_filt2` | $\omega_k=1/m_b$ |

**Columns** are three p231 simulations (Chr1, seed s42), all with known simulation
truth $h$: `n231_g0` (uniform over all 231 founders), `n50_g0` (50 random founders),
and a **short-read-heavy n=50 pool** (45 short-read + 5 long-read; long-read truth
mass 0.10) that stresses the bias — long-read founders are nearly absent, yet the
unweighted EM still credits them. In each panel: **dark green = long-read est, light
green = short-read est**, gray = truth. For the n=50 columns the founders are ordered
by **true frequency**, so the selected (present) founders sit at the right and the
spurious mass placed on absent founders is visible at the left.
"""

CODE_SETUP = r'''
import os, json, csv
import numpy as np
import matplotlib
matplotlib.use("Agg")                 # headless build; Jupyter uses inline automatically
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

ROOT  = "/global/scratch/users/tbellg/kmate"
P231  = f"{ROOT}/benchmarks/p231/results"
SIMS  = f"{ROOT}/benchmarks/p231/sims"
PLOTS = f"{ROOT}/notebooks/plots/h_imbalance"
os.makedirs(PLOTS, exist_ok=True)

# Long-read (cactus) vs short-read (PG) founder split, membership by founder ID.
_split   = json.load(open(f"{ROOT}/data/founder_split_cactus_pg.json"))
LONGREAD = set(map(str, _split["cactus"]))          # 80 long-read founders

# Greenish palette.
C_LONG, C_SHORT, C_TRUTH, C_LVL = "#006d2c", "#a1d76a", "0.45", "#b30000"
LAB_LONG, LAB_SHORT = "long-read founders", "short-read founders"

# (row label, sim-arm dir/tag)
APPROACHES = [
    ("raw\n(no filter, $\\omega\\equiv1$)",                              "raw_raw"),
    ("private k-mer filter\n($\\omega\\equiv1$)",                        "filt2u_raw"),
    ("private k-mer filter\n+ bubble weighting ($\\omega_k{=}1/m_b$)",   "filt2mb_raw"),
]
# (col label, regime, sort-by-truth?)   regime n231_g0 keeps index order (uniform truth)
COLUMNS = [
    ("n231_g0\n(uniform, all 231)",              "n231_g0",       False),
    ("n50 long-read-heavy\n(40 LR + 10 SR)",     "n50_cactheavy", True),
    ("n50 short-read-heavy\n(45 SR + 5 LR)",     "n50_pgheavy",   True),
]

def _sim_dir(regime):
    return f"{SIMS}/cov10_{regime}_s42_hotspots_p231_chr1"

def load_sim_h(arm_tag, regime):
    f = f"{P231}/kmate_global_{arm_tag}/{regime}/p231_{arm_tag}_{regime}_cov10_s42.h_per_chrom.npz"
    if not os.path.exists(f):
        return None, None
    d = np.load(f, allow_pickle=True)
    return d["founders"].astype(str), d["Chr1"].astype(float)

def load_sim_truth(founders, regime):
    th = {}
    for r in csv.DictReader(open(f"{_sim_dir(regime)}/pool_weights.tsv"), delimiter="\t"):
        th[str(r["founder"])] = float(r["weight"])
    return np.array([th.get(f, 0.0) for f in founders])

def longmask(founders):
    return np.array([f in LONGREAD for f in founders])

print("setup ok | long-read founders:", len(LONGREAD))
'''

MD_GRID = r"""## The grid — over-credit on long-read founders and its correction

Top **raw** row: dark-green (long-read) points sit above the gray truth — the EM
over-credits long-read founders. The **private k-mer filter** removes the worst
singleton tags; adding **bubble weighting $\omega_k=1/m_b$** (bottom, production)
flattens the residual over-credit so green tracks the truth. Each title reports the
long-read $h$-mass (truth $\to$ estimate), the RMSE, and the spurious mass on *absent*
founders.
"""

CODE_GRID = r'''
GREY = "#4d4d4d"
plt.rcParams.update({
    "text.color": GREY, "axes.labelcolor": GREY, "axes.titlecolor": GREY,
    "xtick.color": GREY, "ytick.color": GREY, "font.family": "sans-serif",
    # paper style: grey gridlines, no spines
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.spines.left": False, "axes.spines.bottom": False,
    "axes.grid": True, "grid.color": "#dddddd", "grid.linewidth": 0.8,
    "axes.axisbelow": True, "xtick.bottom": False, "ytick.left": False,
})

fig, axes = plt.subplots(len(APPROACHES), len(COLUMNS),
                         figsize=(4.8*len(COLUMNS), 3.7*len(APPROACHES)),
                         squeeze=False, sharey=True)

for ci, (clab, regime, by_truth) in enumerate(COLUMNS):
    for ri, (rlab, arm_tag) in enumerate(APPROACHES):
        ax = axes[ri, ci]
        fo, h = load_sim_h(arm_tag, regime)
        if fo is None:
            ax.text(0.5, 0.5, "pending\n(run 07d / 08)", ha="center", va="center",
                    transform=ax.transAxes, color="0.5")
            ax.set_xticks([]); ax.set_yticks([]);
            if ci == 0: ax.set_ylabel(rlab + "\n\n$h$", fontsize=9, color=GREY)
            if ri == 0: ax.annotate(clab, xy=(0.5, 1.16), xycoords="axes fraction",
                                    ha="center", va="bottom", fontsize=11, color=GREY)
            continue
        truth = load_sim_truth(fo, regime)
        lm = longmask(fo)
        # order founders along x: by true frequency (n50) or native index (n231 uniform)
        order = np.argsort(truth, kind="stable") if by_truth else np.arange(len(fo))
        x = np.arange(len(fo))
        t_o, h_o, lm_o = truth[order], h[order], lm[order]
        # truth (grey dots: 0 for absent founders, 1/n for the selected ones)
        ax.scatter(x, t_o, s=6, color=C_TRUTH, zorder=1)
        # estimates (dark green long-read over light green short-read)
        ax.scatter(x[~lm_o], h_o[~lm_o], marker="x", s=15, lw=.8, color=C_SHORT, zorder=2)
        ax.scatter(x[lm_o],  h_o[lm_o],  marker="x", s=15, lw=.9, color=C_LONG,  zorder=3)
        # metric: RMSE
        rmse = float(np.sqrt(np.mean((h - truth)**2)))
        ax.set_title(f"RMSE {rmse:.4f}", fontsize=10, color=GREY)
        ax.tick_params(colors=GREY)
        if ci == 0: ax.set_ylabel(rlab + "\n\n$h$", fontsize=9, color=GREY)
        if ri == len(APPROACHES) - 1: ax.set_xlabel("founder (ordered by true freq)" if by_truth else "founder index")
        if ri == 0: ax.annotate(clab, xy=(0.5, 1.16), xycoords="axes fraction",
                                ha="center", va="bottom", fontsize=11, color=GREY)

handles = [Line2D([], [], marker='x', ls='', color=C_LONG,  label=LAB_LONG),
           Line2D([], [], marker='x', ls='', color=C_SHORT, label=LAB_SHORT),
           Line2D([], [], marker='o', ls='', color=C_TRUTH, label="truth")]
fig.tight_layout()
fig.legend(handles=handles, loc="lower center", ncol=3, fontsize=9,
           bbox_to_anchor=(0.5, -0.02), frameon=False)

out = f"{PLOTS}/h_imbalance_grid_p231.png"
fig.savefig(out, dpi=130, bbox_inches="tight")
print("saved", out)
plt.show()   # inline display in Jupyter
'''

MD_OUTRO = r"""## Notes

- PNG → `notebooks/plots/h_imbalance/h_imbalance_grid_p231.png`, synced to the paper via
  `paper_svs/sync_figures.sh`.
- Data arms: sims `raw` from `benchmarks/p231/scripts/07d_run_kmate_raw_p231.sh`; the
  short-read-heavy (and cactus-heavy) n50 pools from `08_build_run_skewed_p231.sh`
  (wgsim over the per-founder Chr1 haplotypes — no VISOR). `filt2u`/`filt2mb` arms for the
  random regimes were already on disk.
- To switch column 3 to the **long-read-heavy** pool, change the COLUMNS regime to
  `n50_cactheavy`. All Chr1, single seed s42.
"""

def md(src):  return {"cell_type": "markdown", "metadata": {}, "source": src.splitlines(keepends=True)}
def code(src):
    s = src.strip("\n") + "\n"
    return {"cell_type": "code", "metadata": {}, "execution_count": None,
            "outputs": [], "source": s.splitlines(keepends=True)}

cells = [md(MD_INTRO), code(CODE_SETUP), md(MD_GRID), code(CODE_GRID), md(MD_OUTRO)]
nb = {"cells": cells,
      "metadata": {"kernelspec": {"display_name": "Python 3 (basic)", "language": "python", "name": "python3"},
                   "language_info": {"name": "python", "version": "3"}},
      "nbformat": 4, "nbformat_minor": 5}

with open(NB_PATH, "w") as fh:
    json.dump(nb, fh, indent=1)
print("wrote", NB_PATH)

ns = {}
for c in cells:
    if c["cell_type"] == "code":
        exec("".join(c["source"]), ns)
print("\nALL CELLS EXECUTED OK")
