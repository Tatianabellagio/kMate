#!/usr/bin/env python
"""Build + execute the site-4 SFS-across-time notebook (user 2026-07-10): how the
allele-frequency spectrum of SNPs / small indels / SVs starts at the seedmix (p0)
and shifts across generations 1-3 at site 4 (hot pilot site). Load-only (reads
sfs_time_site4.npz / _summary.csv built by _compute_sfs_time_site4.py) -> runs in
the `basic` env.
"""
import os
import nbformat as nbf
from nbformat.v4 import new_notebook, new_markdown_cell, new_code_cell
from nbconvert.preprocessors import ExecutePreprocessor

ROOT = "/global/scratch/users/tbellg/kmate"
OUT = f"{ROOT}/analysis/grenenet_selection/notebooks/sfs_time_site4.ipynb"
os.makedirs(os.path.dirname(OUT), exist_ok=True)

md_title = """# Allele-frequency spectrum over time at site 4 — SNP / small indel / SV

Site 4 is the "hot" pilot site (`SITE_CLIMATE`). For each variant class — **SNP**, **small
indel** (≤50 bp), **SV** (>50 bp) — this tracks the distribution of alt-allele frequency
(the SFS) starting from the founding **seedmix (p0**, mean over the 8 SEEDMIX kMate reps)
through the site's evolved generations (**gen 1 → 2 → 3**, `site_gen_plot` pool means).

Built by `_compute_sfs_time_site4.py` from the compact `af_store` (per-sample AF vectors);
this notebook only loads the cache and plots. SV threshold and class definitions are
`lib.py`'s (`SV_MIN_BP=50`).

**Reading the log-x histograms below:** the leftmost bin edge is a clip floor (3e-4), so it
pools every record below that frequency (not just ones exactly at that value) — its height
tracks how much mass is piling up near zero, not a literal density at 3e-4. Exact zero (fully
lost) and exact one (fixed) are excluded from the histograms and shown separately as
`frac_lost` / `frac_fixed` in the trend panel."""

code_load = r"""
import os
import numpy as np, pandas as pd, matplotlib.pyplot as plt
from matplotlib.lines import Line2D
plt.rcParams.update({'figure.dpi':110,'font.size':9,'axes.linewidth':0.6})
G = "/global/scratch/users/tbellg/kmate/analysis/grenenet_selection"
PLOTS = f"{G}/notebooks/plots"
os.makedirs(PLOTS, exist_ok=True)
npz = np.load(f"{G}/sfs_time_site4.npz")
summary = pd.read_csv(f"{G}/sfs_time_site4_summary.csv")
COL = {"snp":"#888888","indel":"#0072B2","sv":"#D55E00"}
LABEL = {"snp":"SNP","indel":"indel","sv":"SV"}
GENS = ["p0","g1","g2","g3"]
GEN_LABEL = {"p0":"seedmix (p0)","g1":"gen 1","g2":"gen 2","g3":"gen 3"}
print(summary.to_string(index=False))
"""

code_grid = r"""
# SFS per generation, log-x (most variants are rare) -- one panel per generation,
# classes overlaid as step-density histograms.
XR = (3e-4, 1.0)
bins = np.geomspace(*XR, 40)
fig, axes = plt.subplots(1, 4, figsize=(16, 3.6), sharey=True)
for ax, g in zip(axes, GENS):
    for cls in ["snp", "indel", "sv"]:
        v = npz[f"{cls}_{g}"]
        v = v[v > 0]                       # log-x: drop exact 0 (shown separately as frac_lost)
        ax.hist(np.clip(v, *XR), bins=bins, density=True, histtype="step",
                color=COL[cls], lw=1.8 if cls == "sv" else 1.2, alpha=0.9)
    ax.set_xscale("log")
    ax.set_xlim(*XR)
    n_samp = summary.loc[(summary.gen == g) & (summary.cls == "snp"), "n_samples"].iloc[0]
    ax.set_title(f"{GEN_LABEL[g]}  (n={int(n_samp)} pools)" if g != "p0"
                 else f"{GEN_LABEL[g]}  (n=8 reps)", fontsize=9.5)
    ax.tick_params(labelsize=7)
fig.supxlabel("alt-allele frequency (log scale)", y=0.02, fontsize=10)
fig.supylabel("density", x=0.008, fontsize=10)
fig.legend(handles=[Line2D([0],[0],color=COL[c],lw=2,label=LABEL[c]) for c in ["snp","indel","sv"]],
           loc="upper right", ncol=3, fontsize=9, frameon=False, bbox_to_anchor=(0.99,1.02))
fig.suptitle("Site 4 (hot): allele-frequency spectrum, seedmix → gen 3 — SNP / indel / SV",
             fontsize=12, y=1.06)
fig.tight_layout()
fig.savefig(f"{PLOTS}/sfs_time_site4_grid.png", dpi=130, bbox_inches="tight")
plt.show()
print("saved sfs_time_site4_grid.png")
"""

code_trend = r"""
# Summary trend across generations: median AF, and the tails (fraction lost / fixed)
fig, axes = plt.subplots(1, 3, figsize=(13, 3.6))
x = np.arange(len(GENS))
for cls in ["snp", "indel", "sv"]:
    s = summary[summary.cls == cls].set_index("gen").loc[GENS]
    axes[0].plot(x, s["median"], "o-", color=COL[cls], label=LABEL[cls], lw=1.8, ms=5)
    axes[1].plot(x, s["frac_lost"] * 100, "o-", color=COL[cls], lw=1.8, ms=5)
    axes[2].plot(x, s["frac_fixed"] * 100, "o-", color=COL[cls], lw=1.8, ms=5)
for ax, ttl, ylab in zip(
        axes, ["median alt-AF", "% records lost (AF≈0)", "% records fixed (AF≈1)"],
        ["median AF", "% lost", "% fixed"]):
    ax.set_xticks(x); ax.set_xticklabels([GEN_LABEL[g] for g in GENS], rotation=20, fontsize=8)
    ax.set_title(ttl, fontsize=10); ax.set_ylabel(ylab, fontsize=9)
    ax.tick_params(labelsize=8)
axes[0].legend(fontsize=8, frameon=False)
fig.suptitle("Site 4: SFS summary statistics vs time, by class", fontsize=11, y=1.03)
fig.tight_layout()
fig.savefig(f"{PLOTS}/sfs_time_site4_trend.png", dpi=130, bbox_inches="tight")
plt.show()
print("saved sfs_time_site4_trend.png")
"""

md_take = """### Takeaway
- **All three classes start from a similar shape at the seedmix**: median alt-AF ~0.02, most
  founder-segregating variants rare (as expected for a 231-founder mixture), SV spectrum shifted
  slightly lower than SNP/indel (median 0.011 vs 0.021) — SVs are already rarer in the founder
  panel before any selection acts.
- **Across gen 1→3 the bulk (SNP, indel) spectra barely move** — same median, same shape —
  while the **fraction lost creeps up for all classes but fastest for SVs**: by gen 3, SVs are
  lost at ~0.23% of records vs ~0.07-0.09% for SNP/indel (see `frac_lost` panel), consistent with
  SVs experiencing more drift/purging pressure at this site.
- Read this alongside `SV_TEMPORAL_PURGING_SUMMARY.md`: the de-trended per-variant test found the
  genome-wide SV purging excess is ≈0 (n.s.) once frequency is properly accounted for — the
  raw SFS shift shown here is consistent with that (small, and confounded with SVs' lower starting
  frequency) rather than a new independent signal."""

nb = new_notebook(cells=[
    new_markdown_cell(md_title),
    new_code_cell(code_load),
    new_code_cell(code_grid),
    new_code_cell(code_trend),
    new_markdown_cell(md_take),
])
ep = ExecutePreprocessor(timeout=600, kernel_name="basic", startup_timeout=180)
ep.preprocess(nb, {"metadata": {"path": os.path.dirname(OUT)}})
with open(OUT, "w") as f:
    nbf.write(nb, f)
print(f"[built+executed] {OUT}")
