#!/usr/bin/env python
"""Build + execute notebooks/09_wza_manhattan.ipynb (run in the `basic` env).

Standalone WZA block-level view of the two-stage SV climate-GEA: the GrENE-net
WZA (Booker weighted-Z + SNP-number spline correction, local wza_script.py) run
over the hapFIRE LD blocks (no MAF filter), for both Stage-1 statistics
(A=Δp endpoint, B=selection coefficient). Big, readable block-Manhattan + the
gene-annotated top blocks. Reads build_wza.py outputs in results/.../gea/wza/.

  BPY=/global/home/users/tbellg/miniforge3/envs/basic/bin/python
  $BPY analysis/grenenet_gea/_build_wza_nb.py
"""
import os
import nbformat as nbf
from nbconvert.preprocessors import ExecutePreprocessor

PROJ = "/global/scratch/users/tbellg/kmate"
NBDIR = f"{PROJ}/analysis/grenenet_gea/notebooks"
OUT = f"{NBDIR}/09_wza_manhattan.ipynb"

md_intro = r"""# WZA block-Manhattan — SV climate-GEA on hapFIRE LD blocks

The per-SV two-stage permutation p-values are noisy; **WZA** aggregates them to the
phase-1 **hapFIRE LD blocks** (1.05 M SNPs → 16,674 blocks), exactly as the GrENE-net
project did — Booker weighted-Z (weight = MAF·(1−MAF)) + the SNP-number spline
correction → one `Z_pVal` per block. Run here with **no MAF filter** (all SVs; low-MAF
SVs are auto-downweighted), for both Stage-1 statistics:

- **A — `dp`**: gen-3 endpoint Δp = p₃ − p₀
- **B — `scoef`**: gen0→3 logit-slope selection coefficient

Small `Z_pVal` = block whose SVs are more climate-associated than the genome-wide
background (the spline correction is genomic-control-like, so this tests **outlier
blocks**, not the polygenic background itself).
"""

code_manhattan = r"""
import pandas as pd, numpy as np
import matplotlib.pyplot as plt
WZ = "/global/scratch/users/tbellg/kmate/results/grenenet_gea/gea/wza"
order = [f"Chr{i}" for i in range(1, 6)]
titles = {"dp": "A — endpoint Δp (gen3 − p0)",
          "scoef": "B — selection coefficient (gen0→3 trajectory)"}
fig, axes = plt.subplots(2, 1, figsize=(16, 9), sharex=True)
for ax, s in zip(axes, ("dp", "scoef")):
    w = pd.read_csv(f"{WZ}/wza_{s}_bio1.csv"); w = w[w.Z_pVal.notna()].copy()
    w["chrom"] = "Chr" + w.block.str.split("_").str[0]
    w = w[w.chrom.isin(order)]; w["chrom"] = pd.Categorical(w.chrom, order, ordered=True)
    w = w.sort_values(["chrom", "mid_pos"]); off, cum, ticks = {}, 0, []
    for c in order:
        off[c] = cum; cmax = w.loc[w.chrom == c, "mid_pos"].max()
        ticks.append(cum + cmax / 2); cum += cmax + 1e6
    w["g"] = w.mid_pos + w.chrom.map(off).astype(float)
    w["mlp"] = -np.log10(w.Z_pVal.clip(lower=1e-12))
    nb = len(w); bonf = -np.log10(0.05 / nb)
    for i, c in enumerate(order):
        sd = w[w.chrom == c]
        ax.scatter(sd.g, sd.mlp, s=22, color=["#2c5aa0", "#7fb3d5"][i % 2],
                   edgecolor="none", alpha=.9)
    ax.axhline(bonf, color="r", lw=1.2, ls="--", label="Bonferroni 0.05")
    ax.axhline(-np.log10(0.05), color="gray", lw=1, ls=":", label="p = 0.05")
    t = pd.read_csv(f"{WZ}/wza_{s}_bio1.top.csv").head(6)
    t["chrom"] = "Chr" + t.block.str.split("_").str[0]
    for _, r in t.iterrows():
        gx = r.mid_pos + off.get(r.chrom, 0)
        lab = r.gene_name if isinstance(r.gene_name, str) else r.gene
        ax.annotate(f"{lab}\n({r.dir})", (gx, -np.log10(max(r.Z_pVal, 1e-12))),
                    fontsize=8, ha="center", xytext=(0, 7), textcoords="offset points")
    ax.set_xticks(ticks); ax.set_xticklabels(order, fontsize=12)
    ax.set_ylabel("-log10(WZA Z_pVal)", fontsize=12)
    nsig = (w.Z_pVal < 0.05).sum(); nbonf = (w.Z_pVal < 0.05 / nb).sum()
    ax.set_title(f"{titles[s]}   |   {nb:,} LD blocks, {nsig} at p<0.05 "
                 f"(expected {int(0.05*nb)}), {nbonf} at Bonferroni", fontsize=11)
    ax.legend(fontsize=10, loc="upper right"); ax.margins(x=0.01)
fig.suptitle("WZA block-level SV climate-GEA on hapFIRE LD blocks "
             "(bio1, all SVs, no MAF filter)", fontsize=14, y=0.995)
plt.tight_layout(); plt.show()
"""

code_tables = r"""
# Top blocks per statistic (gene-annotated; dir = mean Δp direction, + = up-in-warm).
for s in ("dp", "scoef"):
    t = pd.read_csv(f"{WZ}/wza_{s}_bio1.top.csv")
    print(f"\n===== {s}: top blocks by Z_pVal =====")
    cols = ["block","chrom","mid_pos","SNPs","Z","Z_pVal","dir","n_up","n_dn","gene_name","flower_locus"]
    print(t[cols].head(12).to_string(index=False))
"""

code_qq = r"""
# QQ of block Z_pVal (vs uniform): A hugs the line (polygenic, no outliers);
# B lifts off in the tail (a few real outlier blocks).
fig, ax = plt.subplots(1, 2, figsize=(12, 4))
for a, s in zip(ax, ("dp", "scoef")):
    p = pd.read_csv(f"{WZ}/wza_{s}_bio1.csv").Z_pVal.dropna().values
    obs = -np.log10(np.sort(p)); exp = -np.log10(np.linspace(1/len(p), 1, len(p)))
    a.plot(exp, obs, ".", ms=3, color="#2c5aa0"); mx = exp.max()
    a.plot([0, mx], [0, mx], "r-", lw=.8); a.set_title(s)
    a.set_xlabel("expected -log10(p)"); a.set_ylabel("observed")
plt.tight_layout(); plt.show()
"""

cells = [
    nbf.v4.new_markdown_cell(md_intro),
    nbf.v4.new_markdown_cell("## Block-Manhattan (A = Δp, B = selection coefficient)"),
    nbf.v4.new_code_cell(code_manhattan.strip()),
    nbf.v4.new_markdown_cell("## QQ of block-level Z_pVal"),
    nbf.v4.new_code_cell(code_qq.strip()),
    nbf.v4.new_markdown_cell("## Top blocks (gene-annotated)"),
    nbf.v4.new_code_cell(code_tables.strip()),
]
nb = nbf.v4.new_notebook(cells=cells, metadata={"kernelspec":
        {"name": "python3", "display_name": "Python 3"}})

if __name__ == "__main__":
    os.makedirs(NBDIR, exist_ok=True)
    ep = ExecutePreprocessor(timeout=600, kernel_name="python3", startup_timeout=180)
    ep.preprocess(nb, {"metadata": {"path": NBDIR}})
    with open(OUT, "w") as f:
        nbf.write(nb, f)
    print("wrote + executed", OUT)
