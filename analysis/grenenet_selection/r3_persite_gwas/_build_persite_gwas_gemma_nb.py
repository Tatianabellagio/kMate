#!/usr/bin/env python
"""Build+execute the GEMMA per-garden GWAS notebook: Manhattans + QQ plots.

The results view for the GEMMA per-garden scan (gemma_persite_gwas.py). Replaces
class_gwas_persite.ipynb / class_gwas_multitrait.ipynb, which were built on the in-house
EMMAX/P3D arm and whose multitrait content is obsolete under the per-garden-only scope.

Panels, in the order they answer questions:
  1. per-garden lambda + hit counts, all three classes
  2. QQ per garden x class -- 3 grid figures, 30 panels each
  3. MAC-stratified QQ -- THE decision panel: lambda ~1.0 in every stratum while the tail runs
     5-60x hot, worst where the hits are (see mac_calibration.csv)
  4. Manhattan per garden x class -- 3 grid figures
  5. Manhattan zoom for the gardens carrying Bonferroni hits

Plot convention (CLAUDE.md): NO chart or subplot titles anywhere. Panel identity goes in an
in-panel corner annotation in axes-fraction coords.

Load-only (results/gemma_gwas/persite_gwas_{class}.npz + mac_calibration.csv). Figures ->
results/gemma_gwas/plots/. Runs in `basic` env (matplotlib hangs in `plotting`).
"""
import os
import nbformat as nbf
from nbformat.v4 import new_notebook, new_markdown_cell, new_code_cell
from nbconvert.preprocessors import ExecutePreprocessor

ROOT = "/global/scratch/projects/fc_moilab/tbellg/kmate"
GEA = f"{ROOT}/analysis/grenenet_selection"
OUT = f"{GEA}/notebooks/persite_gwas.ipynb"   # REPLACES class_gwas_{persite,multitrait}.ipynb
os.makedirs(os.path.dirname(OUT), exist_ok=True)

setup = f'''
import os, json
import numpy as np, pandas as pd
from scipy import stats
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({{"figure.dpi": 110, "font.size": 8, "axes.grid": False}})

RES = "{GEA}/r3_persite_gwas/results/gemma_gwas"
PLOTS = f"{{RES}}/plots"; os.makedirs(PLOTS, exist_ok=True)
CLASSES = ["snp", "nonsnp", "sv"]
COL = {{"snp": "#4C72B0", "nonsnp": "#DD8452", "sv": "#55A868"}}

D = {{c: np.load(f"{{RES}}/persite_gwas_{{c}}.npz", allow_pickle=True) for c in CLASSES}}
summary = json.load(open(f"{{RES}}/persite_gwas_summary.json"))
gardens = np.array(summary["gardens"]); bio1 = np.array(summary["bio1"])
order = np.argsort(bio1)          # cold -> hot, the standard ordering in this project
print(f"{{len(gardens)}} gardens, bio1 {{bio1.min():.1f}}-{{bio1.max():.1f}} C")
for c in CLASSES:
    print(f"  {{c:7s}} {{D[c]['Z'].shape[0]:>9,d}} markers x {{D[c]['Z'].shape[1]}} gardens")
'''

cell_table = '''
rows = []
for c in CLASSES:
    P = D[c]["P"]; M = P.shape[0]; thr = 0.05 / M
    for si in range(len(gardens)):
        p = P[:, si]; ok = np.isfinite(p)
        rows.append(dict(cls=c, garden=int(gardens[si]), bio1=float(bio1[si]),
                         lam=float(np.median(stats.norm.ppf(1 - p[ok] / 2) ** 2) / stats.chi2.ppf(0.5, 1)),
                         n_bonf=int((p[ok] < thr).sum()),
                         best_nlp=float(-np.log10(np.nanmin(p[ok]))) if ok.any() else np.nan))
tab = pd.DataFrame(rows)
piv = tab.pivot(index="garden", columns="cls", values="n_bonf").reindex(gardens[order])
piv["bio1"] = bio1[order].round(1)
lam = tab.pivot(index="garden", columns="cls", values="lam").reindex(gardens[order]).round(3)
print("Bonferroni hits per garden (cold -> hot):"); print(piv.to_string())
print(f"\\nTotal garden-hits: {tab.groupby('cls').n_bonf.sum().to_dict()}")
uniq = {c: int(((D[c]["P"] < 0.05 / D[c]["P"].shape[0]).any(1)).sum()) for c in CLASSES}
print(f"UNIQUE markers hit in >=1 garden: {uniq}")
print(f"\\nlambda per garden -- median {lam.median().round(3).to_dict()}")
'''

cell_qq_grid = '''
def qq_ax(ax, p, color, label):
    p = p[np.isfinite(p)]; p = np.sort(p)
    n = len(p)
    exp = -np.log10((np.arange(1, n + 1) - 0.5) / n)
    obs = -np.log10(np.clip(p, 1e-300, 1))
    # thin the dense null end so the figure stays light without distorting the tail
    keep = np.concatenate([np.arange(0, min(n, 5000)),
                           np.arange(5000, n, max(1, n // 3000))]) if n > 8000 else np.arange(n)
    ax.plot(exp[keep], obs[keep], ".", ms=1.5, color=color, rasterized=True)
    m = max(exp.max(), obs.max())
    ax.plot([0, m], [0, m], "-", lw=0.6, color="0.5")
    lamv = np.median(stats.norm.ppf(1 - p / 2) ** 2) / stats.chi2.ppf(0.5, 1)
    ax.annotate(f"{label}\\n$\\\\lambda$={lamv:.2f}", xy=(0.04, 0.96), xycoords="axes fraction",
                va="top", ha="left", fontsize=6.5)
    ax.set_xticks([]); ax.set_yticks([])

for c in CLASSES:
    P = D[c]["P"]
    fig, axes = plt.subplots(5, 6, figsize=(13, 10.5), sharex=True, sharey=True)
    for k, si in enumerate(order):
        qq_ax(axes.flat[k], P[:, si], COL[c], f"g{gardens[si]} ({bio1[si]:.0f}C)")
    for ax in axes.flat[len(order):]:
        ax.axis("off")
    fig.supxlabel("expected $-\\\\log_{10}p$"); fig.supylabel("observed $-\\\\log_{10}p$")
    fig.tight_layout()
    fig.savefig(f"{PLOTS}/qq_grid_{c}.png", dpi=140, bbox_inches="tight")
    print(f"wrote qq_grid_{c}.png  (30 gardens, cold->hot)")
    plt.close(fig)
'''

cell_qq_mac = '''
# THE decision panel: QQ split by MAC stratum. lambda is ~1.0 in every stratum while the tail
# runs 5-60x hot -- lambda is a median statistic and never touches the quantiles hits live in.
mc = pd.read_csv(f"{RES}/mac_calibration.csv")
mac = {c: np.load(f"{RES}/mac_by_class.npz")[c] for c in CLASSES}
BINS = [(5, 7), (8, 11), (12, 22), (23, 45), (46, 10**9)]
fig, axes = plt.subplots(1, 3, figsize=(13, 4.4))
for ax, c in zip(axes, CLASSES):
    P = D[c]["P"]; m_all = mac[c]
    for (lo, hi), shade in zip(BINS, np.linspace(0.15, 0.85, len(BINS))):
        sel = (m_all >= lo) & (m_all <= hi)
        if sel.sum() < 100:
            continue
        p = P[sel].ravel(); p = p[np.isfinite(p)]
        n = len(p); p = np.sort(p)
        exp = -np.log10((np.arange(1, n + 1) - 0.5) / n)
        obs = -np.log10(np.clip(p, 1e-300, 1))
        keep = np.concatenate([np.arange(0, min(n, 4000)),
                               np.arange(4000, n, max(1, n // 2500))]) if n > 6000 else np.arange(n)
        lab = f"MAC {lo}-{hi}" if hi < 10**9 else f"MAC$\\\\geq${lo}"
        ax.plot(exp[keep], obs[keep], ".", ms=1.5, rasterized=True,
                color=plt.cm.viridis(shade), label=lab)
    m = 7.5
    ax.plot([0, m], [0, m], "-", lw=0.7, color="0.4")
    ax.annotate(c, xy=(0.04, 0.96), xycoords="axes fraction", va="top", fontsize=9, weight="bold")
    ax.set_xlabel("expected $-\\\\log_{10}p$"); ax.set_xlim(0, m)
    ax.legend(fontsize=6, loc="lower right", markerscale=4, frameon=False)
axes[0].set_ylabel("observed $-\\\\log_{10}p$")
fig.tight_layout(); fig.savefig(f"{PLOTS}/qq_by_mac_stratum.png", dpi=150, bbox_inches="tight")
print("wrote qq_by_mac_stratum.png")
plt.close(fig)
print(mc[["cls", "mac_lo", "mac_hi", "lambda_gc", "ratio_1e4", "ratio_1e6"]].to_string(index=False))
'''

cell_manhattan = '''
CH = ["chr1", "chr2", "chr3", "chr4", "chr5"]
def genome_coords(chrom, pos):
    off, cum = {}, 0
    for c in CH:
        off[c] = cum; cum += pos[chrom == c].max() if (chrom == c).any() else 0
    x = np.array([off[c] for c in chrom]) + pos
    centers = [off[c] + (pos[chrom == c].max() / 2 if (chrom == c).any() else 0) for c in CH]
    return x, centers, cum

for c in CLASSES:
    chrom = D[c]["chrom"].astype(str); pos = D[c]["pos"].astype(float); P = D[c]["P"]
    x, centers, tot = genome_coords(chrom, pos)
    band = np.isin(chrom, CH[::2])
    thr = -np.log10(0.05 / P.shape[0])
    fig, axes = plt.subplots(10, 3, figsize=(14, 18), sharex=True, sharey=True)
    for k, si in enumerate(order):
        ax = axes.flat[k]
        nlp = -np.log10(np.clip(P[:, si], 1e-300, 1))
        # keep every informative point, thin the dense null floor 20x -- plotting all ~550k
        # sub-threshold markers x 30 panels is 16M glyphs per figure and renders nothing extra
        strong = nlp > 2.0
        floor = np.where(nlp <= 2.0)[0][::20]
        keep = np.zeros(len(nlp), bool); keep[strong] = True; keep[floor] = True
        ax.scatter(x[keep & band], nlp[keep & band], s=0.7, c="0.62", rasterized=True, lw=0)
        ax.scatter(x[keep & ~band], nlp[keep & ~band], s=0.7, c=COL[c], rasterized=True, lw=0)
        ax.axhline(thr, color="crimson", lw=0.6, ls="--")
        nb = int((P[:, si] < 0.05 / P.shape[0]).sum())
        ax.annotate(f"g{gardens[si]}  {bio1[si]:.0f}$^\\\\circ$C  n={nb}", xy=(0.01, 0.94),
                    xycoords="axes fraction", va="top", fontsize=6.5)
        ax.set_xticks([]); ax.set_ylim(0.5, max(thr + 1.5, np.nanmax(nlp) * 1.05))
    for ax in axes.flat[len(order):]:
        ax.axis("off")
    fig.supylabel("$-\\\\log_{10}p$"); fig.supxlabel("genome position (chr1-chr5)")
    fig.tight_layout()
    fig.savefig(f"{PLOTS}/manhattan_grid_{c}.png", dpi=130, bbox_inches="tight")
    print(f"wrote manhattan_grid_{c}.png")
    plt.close(fig)
'''

cell_hits = '''
# the actual hit list, with MAC attached -- the column that decides whether to trust a row
recs = []
for c in CLASSES:
    P = D[c]["P"]; thr = 0.05 / P.shape[0]
    chrom = D[c]["chrom"].astype(str); pos = D[c]["pos"]
    m_all = mac[c]
    hit_i, hit_g = np.where(P < thr)
    for i, g in zip(hit_i, hit_g):
        recs.append(dict(cls=c, chrom=chrom[i], pos=int(pos[i]), garden=int(gardens[g]),
                         bio1=round(float(bio1[g]), 1), mac=int(m_all[i]),
                         maf=round(m_all[i] / 231, 4),
                         nlp=round(float(-np.log10(max(P[i, g], 1e-300))), 2)))
hits = pd.DataFrame(recs).sort_values("nlp", ascending=False)
hits.to_csv(f"{RES}/persite_bonferroni_hits.csv", index=False)
print(f"{len(hits)} garden-hits; {hits.groupby('cls').size().to_dict()}")
print(f"\\nMAC distribution of hits vs panel:")
for c in CLASSES:
    h = hits[hits.cls == c]
    print(f"  {c:7s} hits MAC<12: {(h.mac<12).mean():5.1%}  |  panel MAC<12: {(mac[c]<12).mean():5.1%}")
print("\\nTop 15 by -log10 p:")
print(hits.head(15).to_string(index=False))
'''

md_top = r"""# Per-garden GWAS on the ecotype selection coefficient — GEMMA results

30 gardens × 3 marker classes, each an **independent** GEMMA v0.98.5 univariate LMM
(`-lmm 1`, Wald) on the rank-inverse-normalized per-founder selection coefficient, with a
leave-one-chromosome-out GRM built from all classes pooled.

There is **no cross-site meta** — the all-sites question belongs to the GEA track. So a marker
significant in one garden is a per-garden result, and cross-garden recurrence below is
**descriptive, not a tested contrast**.

Methods of record: `r3_persite_gwas/METHODS_PERSITE_GWAS.md`.
"""

md_warn = r"""## ⚠️ Read the MAC-stratified QQ before interpreting any hit

Genomic-control λ is ~1.0 in every garden and every MAC stratum, which looks like clean
calibration. It is not sufficient: λ is a **median** statistic and never touches the extreme
quantiles that Bonferroni hits are drawn from. Stratified by MAC, the tail runs **5–60×** hot,
and **90% of significant SNP markers sit in MAC 5–11**, the worst-calibrated stratum.

Not all artifact — a 3–4× excess at p<1e-4 persists in clean strata (expected for a polygenic
trait plus LD), and the inflation is non-monotone, rising again at MAC ≥46 where power says real
signal belongs. Rare-end artifact superimposed on common-end signal; these numbers cannot
separate them. **Open decision: permutation null per stratum, or raise the MAC floor.**
"""

nb = new_notebook(cells=[
    new_markdown_cell(md_top),
    new_code_cell(setup),
    new_markdown_cell("## 1. Per-garden λ and Bonferroni hit counts"),
    new_code_cell(cell_table),
    new_markdown_cell("## 2. QQ per garden — 30 panels, ordered cold → hot"),
    new_code_cell(cell_qq_grid),
    new_markdown_cell(md_warn),
    new_markdown_cell("## 3. QQ by MAC stratum — the decision panel"),
    new_code_cell(cell_qq_mac),
    new_markdown_cell("## 4. Manhattan per garden — 30 panels, ordered cold → hot"),
    new_code_cell(cell_manhattan),
    new_markdown_cell("## 5. The hit list, with MAC attached"),
    new_code_cell(cell_hits),
])

ep = ExecutePreprocessor(timeout=3600, kernel_name="basic", startup_timeout=180)
ep.preprocess(nb, {"metadata": {"path": os.path.dirname(OUT)}})
with open(OUT, "w") as f:
    nbf.write(nb, f)
print(f"wrote {OUT}")
