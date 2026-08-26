#!/usr/bin/env python
"""Build+execute the two GEMMA per-garden GWAS notebooks, in the house styles.

Requested 2026-08-25: match the two reference notebooks rather than invent a layout.
  1. notebooks/persite_gwas_manhattan_qq.ipynb  -- style of raw_manhattan_snp_vs_sv_lfmm_tile:
     ONE figure per garden per class contrast, mirror Manhattan (SNP up / class down, alternating
     chromosome shading, per-class Bonferroni dashed) on the left, QQ with both classes and their
     lambdas on the right. Manhattan and QQ in one place, as asked.
  2. notebooks/persite_gwas_dotgrid.ipynb       -- style of newpeak_dotgrid_lfmm_nonsnp:
     gene x garden dot grid of class-only ("new") peaks, gardens as columns ordered cold->hot,
     dot size = -log10 p, dot colour = TAIR/UniProt functional category, black-edged = class-only
     (SNP-blind) vs faded/white-edged = SNP also significant, right-hand bar = |variant size|.

Both read the GEMMA per-garden results (results/gemma_gwas/persite_gwas_{class}.npz), NOT the
retired in-house class_gwas_*.npz.

**One deliberate departure from `_build_persite_new_peaks_nb.py`** (now retired): that builder
called `merge_small_blocks` SEPARATELY for each class. merge_small_blocks remaps a block id to
its predecessor based on how many records the class puts in each block, so the dense SNP set and
the sparse SV set get DIFFERENT partitions -- and the new-peak logic then joins them on block id,
silently comparing non-identical intervals (this is the mechanism behind the known SV
block-merge artifact). Here the partition is the UNMERGED tiling assignment, which is fixed by
the panel's LD structure and therefore identical for every class, so a block id means the same
interval on both sides of the join.

Env: `basic` (matplotlib hangs in `plotting`). Gene annotation needs outbound HTTPS (works from
savio4 compute nodes) and degrades gracefully to symbols-only if unavailable.
"""
import os
import nbformat as nbf
from nbformat.v4 import new_notebook, new_markdown_cell, new_code_cell
from nbconvert.preprocessors import ExecutePreprocessor

ROOT = "/global/scratch/projects/fc_moilab/tbellg/kmate"
GEA = f"{ROOT}/analysis/grenenet_selection"
NBDIR = f"{GEA}/notebooks"
os.makedirs(NBDIR, exist_ok=True)

# ---------------------------------------------------------------- shared setup cell
SETUP = '''
import os, sys
import numpy as np, pandas as pd
from scipy import stats
import matplotlib as mpl, matplotlib.pyplot as plt
import matplotlib.lines as mlines
os.chdir("__ROOT__")
sys.path.insert(0, "analysis/grenenet_selection")
sys.path.insert(0, "analysis/grenenet_selection/blocks")
import lib, blocks_tiling as bt
plt.rcParams.update({"figure.dpi": 120, "font.size": 10})

RES = "analysis/grenenet_selection/r3_persite_gwas/results/gemma_gwas"
PLOTS = f"{RES}/plots"; os.makedirs(PLOTS, exist_ok=True)
CH = ["Chr1", "Chr2", "Chr3", "Chr4", "Chr5"]

def load(name):
    d = np.load(f"{RES}/persite_gwas_{name}.npz", allow_pickle=True)
    chrom = np.array([c.replace("chr", "Chr") for c in d["chrom"].astype(str)])
    pos = d["pos"].astype(np.int64)
    # UNMERGED tiling: fixed by panel LD, so identical across classes -- see builder docstring
    blk = np.asarray(bt.assign_tiling(chrom, pos, r2=0.9), dtype=object)
    return dict(chrom=chrom, pos=pos, P=d["P"], blk=blk,
                sites=[int(s) for s in d["sites"]], bio1=d["bio1"].astype(float))

SNP = load("snp")
sites = SNP["sites"]; bio1 = SNP["bio1"]
site_bio1 = {s: float(b) for s, b in zip(sites, bio1)}
site_order = sorted(sites, key=lambda s: site_bio1[s])          # cold -> hot
GENES = lib.load_genes()

CHROM_LEN = {c: int(SNP["pos"][SNP["chrom"] == c].max()) for c in CH}
OFF, _cum = {}, 0
for c in CH:
    OFF[c] = _cum; _cum += CHROM_LEN[c]
GENOME = _cum

def gpos(chrom, pos):
    return np.array([OFF[c] for c in chrom]) + pos

def genes_on(blk_id, d):
    """Gene symbols overlapping a block's span (+/-2 kb promoter)."""
    m = d["blk"] == blk_id
    if not m.any():
        return []
    ch = d["chrom"][m][0]; lo, hi = d["pos"][m].min() - 2000, d["pos"][m].max() + 2000
    g = GENES[(GENES.chrom.astype(str).str.replace("chr", "Chr") == ch) &
              (GENES.end > lo) & (GENES.start < hi)]
    col = "symbol" if "symbol" in g.columns else g.columns[-1]
    return [str(x) for x in g[col].dropna().unique()[:4]]

print(f"{len(sites)} gardens, bio1 {bio1.min():.1f}-{bio1.max():.1f} C; "
      f"snp {len(SNP['pos']):,} markers")
'''.replace("__ROOT__", ROOT)

HELPERS = '''
def lam_gc(p):
    p = np.clip(np.asarray(p, float), 1e-300, 1.0)
    p = p[np.isfinite(p)]
    return float(np.median(stats.chi2.isf(p, 1)) / stats.chi2.isf(0.5, 1))

def qq_xy(p, keep_head=8000, n_bulk=2500):
    """Observed vs expected -log10p; keep the whole significant tail, thin the dense bulk."""
    p = np.sort(np.clip(np.asarray(p, float)[np.isfinite(p)], 1e-300, 1.0))
    n = len(p)
    exp = -np.log10((np.arange(n) + 0.5) / n); obs = -np.log10(p)
    if n <= keep_head + n_bulk:
        idx = np.arange(n)
    else:
        idx = np.concatenate([np.arange(keep_head),
                              np.unique(np.geomspace(keep_head, n - 1, n_bulk).astype(int))])
    return exp[idx], obs[idx]

def leads(d, si):
    """Per block, the record with the smallest p in garden si. Returns DataFrame."""
    p = d["P"][:, si]
    ok = np.isfinite(p)
    df = pd.DataFrame(dict(blk=d["blk"][ok], chrom=d["chrom"][ok], pos=d["pos"][ok],
                           nlp=-np.log10(np.clip(p[ok], 1e-300, 1))))
    return df.loc[df.groupby("blk")["nlp"].idxmax()]

def new_peaks(cd, si, snp_lead=None):
    """Blocks whose CLASS lead clears Bonferroni while the SNP lead in the SAME block does not."""
    sl = leads(SNP, si) if snp_lead is None else snp_lead
    cl = leads(cd, si)
    bs = -np.log10(0.05 / np.isfinite(SNP["P"][:, si]).sum())
    bc = -np.log10(0.05 / np.isfinite(cd["P"][:, si]).sum())
    m = cl.merge(sl[["blk", "nlp"]], on="blk", how="left", suffixes=("", "_snp"))
    m["nlp_snp"] = m["nlp_snp"].fillna(0.0)
    m["snp_also"] = m["nlp_snp"] > bs
    return m[m["nlp"] > bc].copy(), bs, bc

SNP_DARK, SNP_LIGHT = "#3B3B3B", "#9A9A9A"
CLS_COL = {"nonsnp": ("#C1560F", "#E8A87C"), "sv": ("#2E7D32", "#88C68B")}
'''

# ---------------------------------------------------------------- notebook 1
MANHATTAN = '''
CLS = "__CLS__"
CD = load(CLS)
cdark, clight = CLS_COL[CLS]
rows = []
for s in site_order:
    si = sites.index(s)
    sp, cp = SNP["P"][:, si], CD["P"][:, si]
    npk, bs, bc = new_peaks(CD, si)
    only = npk[~npk.snp_also]

    fig = plt.figure(figsize=(19, 5.2))
    gs = fig.add_gridspec(1, 2, width_ratios=[3.5, 1], wspace=0.16)

    # ---- mirror Manhattan: SNP up, class down ----
    ax = fig.add_subplot(gs[0])
    for i, ch in enumerate(CH):
        ms, mc = SNP["chrom"] == ch, CD["chrom"] == ch
        ns = -np.log10(np.clip(sp[ms], 1e-300, 1)); nc = -np.log10(np.clip(cp[mc], 1e-300, 1))
        xs, xc = gpos(SNP["chrom"][ms], SNP["pos"][ms]), gpos(CD["chrom"][mc], CD["pos"][mc])
        ks = (ns > 1.0) | (np.arange(len(ns)) % 12 == 0)      # thin the null floor only
        kc = (nc > 1.0) | (np.arange(len(nc)) % 4 == 0)
        ax.scatter(xs[ks], ns[ks], s=2, c=(SNP_DARK if i % 2 == 0 else SNP_LIGHT),
                   rasterized=True, linewidths=0)
        ax.scatter(xc[kc], -nc[kc], s=2, c=(cdark if i % 2 == 0 else clight),
                   rasterized=True, linewidths=0)
    ax.axhline(0, c="k", lw=.7)
    ax.axhline(bs, ls="--", c="0.35", lw=.9); ax.axhline(-bc, ls="--", c=cdark, lw=.9)
    ymax = max(np.nanmax(-np.log10(np.clip(sp, 1e-300, 1))),
               np.nanmax(-np.log10(np.clip(cp, 1e-300, 1)))) * 1.15
    for _, r in only.iterrows():
        g = genes_on(r.blk, CD)
        lab = g[0] if g else r.blk
        x = OFF[r.chrom] + r.pos
        ax.axvline(x, ls=":", c="0.55", lw=.7)
        ax.annotate(lab, xy=(x, -ymax * 0.96), ha="center", va="bottom", fontsize=6,
                    color=cdark, rotation=90)
        rows.append(dict(garden=s, bio1=site_bio1[s], blk=r.blk, chrom=r.chrom, pos=int(r.pos),
                         nlp=round(r.nlp, 2), nlp_snp=round(r.nlp_snp, 2), genes=";".join(g)))
    ax.set_ylim(-ymax, ymax)
    ax.set_xticks([OFF[c] + CHROM_LEN[c] / 2 for c in CH]); ax.set_xticklabels(CH, fontsize=8)
    ax.set_ylabel(f"-log10 p  (SNP \\u2191   {CLS} \\u2193)", fontsize=9)
    ax.annotate(f"garden {s} \\u00b7 bio1 {site_bio1[s]:.1f}\\u00b0C \\u00b7 GEMMA per-garden LMM "
                f"\\u00b7 {len(only)} {CLS}-only Bonferroni new peaks",
                xy=(0.004, 0.97), xycoords="axes fraction", va="top", fontsize=8,
                bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="0.7", alpha=.85))
    for sp_ in ["top", "right"]:
        ax.spines[sp_].set_visible(False)

    # ---- QQ, both classes, lambdas in the legend ----
    axq = fig.add_subplot(gs[1])
    xs_, ys_ = qq_xy(sp); xc_, yc_ = qq_xy(cp)
    hi = max(xs_.max(), ys_.max(), xc_.max(), yc_.max())
    axq.plot([0, hi], [0, hi], c="k", lw=.8, ls="--", zorder=1)
    axq.scatter(xs_, ys_, s=4, c=SNP_DARK, rasterized=True, linewidths=0, zorder=2,
                label=f"SNP  (\\u03bb={lam_gc(sp):.2f})")
    axq.scatter(xc_, yc_, s=4, c=cdark, rasterized=True, linewidths=0, zorder=2,
                label=f"{CLS}  (\\u03bb={lam_gc(cp):.2f})")
    axq.set_xlabel("expected -log10 p", fontsize=9); axq.set_ylabel("observed -log10 p", fontsize=9)
    axq.legend(loc="upper left", frameon=False, fontsize=8, markerscale=2, handletextpad=0.3)
    for sp_ in ["top", "right"]:
        axq.spines[sp_].set_visible(False)

    fig.tight_layout()
    fig.savefig(f"{PLOTS}/manhattan_qq_{CLS}_garden{s}.png", dpi=120, bbox_inches="tight")
    plt.show(); plt.close(fig)

NEW___CLS__ = pd.DataFrame(rows)
NEW___CLS__.to_csv(f"{RES}/persite_newpeaks___CLS__.csv", index=False)
print(f"{len(NEW___CLS__)} {CLS}-only new peaks across {NEW___CLS__.garden.nunique()} gardens"
      if len(NEW___CLS__) else f"no {CLS}-only new peaks")
'''

# ---------------------------------------------------------------- notebook 2 (dot grid)
DOTGRID_DATA = '''
# Collect every class-only ("new") peak across gardens, for BOTH non-SNP and SV.
ALL = []
for CLS in ["nonsnp", "sv"]:
    CD = load(CLS)
    for s in site_order:
        si = sites.index(s)
        npk, bs, bc = new_peaks(CD, si)
        for _, r in npk.iterrows():
            for g in (genes_on(r.blk, CD) or [r.blk]):
                ALL.append(dict(cls=CLS, garden=s, bio1=site_bio1[s], blk=r.blk, gene=g,
                                chrom=r.chrom, pos=int(r.pos), nlp=float(r.nlp),
                                snp_also=bool(r.snp_also)))
ALL = pd.DataFrame(ALL)
ALL.to_csv(f"{RES}/persite_newpeak_genes.csv", index=False)
print(f"{len(ALL)} (gene, garden) cells | {ALL.gene.nunique()} genes | "
      f"class-only {int((~ALL.snp_also).sum())} vs SNP-also {int(ALL.snp_also.sum())}")
print(ALL.groupby("cls").agg(cells=("gene", "size"), genes=("gene", "nunique")).to_string())
'''

DOTGRID_SIZE = '''
# variant size (bp) per position, from the SAME panel meta the scan was built from
SIZE = {}
for _ci in range(1, 6):
    _m = np.load(f"{lib.PROJ}/panel/arch3/chr{_ci}/var_pa_231_arch3_chr{_ci}.meta.npz",
                 allow_pickle=True)
    _p = _m["pos"].astype(np.int64)
    _d = np.abs(_m["alt_len"].astype(int) - _m["ref_len"].astype(int))
    for q, dd in zip(_p, _d):                       # max per pos (multiallelic)
        k = (f"Chr{_ci}", int(q))
        if dd > SIZE.get(k, -1):
            SIZE[k] = int(dd)
ALL["size_bp"] = [SIZE.get((c, p), np.nan) for c, p in zip(ALL.chrom, ALL.pos)]
print(f"size attached for {ALL.size_bp.notna().mean():.1%} of cells; "
      f"median {ALL.size_bp.median():.0f} bp, max {ALL.size_bp.max():.0f} bp")
'''

DOTGRID_ANNOT = '''
# functional category per gene -- TAIR GO + UniProt (the mandated annotator). Degrades to
# "unclassified" if outbound HTTPS is unavailable, so the grid still renders.
CATS, ann = {}, None
try:
    import importlib.util as _ilu
    _p = "analysis/grenenet_selection/r2_gea_nonsnp/phase1_replication/annotate_genes_tair_uniprot.py"
    if not os.path.exists(_p):
        _p = "analysis/grenenet_selection/archive_gea/phase1_replication/annotate_genes_tair_uniprot.py"
    _s = _ilu.spec_from_file_location("ann", _p); ann = _ilu.module_from_spec(_s)
    _s.loader.exec_module(ann)
    print("annotator loaded:", _p)
except Exception as e:
    print("annotator unavailable ->", type(e).__name__, e)

genes = sorted(ALL.gene.unique())
if ann is not None:
    try:
        # annotate() returns columns gene / symbol / ... / `categories` (PLURAL, comma-joined,
        # empty string when nothing matched) -- not `category`.
        res = ann.annotate(genes)
        res.to_csv(f"{RES}/persite_newpeak_gene_annotation.csv", index=False)
        CATS = {r.gene: (str(r.categories).split(",")[0] if str(r.categories).strip() else
                         "unclassified") for r in res.itertuples()}
        SYM = {r.gene: (r.symbol or r.gene) for r in res.itertuples()}
        NAME = {r.gene: r.protein_name for r in res.itertuples()}
    except Exception as e:
        print(f"  annotate() failed: {type(e).__name__}: {e}")
        SYM, NAME = {}, {}
else:
    SYM, NAME = {}, {}
ALL["category"] = [CATS.get(g, "unclassified") for g in ALL.gene]
ALL["symbol"] = [SYM.get(g, g) for g in ALL.gene]
print(f"\\ncategories assigned: {(ALL.category != 'unclassified').sum()}/{len(ALL)} cells, "
      f"{ALL.loc[ALL.category != 'unclassified', 'gene'].nunique()}/{ALL.gene.nunique()} genes")
print(ALL.drop_duplicates('gene').category.value_counts().to_string())
'''

DOTGRID_PLOT = '''
CAT_ORDER = [c for c in ALL.category.value_counts().index]
PAL = plt.cm.tab10.colors + plt.cm.Set3.colors
CCOL = {c: PAL[i % len(PAL)] for i, c in enumerate(CAT_ORDER)}
xpos = {s: i for i, s in enumerate(site_order)}

for CLS in ["nonsnp", "sv"]:
    sub = ALL[ALL.cls == CLS]
    if not len(sub):
        print(f"no {CLS} new peaks to plot"); continue
    # rows ordered by |variant size| (largest at top), matching the reference dot grid
    gsz = sub.groupby("gene").size_bp.max().fillna(0).sort_values()
    gord = list(gsz.index)
    yidx = {g: i for i, g in enumerate(gord)}
    h = max(3.2, 0.23 * len(gord))
    fig, (ax, axr) = plt.subplots(1, 2, figsize=(13.5, h), sharey=True,
                                  gridspec_kw={"width_ratios": [4.2, 1], "wspace": 0.03})
    for _, r in sub.iterrows():
        ax.scatter(xpos[r.garden], yidx[r.gene], s=10 + r.nlp * 2.0,
                   c=[CCOL[r.category]], alpha=0.45 if r.snp_also else 1.0,
                   edgecolors="white" if r.snp_also else "black",
                   linewidths=0.5, zorder=3)
    ax.set_xticks(range(len(site_order)))
    ax.set_xticklabels([f"{s}\\n{site_bio1[s]:.0f}\\u00b0" for s in site_order], fontsize=6)
    ax.set_yticks(range(len(gord))); ax.set_yticklabels(gord, fontsize=6)
    ax.set_xlabel("garden (cold $\\\\rightarrow$ hot)", fontsize=9)
    ax.grid(axis="y", lw=0.3, color="0.9", zorder=0)
    ax.set_axisbelow(True)
    ax.annotate(f"{CLS}-only new peaks", xy=(0.005, 1.005), xycoords="axes fraction",
                fontsize=9, weight="bold")
    sizes = gsz.reindex(gord).fillna(0).to_numpy()
    axr.barh(range(len(gord)), np.maximum(sizes, 1), height=0.72, color="0.6")
    axr.set_xscale("log"); axr.set_xlabel("|indel/SV| size (bp)", fontsize=8)
    axr.grid(axis="x", lw=0.3, color="0.9"); axr.set_axisbelow(True)
    handles = [mlines.Line2D([], [], marker="o", ls="", ms=6, color=CCOL[c], label=c)
               for c in CAT_ORDER]
    handles += [mlines.Line2D([], [], marker="o", ls="", ms=6, mfc="0.7", mec="black",
                              label=f"{CLS}-only (SNP-blind)"),
                mlines.Line2D([], [], marker="o", ls="", ms=6, mfc="0.85", mec="white",
                              label="SNP also significant")]
    ax.legend(handles=handles, fontsize=6, loc="upper left", bbox_to_anchor=(1.28, 1.0),
              frameon=False)
    fig.savefig(f"{PLOTS}/newpeak_dotgrid_{CLS}.png", dpi=150, bbox_inches="tight")
    plt.show(); plt.close(fig)
'''


def build_manhattan():
    C = [new_markdown_cell("""# Per-garden GWAS — mirror Manhattans + QQ, garden by garden

One figure per garden: **mirror Manhattan** on the left (SNP up, non-SNP/SV down, alternating
chromosome shading, dashed line = that class's Bonferroni 0.05/n) and the **QQ for the same
garden** on the right, both classes overlaid with their genomic-inflation λ in the legend.
Layout follows `raw_manhattan_snp_vs_sv_lfmm_tile.ipynb`.

Estimator is **GEMMA v0.98.5 `-lmm 1`** (Wald) per garden, LOCO GRM, 231 founders, rank-inverse-
normalised selection coefficient. There is no cross-garden meta — each panel is an independent
scan (see `METHODS_PERSITE_GWAS.md`).

Vertical dotted line + rotated gene label = a **class-only new peak**: a tiling block whose lead
non-SNP/SV record clears Bonferroni while the lead SNP record *in the same block* does not.

> ⚠️ **Read the λ in each QQ panel with care.** λ near 1.0 does not certify these scans: stratified
> by minor-allele count the *tail* runs 5–60× hot, and ~90% of significant markers sit in the
> lowest-frequency stratum. See `persite_gwas.ipynb` §3 and `mac_calibration.csv`. Peaks below are
> a candidate set, not a calibrated hit list."""),
         new_code_cell(SETUP), new_code_cell(HELPERS)]
    for cls in ("nonsnp", "sv"):
        C.append(new_markdown_cell(f"## SNP vs {cls} — all 30 gardens, cold → hot"))
        C.append(new_code_cell(MANHATTAN.replace("__CLS__", cls)))
    return C


def build_dotgrid():
    C = [new_markdown_cell("""# Per-garden GWAS — class-only new peaks across gardens (dot grid)

Which genes does the non-SNP / SV scan flag that the SNP scan misses, and in which gardens?
One dot per (gene, garden) where that gene's tiling block carries a Bonferroni-significant
non-SNP/SV lead. Layout follows `newpeak_dotgrid_lfmm_nonsnp.ipynb`, with **gardens as columns
ordered cold → hot** instead of climate axes.

- **dot size** = −log10 p of the block's lead record
- **dot colour** = TAIR-GO + UniProt functional category
- **solid, black-edged** = class-only (the SNP-blind peaks); **faded, white-edged** = the SNP
  scan is also significant in that block, kept for context rather than filtered away
- **right-hand bar** = |alt_len − ref_len| for the lead variant, log scale; rows ordered by size

> ⚠️ Cross-garden recurrence here is **descriptive, not a tested contrast** — the cross-garden
> test was the meta this analysis deliberately dropped, and all 30 gardens share one founder
> panel and one `p0`, so their scans are correlated by construction. Low-frequency markers also
> dominate the underlying hit list (see the MAC caveat in `persite_gwas.ipynb`)."""),
         new_code_cell(SETUP), new_code_cell(HELPERS),
         new_markdown_cell("## Collect class-only peaks and their genes"),
         new_code_cell(DOTGRID_DATA),
         new_markdown_cell("## Variant size and functional category"),
         new_code_cell(DOTGRID_SIZE), new_code_cell(DOTGRID_ANNOT),
         new_markdown_cell("## The dot grid"),
         new_code_cell(DOTGRID_PLOT)]
    return C


def run(cells, out):
    nb = new_notebook(cells=cells)
    ep = ExecutePreprocessor(timeout=7200, kernel_name="basic", startup_timeout=180)
    ep.preprocess(nb, {"metadata": {"path": NBDIR}})
    with open(out, "w") as f:
        nbf.write(nb, f)
    print(f"wrote {out}")


if __name__ == "__main__":
    import sys
    which = sys.argv[1] if len(sys.argv) > 1 else "both"
    if which in ("both", "manhattan"):
        run(build_manhattan(), f"{NBDIR}/persite_gwas_manhattan_qq.ipynb")
    if which in ("both", "dotgrid"):
        run(build_dotgrid(), f"{NBDIR}/persite_gwas_dotgrid.ipynb")
