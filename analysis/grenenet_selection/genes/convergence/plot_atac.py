#!/usr/bin/env python
"""ATAC figures: the GPX6 locus zoom, and the enrichment result behind it.

Two outputs.

`atac_locus_<sym>.png` -- our version of the MOI-LAB zoom-Manhattan, with the tracks we
can actually compute for a non-SNP candidate. Four panels on a shared genomic axis:
  A  GEA -log10 p per record on the candidate's own climate axis (from wza_in_clq09_tile),
     by variant class, Bonferroni lines per class.
  B  TFBS turnover from tfbs_turnover.py -- sites LOST and GAINED on the ALT allele.
  C  ATAC peaks, one row per tissue (flower / leaf / root / shoot), from the multi-tissue
     peak union; the pale band behind every panel is the candidate variant's REF footprint.
  D  TAIR10 gene models.

`atac_enrichment.png` -- why panel C is not evidence on its own. Overlap rate by region for
candidates against the region-annotated testable background, and the shortlist against its
region-matched expectation. The pool sits on the diagonal (1.04x overall); only the
shortlist leaves it.

No titles anywhere -- panels are identified by an in-panel corner annotation (house rule).

Usage:  PY plot_atac.py [SYM ...]      default GPX6
env: kmate. Compute node.
"""
from __future__ import annotations
import os, sys
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

HERE = os.path.dirname(os.path.abspath(__file__))
GEA_DIR = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, GEA_DIR); sys.path.insert(0, HERE)
sys.path.insert(0, f"{GEA_DIR}/r2_gea_nonsnp/phase1_replication")
import lib                                                        # noqa: E402
import atac_overlap as AO                                         # noqa: E402

OUT = f"{HERE}/results"
FIG = f"{OUT}/plots/atac"          # raw render dir; organize_figures.py links it
                                   # into results/figures/atac (that tree is rebuilt
                                   # from scratch, so nothing may be written there)
WZAIN = f"{lib.GEA}/r2_gea_nonsnp/phase1_replication/results/multiaxis/wza_in_clq09_tile"
TISSUES = AO.TISSUES
BONF = {"snp": 7.55, "smallindel": 7.02, "sv": 5.41}       # per-class, as plot_locus_combined
CLS_C = {"snp": "#b8b8b8", "smallindel": "#4c72b0", "sv": "#c44e52"}
TIS_C = {"flower": "#b07aa1", "leaf": "#59a14f", "root": "#9c755f", "shoot": "#4e79a7"}

# candidate -> the variant to zoom on; axis is the record's own best GEA axis
CANDIDATES = {
    # alt_len is not decoration: 2.14% of arch3 positions carry several records, and
    # circling on chrom+pos alone outlines the wrong allele (memory
    # `panel-multiallelic-pos-key-trap`). AT4G13200 has two records at 7,669,175.
    "GPX6":      dict(gene="AT4G11600", chrom="Chr4", pos=7011705, ref_len=1225,
                      alt_len=61, axis="pc1"),
    "AT4G13200": dict(gene="AT4G13200", chrom="Chr4", pos=7669175, ref_len=1,
                      alt_len=3, axis=None),
    "AT2G30000": dict(gene="AT2G30000", chrom="Chr2", pos=12805666, ref_len=1,
                      alt_len=2, axis=None),
    "FUS3":      dict(gene="AT3G26790", chrom="Chr3", pos=9856500, ref_len=3,
                      alt_len=1, axis=None),
    "CRK14":     dict(gene="AT4G23220", chrom="Chr4", pos=12157244, ref_len=64,
                      alt_len=1, axis=None),
}


def gea_axis(cfg):
    if cfg["axis"]:
        return cfg["axis"]
    G = pd.read_csv(f"{OUT}/gea_pool.csv")
    m = G[(G.chrom == cfg["chrom"]) & (G.pos == cfg["pos"])
          & (G.ref_len == cfg["ref_len"]) & (G.alt_len == cfg["alt_len"])]
    return m.best_axis.iloc[0] if len(m) else None


def load_gea(cfg, lo, hi, axis):
    frames = []
    for cls in ("snp", "smallindel", "sv"):
        f = f"{WZAIN}/lfmm_{cls}_gen9_{axis}.csv"
        if not os.path.exists(f):
            continue
        d = pd.read_csv(f)
        d = d[(d.chrom == cfg["chrom"]) & (d.pos >= lo) & (d.pos <= hi) & (d.MAF > 0.05)].copy()
        d["cls"] = cls
        frames.append(d)
    if not frames:
        return pd.DataFrame()
    g = pd.concat(frames, ignore_index=True)
    g["nlp"] = -np.log10(g.pval.clip(lower=1e-300))
    return g


def plot_locus(sym, pad=8000):
    cfg = CANDIDATES[sym]
    ch, vpos, rl = cfg["chrom"], cfg["pos"], cfg["ref_len"]
    assert "alt_len" in cfg, "alt_len is required -- see the multiallelic note above"
    vend = vpos + rl - 1
    lo, hi = vpos - pad, vend + pad
    axis = gea_axis(cfg)

    genes = lib.load_genes()
    gw = genes[(genes.chrom == ch) & (genes.end >= lo) & (genes.start <= hi)]
    peaks = AO.load_peaks()
    ps, pe, num, tm = peaks[ch]
    k = (pe >= lo) & (ps <= hi)
    ps, pe, num, tm = ps[k], pe[k], num[k], tm[k]

    G = load_gea(cfg, lo, hi, axis) if axis else pd.DataFrame()
    tf = f"{OUT}/tfbs/{sym}_turnover.csv"
    T = pd.read_csv(tf) if os.path.exists(tf) else pd.DataFrame()
    if len(T):
        T = T[T.effect.isin(["LOST", "GAINED"])]

    fig, axs = plt.subplots(4, 1, figsize=(11, 8.6), sharex=True,
                            gridspec_kw=dict(height_ratios=[3.0, 1.5, 1.3, 1.1], hspace=0.12))
    # The variant footprint, behind every panel. A 2 bp indel is sub-pixel across a 16 kb
    # window, so the band is widened to a visible minimum and the true position is carried
    # by the dashed line -- the band marks "here", the line marks "exactly here".
    band = max(vend - vpos, (hi - lo) / 400)
    mid = (vpos + vend) / 2
    for ax in axs:
        ax.axvspan(mid - band / 2, mid + band / 2, color="#c44e52", alpha=0.13, lw=0,
                   zorder=0)
        ax.axvline(vpos, color="#c44e52", ls="--", lw=0.7, alpha=0.55, zorder=1)

    # ---- A: GEA Manhattan ---------------------------------------------------------------
    ax = axs[0]
    if len(G):
        for cls, d in G.groupby("cls"):
            ax.scatter(d.pos, d.nlp, s=np.where(d.cls == "snp", 11, 34), c=CLS_C[cls],
                       marker="o" if cls == "snp" else ("D" if cls == "sv" else "s"),
                       edgecolor="none" if cls == "snp" else "k", linewidth=0.4,
                       label=f"{cls} (n={len(d)})", zorder=3, alpha=0.85)
            ax.axhline(BONF[cls], color=CLS_C[cls], ls=":", lw=0.9, zorder=1)
        lead = G[(G.pos == vpos) & (G.ref_len == rl) & (G.alt_len == cfg["alt_len"])]
        if len(lead):
            ax.scatter(lead.pos, lead.nlp, s=190, facecolor="none", edgecolor="#c44e52",
                       linewidth=1.8, zorder=5)
        ax.legend(loc="upper right", fontsize=7.5, frameon=False, ncol=3)
    ax.set_ylabel(f"$-\\log_{{10}}p$  GEA ({axis})")
    ax.annotate("A   climate association", (0.012, 0.93), xycoords="axes fraction",
                fontsize=9, fontweight="bold", va="top")

    # ---- B: TFBS turnover ----------------------------------------------------------------
    ax = axs[1]
    if len(T):
        for eff, y, c, mk in [("LOST", 0.68, "#c44e52", "v"), ("GAINED", 0.28, "#2e7d32", "^")]:
            d = T[T.effect == eff]
            if not len(d):
                continue
            ax.scatter(d.gstart, np.full(len(d), y), marker=mk, s=46, c=c,
                       edgecolor="k", linewidth=0.3, zorder=3)
            ax.annotate(f"{eff.lower()}  n={len(d)}", (0.012, y + 0.14),
                        xycoords=("axes fraction", "data"), fontsize=8, color=c, va="center")
    ax.set_ylim(0, 1); ax.set_yticks([])
    ax.set_ylabel("TFBS")
    ax.annotate("B   motif turnover, REF vs ALT", (0.012, 0.93), xycoords="axes fraction",
                fontsize=9, fontweight="bold", va="top")

    # ---- C: ATAC -------------------------------------------------------------------------
    ax = axs[2]
    for i, t in enumerate(TISSUES):
        y = len(TISSUES) - 1 - i
        ax.axhline(y, color="#eeeeee", lw=0.6, zorder=0)
        sel = tm[:, i] == 1
        for s, e in zip(ps[sel], pe[sel]):
            ax.add_patch(Rectangle((s, y - 0.3), e - s, 0.6, color=TIS_C[t], alpha=0.85,
                                   lw=0, zorder=2))
    ax.set_ylim(-0.7, len(TISSUES) - 0.3)
    ax.set_yticks(range(len(TISSUES)))
    ax.set_yticklabels(TISSUES[::-1], fontsize=8)
    ax.set_ylabel("ATAC")
    # "overlaps a peak" undersells a deletion: report how much accessible sequence the REF
    # footprint actually covers, and whether any peak is removed outright. GPX6 spans two.
    hit = np.where((pe >= vpos) & (ps <= vend))[0]
    ntis = int(num[hit].max()) if len(hit) else 0
    ov = int(sum(min(pe[j], vend) - max(ps[j], vpos) + 1 for j in hit))
    whole = int(sum((ps[j] >= vpos) and (pe[j] <= vend) for j in hit))
    if not len(hit):
        msg = "no ATAC peak overlaps the variant"
    elif rl > 50:                                   # a deletion: it removes peak sequence
        msg = (f"deletion covers {ov:,} bp of peak across {len(hit)} peak"
               f"{'s' if len(hit) > 1 else ''}"
               + (f", {whole} removed outright" if whole else "")
               + f"; up to {ntis} of 4 tissues")
    else:
        msg = f"variant sits inside a peak called in {ntis} of 4 tissues"
    ax.annotate(f"C   open chromatin — {msg}",
                (0.012, 0.93), xycoords="axes fraction", fontsize=9, fontweight="bold",
                va="top")

    # ---- D: gene models ------------------------------------------------------------------
    ax = axs[3]
    for j, g in enumerate(gw.itertuples()):
        y = -(j % 2) * 0.45
        ax.add_patch(Rectangle((g.start, y - 0.10), g.end - g.start, 0.20,
                               color="#37474f", alpha=0.85, lw=0, zorder=2))
        nm = g.name if isinstance(g.name, str) and g.name else g.gene
        ax.annotate(f"{nm} {'▶' if g.strand == '+' else '◀'}",
                    ((max(g.start, lo) + min(g.end, hi)) / 2, y + 0.16),
                    ha="center", fontsize=7.5,
                    fontweight="bold" if g.gene == cfg["gene"] else "normal",
                    color="#c44e52" if g.gene == cfg["gene"] else "#37474f")
    ax.set_ylim(-0.75, 0.42); ax.set_yticks([]); ax.set_ylabel("genes")
    ax.set_xlim(lo, hi)
    ax.annotate("D   genes", (0.012, 1.14), xycoords="axes fraction", fontsize=9,
                fontweight="bold", va="top")
    ax.set_xlabel(f"{ch} position (bp)")
    ax.set_xlim(lo, hi)
    ax.ticklabel_format(axis="x", style="plain", useOffset=False)

    os.makedirs(FIG, exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(f"{FIG}/atac_locus_{sym}.{ext}", dpi=170, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {FIG}/atac_locus_{sym}.png   "
          f"(deletion {vpos:,}-{vend:,}, {rl-1:,} bp; peak in {ntis}/4 tissues; "
          f"TFBS lost {int((T.effect=='LOST').sum()) if len(T) else 0} / "
          f"gained {int((T.effect=='GAINED').sum()) if len(T) else 0})")


# ------------------------------------------------------------------------------------------------
def plot_enrichment():
    C = pd.read_csv(f"{OUT}/atac_overlap.csv")
    BR = pd.read_csv(f"{OUT}/atac_enrich_background.csv.gz")
    rate = BR.groupby("region").in_atac.mean()
    SY = {"AT4G11600": "GPX6", "AT5G48410": "GLR1.3", "AT5G48412": "AT5G48412",
          "AT4G13200": "AT4G13200", "AT2G30000": "AT2G30000", "AT3G26790": "FUS3",
          "AT4G23220": "CRK14", "AT2G01180": "LPP1", "AT3G51450": "SSL7",
          "AT5G23210": "SCPL34"}
    D = C[C.gene.isin(SY)].copy(); D["sym"] = D.gene.map(SY); D["p_bg"] = D.region.map(rate)

    order = ["5'UTR", "promoter", "3'UTR", "intergenic", "CDS", "intron", "TE"]
    fig, axs = plt.subplots(1, 2, figsize=(12.6, 4.8), gridspec_kw=dict(width_ratios=[1.55, 1]))
    fig.subplots_adjust(wspace=0.42)

    # --- A: region by region, candidates vs the region-annotated background ---------------
    ax = axs[0]
    y = np.arange(len(order))
    cv = [100 * C[C.region == r].in_atac.mean() for r in order]
    bv = [100 * rate.get(r, np.nan) for r in order]
    ax.barh(y + 0.19, bv, height=0.36, color="#c9c9c9", label="testable background")
    ax.barh(y - 0.19, cv, height=0.36, color="#4c72b0", label="candidates")
    for i, r in enumerate(order):
        ax.annotate(f"{cv[i]/bv[i]:.2f}x   n={len(C[C.region==r]):,}",
                    (max(cv[i], bv[i]) + 1.4, i), va="center", fontsize=7.5, color="#555555")
    ax.set_yticks(y); ax.set_yticklabels(order)
    ax.set_ylim(len(order) - 0.4, -1.25)                 # blank top row for the annotation
    ax.set_xlim(0, 72)
    ax.set_xlabel("variants overlapping an ATAC peak (%)")
    ax.legend(fontsize=8, frameon=False, loc="lower right", bbox_to_anchor=(1.0, 0.02))
    ax.annotate("A   the pool is explained by region alone — 21.6% observed "
                "vs 20.8% expected = 1.04x",
                (0.0, -1.0), xycoords=("axes fraction", "data"), xytext=(4, 0),
                textcoords="offset points", fontsize=8.8, va="center", ha="left")

    # --- B: the shortlist, against its region-matched expectation -------------------------
    ax = axs[1]
    g = D.groupby("sym").agg(obs=("in_atac", "max"),
                             exp=("p_bg", lambda s: 1 - np.prod(1 - s.values))).reset_index()
    g = g.sort_values(["obs", "exp"], ascending=[False, True]).reset_index(drop=True)
    yy = np.arange(len(g))
    ax.barh(yy, 100 * g.exp, height=0.46, color="#d6d6d6", label="expected (region-matched)")
    ax.scatter(np.where(g.obs, 100, 0), yy, marker="o", s=58,
               c=np.where(g.obs, "#2e7d32", "#c44e52"), zorder=3, edgecolor="k", linewidth=0.4)
    for i, r in g.iterrows():
        ax.annotate("in a peak" if r.obs else "not in a peak",
                    (100 if r.obs else 0, i), xytext=(-6 if r.obs else 6, 0),
                    textcoords="offset points", ha="right" if r.obs else "left",
                    va="center", fontsize=7, color="#2e7d32" if r.obs else "#c44e52")
    ax.set_yticks(yy); ax.set_yticklabels(g["sym"], fontsize=8.5)
    ax.set_ylim(len(g) - 0.4, -1.5)                      # blank top row for the annotation
    ax.set_xlim(-2, 118); ax.set_xticks([0, 25, 50, 75, 100])
    ax.set_xlabel("P(gene has a variant in a peak), %")
    ax.legend(fontsize=7.5, frameon=False, loc="lower right", bbox_to_anchor=(1.0, -0.03))
    ax.annotate(f"B   shortlist: {int(g.obs.sum())}/10 genes, 3.9 expected, p = 0.006",
                (0.0, -1.2), xycoords=("axes fraction", "data"), xytext=(4, 0),
                textcoords="offset points", fontsize=8.8, va="center", ha="left")

    os.makedirs(FIG, exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(f"{FIG}/atac_enrichment.{ext}", dpi=170, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {FIG}/atac_enrichment.png")


if __name__ == "__main__":
    syms = sys.argv[1:] or ["GPX6"]
    for s in syms:
        plot_locus(s)
    plot_enrichment()
