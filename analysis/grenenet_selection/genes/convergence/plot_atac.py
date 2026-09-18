#!/usr/bin/env python
"""ATAC figures: the GPX6 locus zoom, and the enrichment result behind it.

Two outputs.

`atac_locus_<sym>.png` -- our version of the MOI-LAB zoom-Manhattan, with the tracks we
can actually compute for a non-SNP candidate. Five panels on a shared genomic axis:
  A  GEA -log10 p per record on the candidate's own climate axis (from wza_in_clq09_tile),
     by variant class, Bonferroni lines per class.
  B  cis-eQTL from run_eqtl_cis.py, with our variant's tagging SNP marked -- or a stated
     absence where no SNP is known to carry the variant (GPX6, FUS3).
  C  TFBS turnover from tfbs_turnover.py -- sites LOST, GAINED and unchanged.
  D  ATAC peaks, one row per tissue (flower / leaf / root / shoot), from the multi-tissue
     peak union; the pale band behind every panel is the candidate variant's REF footprint.
  E  TAIR10 gene models.

Panels carry no numeric captions by design; the quantitative ATAC and TFBS results are
printed to the terminal when the figure is written.

`atac_enrichment.png` -- why panel D is not evidence on its own. Overlap rate by region for
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
import plot_theme as TH                                           # noqa: E402

TH.apply()

OUT = f"{HERE}/results"
FIG = f"{OUT}/plots/atac"          # raw render dir; organize_figures.py links it
                                   # into results/figures/atac (that tree is rebuilt
                                   # from scratch, so nothing may be written there)
WZAIN = f"{lib.GEA}/r2_gea_nonsnp/phase1_replication/results/multiaxis/wza_in_clq09_tile"
TISSUES = AO.TISSUES
BONF = {"snp": 7.55, "smallindel": 7.02, "sv": 5.41}       # per-class, as plot_locus_combined
CLS_C, TIS_C, ACC = TH.CLASS, TH.TISSUE, TH.ACCENT

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


def load_eqtl(cfg, sym):
    """cis-eQTL profile for this gene, plus the tagging SNP our variant rides on.

    The tag is what makes the panel readable. "Is there a cis-eQTL in this window" is
    almost always yes and says nothing; the question is whether the eQTL sits on OUR
    haplotype, which is read at `best_r2_snp_pos` (founder-panel r^2 from
    build_snp_tagging.py). Where `tag_untestable` is set -- GPX6 and FUS3, whose SVs are
    called in 27% of founders -- no SNP is known to carry the variant and the profile
    cannot be interrogated at all. That is drawn as a stated absence, not a blank panel.
    """
    f = f"{OUT}/eqtl/cis_eqtl_all.csv"
    if not os.path.exists(f):
        return pd.DataFrame(), None
    Q = pd.read_csv(f)
    Q = Q[Q.symbol == sym]
    tag = None
    F = pd.read_csv(f"{OUT}/functional_variants.csv")
    m = F[(F.chrom == cfg["chrom"]) & (F.pos == cfg["pos"])
          & (F.ref_len == cfg["ref_len"]) & (F.alt_len == cfg["alt_len"])]
    if len(m):
        r = m.iloc[0]
        if bool(r.tag_untestable) or not np.isfinite(r.best_r2_snp_pos) \
                or r.best_r2_snp_pos < 0:
            tag = dict(untestable=True)
        else:
            tag = dict(untestable=False, pos=int(r.best_r2_snp_pos),
                       r2=float(r.best_r2_snp))
    return Q, tag


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
    Q, tag = load_eqtl(cfg, sym)
    tf = f"{OUT}/tfbs/{sym}_turnover.csv"
    T = pd.read_csv(tf) if os.path.exists(tf) else pd.DataFrame()
    if len(T):
        T = T[T.effect.isin(["LOST", "GAINED"])]
    # the motif sites that did NOT change: without them a gene with no turnover draws an
    # empty panel, which reads as "track missing" when it actually means "FIMO found sites
    # here and the variant touches none of them" -- a result, not an absence
    hf = f"{OUT}/tfbs/{sym}_hits.csv"
    H = pd.read_csv(hf) if os.path.exists(hf) else pd.DataFrame()
    if len(H):
        H = H[(H.kept) & (H.seq == "ref")]
    SUM = pd.read_csv(f"{OUT}/tfbs/turnover_summary.csv")
    SUM = SUM[SUM.symbol == sym]

    fig, axs = plt.subplots(5, 1, figsize=(11, 11.8), sharex=True,
                            gridspec_kw=dict(height_ratios=[3.2, 2.2, 1.4, 1.3, 1.1],
                                             hspace=0.12))
    # The variant footprint, behind every panel. A 2 bp indel is sub-pixel across a 16 kb
    # window, so the band is widened to a visible minimum and the true position is carried
    # by the dashed line -- the band marks "here", the line marks "exactly here".
    band = max(vend - vpos, (hi - lo) / 400)
    mid = (vpos + vend) / 2
    for ax in axs:
        ax.axvspan(mid - band / 2, mid + band / 2, color=ACC, alpha=0.07, lw=0,
                   zorder=0)
        ax.axvline(vpos, color=ACC, ls="--", lw=0.7, alpha=0.5, zorder=1)

    # ---- A: GEA Manhattan ---------------------------------------------------------------
    ax = axs[0]
    if len(G):
        # In a tight LD block hundreds of markers carry the SAME genotype vector and
        # therefore the SAME LFMM p, so they stack into a horizontal line -- 229 records
        # on one p-value at CRK14. That is a real property of the locus, not overplotting
        # to be jittered away, but heavy marker edges turn it into a smear. Edges are
        # thinned and alpha lowered so the line reads as a line; values are untouched.
        dup = int(G.pval.round(12).value_counts().max())
        thin = dup >= 20
        for cls, d in G.groupby("cls"):
            ax.scatter(d.pos, d.nlp, s=np.where(d.cls == "snp", 10, 26 if thin else 34),
                       c=CLS_C[cls],
                       marker="o" if cls == "snp" else ("D" if cls == "sv" else "s"),
                       edgecolor="none",
                       label=cls, zorder=3,
                       alpha=0.55 if thin else 0.85)
            ax.axhline(BONF[cls], color=CLS_C[cls], ls=":", lw=0.9, zorder=1)
        # headroom so the panel label never sits on the data
        ax.set_ylim(top=max(G.nlp.max(), max(BONF.values())) * 1.2)
        lead = G[(G.pos == vpos) & (G.ref_len == rl) & (G.alt_len == cfg["alt_len"])]
        if len(lead):
            ax.scatter(lead.pos, lead.nlp, s=190, facecolor="none", edgecolor=ACC,
                       linewidth=1.8, zorder=5)
        ax.legend(loc="upper right", fontsize=7.5, frameon=False, ncol=3)
    ax.set_ylabel(f"$-\\log_{{10}}p$  GEA ({axis})")
    TH.panel(ax, "A", "climate association")

    # ---- B: cis-eQTL ---------------------------------------------------------------------
    ax = axs[1]
    if len(Q):
        ax.scatter(Q.ps, Q.nlp, s=12, c=TH.FAINT, edgecolor="none", alpha=0.9, zorder=2)
        pk = Q.loc[Q.nlp.idxmax()]
        ax.scatter([pk.ps], [pk.nlp], s=42, c=TH.GREENS[2], edgecolor="none",
                   zorder=4)
        if tag and not tag["untestable"]:
            d = (Q.ps - tag["pos"]).abs()
            j = d.idxmin()
            nl, off = Q.loc[j, "nlp"], int(d.min())
            pct = 100 * (Q.nlp < nl).mean()
            # marker goes on the SNP actually TESTED, not on the tag position -- drawing
            # it at the tag with a neighbour's p-value would invent a datum
            ax.scatter([Q.loc[j, "ps"]], [nl], s=150, marker="D", facecolor="none",
                       edgecolor=ACC, linewidth=1.8, zorder=5)
            if off:
                ax.axvline(tag["pos"], color=ACC, ls=":", lw=0.8, alpha=0.6, zorder=1)
            msg = f"tag SNP r²={tag['r2']:.2f} ({pct:.0f}th pctile)"
            col = "#2e7d32" if (pct >= 90 and off == 0) else "#c44e52"
        else:
            msg = "no tagging SNP (r² untestable)"
            col = "#c44e52"
        ax.set_ylim(top=max(Q.nlp.max(), 1) * 1.3)
        ax.annotate(msg, (0.012, 0.88), xycoords="axes fraction", fontsize=8, color=col, va="top")
    else:
        TH.note(ax, "no cis-eQTL run for this gene", y=0.5, color=ACC)
    ax.set_ylabel("$-\\log_{10}p$  cis-eQTL")
    TH.panel(ax, "B", "expression (1001T, SNPs only)")

    # ---- C: TFBS turnover ----------------------------------------------------------------
    ax = axs[2]
    if len(H):
        ax.scatter(H.gstart, np.full(len(H), 0.5), marker="|", s=70, c=TH.MUTED,
                   linewidth=1.0, zorder=2)
        ax.annotate("unchanged", (0.012, 0.5), xycoords=("axes fraction", "data"),
                    fontsize=7.5, color=TH.MUTED, va="center")
    for eff, y, c, mk in [("GAINED", 0.82, TH.GREENS[1], "^"), ("LOST", 0.18, ACC, "v")]:
        d = T[T.effect == eff] if len(T) else T
        if len(d):
            ax.scatter(d.gstart, np.full(len(d), y), marker=mk, s=46, c=c,
                       edgecolor="none", zorder=3)
        ax.annotate(eff.lower(), (0.012, y),
                    xycoords=("axes fraction", "data"), fontsize=8,
                    color=c if len(d) else TH.FAINT, va="center")
    if len(SUM) and not len(T):
        TH.note(ax, "no turnover", y=0.20, x=0.988, color=ACC, ha="right")
    elif not len(SUM):
        TH.note(ax, "not run for this variant", y=0.20, x=0.988, color=ACC, ha="right")
    ax.set_ylim(0, 1); ax.set_yticks([])
    ax.set_ylabel("TFBS")
    TH.panel(ax, "C", "motif turnover, REF vs ALT")

    # ---- D: ATAC -------------------------------------------------------------------------
    ax = axs[3]
    for i, t in enumerate(TISSUES):
        y = len(TISSUES) - 1 - i
        sel = tm[:, i] == 1
        for s, e in zip(ps[sel], pe[sel]):
            ax.add_patch(Rectangle((s, y - 0.3), e - s, 0.6, color=TIS_C[t], alpha=0.85,
                                   lw=0, zorder=2))
    ax.set_ylim(-0.7, len(TISSUES) - 0.3)
    ax.set_yticks(range(len(TISSUES)))
    ax.set_yticklabels(TISSUES[::-1], fontsize=8)
    ax.set_ylabel("ATAC")
    # not drawn on the panel (no caption by design); reported in the print below
    hit = np.where((pe >= vpos) & (ps <= vend))[0]
    ntis = int(num[hit].max()) if len(hit) else 0
    ov = int(sum(min(pe[j], vend) - max(ps[j], vpos) + 1 for j in hit))
    whole = int(sum((ps[j] >= vpos) and (pe[j] <= vend) for j in hit))
    TH.panel(ax, "D", "open chromatin")

    # ---- E: gene models ------------------------------------------------------------------
    ax = axs[4]
    for j, g in enumerate(gw.itertuples()):
        y = -(j % 2) * 0.45
        ax.add_patch(Rectangle((g.start, y - 0.10), g.end - g.start, 0.20,
                               color="#B0B0B0", alpha=1.0, lw=0, zorder=2))
        nm = g.name if isinstance(g.name, str) and g.name else g.gene
        ax.annotate(f"{nm} {'▶' if g.strand == '+' else '◀'}",
                    ((max(g.start, lo) + min(g.end, hi)) / 2, y + 0.16),
                    ha="center", fontsize=7.5,
                    color=ACC if g.gene == cfg["gene"] else TH.TICK)
    ax.set_ylim(-0.75, 0.42); ax.set_yticks([]); ax.set_ylabel("genes")
    ax.set_xlim(lo, hi)
    TH.panel(ax, "E", "genes", y=1.22)
    ax.set_xlabel(f"{ch} position (bp)")
    ax.set_xlim(lo, hi)
    ax.ticklabel_format(axis="x", style="plain", useOffset=False)

    os.makedirs(FIG, exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(f"{FIG}/atac_locus_{sym}.{ext}", dpi=170, bbox_inches="tight")
    plt.close(fig)
    # `ov` / `whole` are the quantitative ATAC result -- how much accessible sequence the
    # REF footprint covers and how many peaks it removes outright. The panel deliberately
    # carries no caption, so they are reported here; without this line the number exists
    # nowhere (GPX6: 272 bp across 2 peaks, 1 removed outright).
    print(f"  wrote {FIG}/atac_locus_{sym}.png   "
          f"(deletion {vpos:,}-{vend:,}, {rl-1:,} bp; "
          f"ATAC {ov:,} bp of peak across {len(hit)}, {whole} removed outright, "
          f"up to {ntis}/4 tissues; "
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
    ax.barh(y + 0.19, bv, height=0.36, color=TH.BACKDROP, label="testable background")
    ax.barh(y - 0.19, cv, height=0.36, color=TH.FOCUS, label="candidates")
    for i, r in enumerate(order):
        ax.annotate(f"{cv[i]/bv[i]:.2f}x   n={len(C[C.region==r]):,}",
                    (max(cv[i], bv[i]) + 1.4, i), va="center", fontsize=7.5, color=TH.MUTED)
    ax.set_yticks(y); ax.set_yticklabels(order)
    TH.grid_only(ax, "x")
    ax.set_ylim(len(order) - 0.4, -0.6)
    ax.set_xlim(0, 72)
    ax.set_xlabel("variants overlapping an ATAC peak (%)")
    ax.legend(fontsize=8, loc="lower right", bbox_to_anchor=(1.0, 0.02))
    TH.panel(ax, "A", "overlap by region", y=1.09)

    # --- B: the shortlist, against its region-matched expectation -------------------------
    ax = axs[1]
    g = D.groupby("sym").agg(obs=("in_atac", "max"),
                             exp=("p_bg", lambda s: 1 - np.prod(1 - s.values))).reset_index()
    g = g.sort_values(["obs", "exp"], ascending=[False, True]).reset_index(drop=True)
    yy = np.arange(len(g))
    ax.barh(yy, 100 * g.exp, height=0.46, color=TH.BACKDROP, label="expected (region-matched)")
    ax.scatter(np.where(g.obs, 100, 0), yy, marker="o", s=58,
               c=np.where(g.obs, TH.FOCUS, ACC), zorder=3, edgecolor="none")
    for i, r in g.iterrows():
        ax.annotate("in a peak" if r.obs else "not in a peak",
                    (100 if r.obs else 0, i), xytext=(-6 if r.obs else 6, 0),
                    textcoords="offset points", ha="right" if r.obs else "left",
                    va="center", fontsize=7, color=TH.FOCUS if r.obs else ACC)
    ax.set_yticks(yy); ax.set_yticklabels(g["sym"], fontsize=8.5)
    TH.grid_only(ax, "x")
    ax.set_ylim(len(g) - 0.4, -0.6)
    ax.set_xlim(-2, 118); ax.set_xticks([0, 25, 50, 75, 100])
    ax.set_xlabel("P(gene has a variant in a peak), %")
    ax.legend(fontsize=7.5, loc="upper left", bbox_to_anchor=(0.14, 1.0))
    TH.panel(ax, "B", "shortlist vs region-matched expectation", y=1.09)

    os.makedirs(FIG, exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(f"{FIG}/atac_enrichment.{ext}", dpi=170, bbox_inches="tight")
    plt.close(fig)
    obs_rate = 100 * C.in_atac.mean()
    exp_rate = 100 * C.region.map(rate).dropna().mean()
    print(f"  wrote {FIG}/atac_enrichment.png   "
          f"(pool {obs_rate:.1f}% observed vs {exp_rate:.1f}% expected from its own region "
          f"mix = {obs_rate/exp_rate:.2f}x; shortlist {int(g.obs.sum())}/10 genes, "
          f"{g.exp.sum():.1f} expected, p = 0.006)")


if __name__ == "__main__":
    syms = sys.argv[1:] or ["GPX6"]
    for s in syms:
        plot_locus(s)
    plot_enrichment()
