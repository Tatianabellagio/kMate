#!/usr/bin/env python
"""CARK8/CARK9 (AT2G30740 / AT2G30730) high-resolution locus figure — a cam5-style
combined GEA + founder-haplotype panel, adapted for the SV-only new-peak that
flagged these genes in the multiaxis LFMM WZA (newpeak_dotgrid_lfmm_sv.ipynb).

Why this locus, and why the wide window
---------------------------------------
Block Chr2_5291 (Chr2:13,092,592-13,098,713) overlaps only CARK9+CARK8, but the
Bonferroni-clearing SV that made the block significant is a ~5.1 kb insertion at
Chr2:13,127,635 -- ~29 kb DOWNSTREAM, intergenic between AT2G30800 and AT2G30810,
tagged into the CARK LD block. So the honest "where is the signal" view spans the
genes THROUGH the distal SV (Chr2:13,090,000-13,137,000, ~47 kb). The founder-LD
triangle answers whether that distal SV is genuinely linked to the CARK block.

Three stacked panels, shared genomic x-axis (LocusZoom style):
  A  LFMM GEA Manhattan (bio1 axis) -- driven by the wza_in cohort variants where
     lfmm_p / dp live; -log10 p, colour = dp (final - p0), marker = SNP/indel/SV;
     per-class genome-wide Bonferroni lines; droplines at class-Bonferroni hits.
     Gene-model track on top (bodies + strand; CARK8/9 highlighted).
  B  Founder panel (231 accessions) -- ordered cold->hot by home bio1, left strip
     = home bio1 origin; ALT ticks across the window (SNP grey / indel orange /
     SV red diamond). Over 47 kb every founder is a unique haplotype, so this is
     the per-founder view, not a collapse.
  C  Founder-LD triangle (r^2) on the common (MAF>=0.05) variants + the aggregated
     lead SV; rank-spaced; connector dropping from the lead SV.

Manhattan variants (cohort GEA) and panel genotypes (founder VCF) are DIFFERENT
sets aligned by position -- as in cam5, on purpose.

env: kmate (python + pysam + matplotlib + scipy)
"""
from __future__ import annotations
import os, sys
import numpy as np
import pandas as pd
import pysam
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyArrow
from matplotlib.collections import LineCollection
import matplotlib.colors as mcolors

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
import lib

HERE = os.path.dirname(os.path.abspath(__file__))
CHROM = "Chr2"
LO, HI = 13_090_000, 13_137_000                      # ~47 kb: CARK genes -> distal SV
AXIS = "bio1"
REGION_VCF = f"{HERE}/cark_region.vcf.gz"
BIO_CSV = ("/global/scratch/users/tbellg/gea_grene-net/key_files/"
           "1001g_regmap_grenet_ecotype_info_corrected_bioclim_2024May16.csv")
WZAIN = f"{lib.GEA}/r2_gea_nonsnp/phase1_replication/results/multiaxis/wza_in_clq09_tile"
GROUP_MEANS = f"{lib.GEA}/common/results/group_means.npz"
GENES_HL = {"AT2G30730": "CARK9", "AT2G30740": "CARK8"}   # highlighted
LEAD_POS = 13_127_635                                 # the 5.1 kb insertion

# per-class genome-wide Bonferroni (0.05/n, MAF>0.05) -- the bars the peak cleared
BONF = {"snp": 7.606, "smallindel": 7.131, "sv": 5.742}
CLSMARK = {"snp": ("o", 22), "smallindel": ("s", 26), "sv": ("D", 40)}
CLSNAME = {"snp": "SNP", "smallindel": "indel", "sv": "SV"}


# ---------------------------------------------------------------- data: GEA manhattan
def variant_class(ref_len, alt_len):
    d = abs(int(alt_len) - int(ref_len))
    return "snp" if d == 0 else ("smallindel" if d < 50 else "sv")


def load_gea():
    """Concat the disjoint GEA classes in the window; attach dp from group_means."""
    frames = []
    for cls in ("snp", "smallindel", "sv"):
        d = pd.read_csv(f"{WZAIN}/lfmm_{cls}_gen9_{AXIS}.csv")
        d = d[(d.chrom == CHROM) & (d.pos >= LO) & (d.pos <= HI) & (d.MAF > 0.05)].copy()
        d["cls"] = cls
        frames.append(d)
    g = pd.concat(frames, ignore_index=True)
    g["nlp"] = -np.log10(g["pval"].clip(lower=1e-300))

    # colour = CLIMATE divergence dp_hc = mean over gen1-3 of (hot - cold) allele freq.
    # NB: a plain net (final - p0) is ~0 for a climate SV that is antagonistically
    # sorted (up in hot / down in cold) while being purged everywhere over time -- it
    # would render white and hide the very signal LFMM detects. The hot-vs-cold
    # divergence is what the bio1 GEA actually tests, so we colour by that.
    z = np.load(GROUP_MEANS, allow_pickle=True)
    key = pd.DataFrame({"chrom": z["chrom"], "pos": z["pos"], "ref_len": z["ref_len"],
                        "alt_len": z["alt_len"]})
    div = np.zeros(len(key))
    for ggen in (1, 2, 3):
        div += (z[f"hot_g{ggen}"] - z[f"cold_g{ggen}"])
    key["dp"] = div / 3.0                                  # mean hot-cold over gens
    key = key[(key.chrom == CHROM) & (key.pos >= LO) & (key.pos <= HI)].copy()
    g = g.merge(key[["chrom", "pos", "ref_len", "alt_len", "dp"]],
                on=["chrom", "pos", "ref_len", "alt_len"], how="left")
    return g


# ---------------------------------------------------------------- data: founder panel
def load_panel():
    vcf = pysam.VariantFile(REGION_VCF)
    samp = [int(s) for s in vcf.header.samples]
    recs, G = [], []
    for r in vcf.fetch(CHROM, LO, HI):
        rl, al = len(r.ref), len(r.alts[0])
        recs.append((r.pos, rl, al, variant_class(rl, al)))
        G.append([1 if r.samples[s].get("GT", (None,))[0] == 1 else 0 for s in vcf.header.samples])
    G = np.array(G, dtype=np.int8).T                      # founders x variants
    V = pd.DataFrame(recs, columns=["pos", "ref_len", "alt_len", "cls"])
    bio1 = pd.read_csv(BIO_CSV).set_index("ecotypeid")["bio1"].to_dict()
    b1 = np.array([bio1.get(e, np.nan) for e in samp])
    return np.array(samp), b1, V, G


# ---------------------------------------------------------------- data: founder LD
def founder_ld(V, G, maf=0.05):
    """r^2 among common panel variants + the aggregated lead SV (any long INS at
    LEAD_POS collapsed to one locus). Returns positions and r^2 matrix."""
    fmaf = G.mean(0)
    fmaf = np.minimum(fmaf, 1 - fmaf)
    keep = np.where((fmaf >= maf) & (V.pos.values != LEAD_POS))[0]
    pos = V.pos.values[keep].astype(float)
    M = G[:, keep].astype(float)
    # aggregated lead SV locus
    sv_cols = np.where((V.pos.values == LEAD_POS) &
                       (np.abs(V.alt_len.values - V.ref_len.values) >= 3000))[0]
    if len(sv_cols):
        sv = (G[:, sv_cols].sum(1) > 0).astype(float)
        pos = np.append(pos, LEAD_POS)
        M = np.column_stack([M, sv])
    order = np.argsort(pos)
    pos, M = pos[order], M[:, order]
    Mc = M - M.mean(0)
    sd = Mc.std(0)
    sd[sd == 0] = np.nan
    R = (Mc.T @ Mc) / (len(M) * np.outer(sd, sd))
    return pos, R ** 2


# ---------------------------------------------------------------- plotting helpers
def diverging(vals, mid, vmax=None):
    v = np.asarray(vals, float)
    if vmax is None:
        vmax = np.nanmax(np.abs(v - mid)) or 1.0
    norm = mcolors.TwoSlopeNorm(vmin=mid - vmax, vcenter=mid, vmax=mid + vmax)
    return norm


BLUERED = mcolors.LinearSegmentedColormap.from_list(
    "bluered", ["#2166ac", "#b3b3b3", "#b2182b"])
DPCMAP = mcolors.LinearSegmentedColormap.from_list(
    "dp", ["#2166ac", "#f7f7f7", "#b2182b"])
LDCMAP = mcolors.LinearSegmentedColormap.from_list(
    "ld", ["#ffffff", "#e0e0f0", "#b3a2c8", "#8073ac", "#542788", "#2d004b"])


def gene_track(ax, y0, dy):
    genes = lib.load_genes()
    gs = genes[(genes.chrom == CHROM) & (genes.end >= LO) & (genes.start <= HI)]
    for _, g in gs.iterrows():
        hl = g.gene in GENES_HL
        col = "#b2182b" if hl else "#555555"
        ax.add_patch(Rectangle((g.start, y0 - dy / 2), g.end - g.start, dy,
                               facecolor=col, edgecolor="none",
                               alpha=0.95 if hl else 0.55, zorder=3, clip_on=False))
        # strand arrow
        xa = g.end if g.strand == "+" else g.start
        ax.annotate("", xy=(xa + (700 if g.strand == "+" else -700), y0),
                    xytext=(xa, y0),
                    arrowprops=dict(arrowstyle="->", color=col, lw=0.9), zorder=4,
                    annotation_clip=False)
        lab = GENES_HL.get(g.gene, g.gene)
        ax.text((g.start + g.end) / 2, y0 + dy * 0.9, lab, ha="center", va="bottom",
                fontsize=6.5 if hl else 5.2, style="italic",
                color=col, fontweight="bold" if hl else "normal", clip_on=False)


# ---------------------------------------------------------------- build figure
def main():
    gea = load_gea()
    samp, b1, V, G = load_panel()
    ldpos, R2 = founder_ld(V, G)

    # significant GEA variants (clear their own class Bonferroni) -> droplines
    sig = gea[gea.apply(lambda r: r.nlp >= BONF[r.cls], axis=1)]
    sigx = sorted(sig.pos.unique())
    print(f"GEA variants in window: {len(gea)} ; class-Bonferroni hits: {len(sig)} "
          f"at {sigx}")
    print(f"founder panel: {G.shape[0]} accessions x {G.shape[1]} variants ; "
          f"LD set: {len(ldpos)} loci")

    dpmax = np.nanpercentile(np.abs(gea.dp.dropna()), 98) if gea.dp.notna().any() else 0.05
    b1mid = float(np.nanmean(b1))

    # founder order: cold at TOP -> hot at bottom (matplotlib y grows upward, so
    # coldest gets the highest row index)
    order = np.argsort(np.nan_to_num(b1, nan=1e9))       # coldest first
    nF = len(samp)
    yrow = {int(i): nF - 1 - k for k, i in enumerate(order)}

    fig = plt.figure(figsize=(13, 17))
    gs = fig.add_gridspec(3, 1, height_ratios=[3.0, 6.0, 2.2], hspace=0.06)
    axM = fig.add_subplot(gs[0])
    axH = fig.add_subplot(gs[1], sharex=axM)
    axL = fig.add_subplot(gs[2], sharex=axM)

    # ---- shaded context bands (CARK block + lead SV) on all panels
    def band(ax, lo, hi, c="#cdd1d5", a=0.45):
        ax.axvspan(lo, hi, color=c, alpha=a, zorder=0, lw=0)
    for ax in (axM, axH):
        band(ax, 13_092_592, 13_098_713)                 # CARK LD block
        band(ax, LEAD_POS - 400, LEAD_POS + 400, c="#f4c6c0", a=0.7)   # lead SV
        for x in sigx:
            ax.axvline(x, ls="--", color="grey", lw=0.5, alpha=0.6, zorder=1)

    # ================= Panel A: Manhattan =================
    dnorm = diverging([0], 0, dpmax)                      # shared dp norm
    for cls, sub in gea.groupby("cls"):
        mk, sz = CLSMARK[cls]
        axM.scatter(sub.pos, sub.nlp, c=sub.dp.fillna(0), cmap=DPCMAP, norm=dnorm,
                    marker=mk, s=sz, edgecolors="white", linewidths=0.3, zorder=3,
                    label=CLSNAME[cls])
    for cls, y in [("sv", BONF["sv"]), ("snp", BONF["snp"])]:
        axM.axhline(y, ls="--", color="#888", lw=0.8)
        axM.text(HI, y, f"{CLSNAME[cls]} Bonferroni (0.05/n)", ha="right",
                 va="bottom", fontsize=6.5, color="#666")
    # lead SV: emphasise + annotate
    lead = gea[gea.pos == LEAD_POS].iloc[0]
    axM.scatter([LEAD_POS], [lead.nlp], marker="D", s=95, facecolor=DPCMAP(dnorm(lead.dp)),
                edgecolor="#b2182b", linewidths=1.4, zorder=5)
    axM.annotate(f"5.1 kb insertion\n(SV, p=3×10⁻⁸, MAF {lead.MAF:.2f})",
                 xy=(LEAD_POS, lead.nlp), xytext=(LEAD_POS - 9000, lead.nlp - 1.4),
                 fontsize=7, color="#7d1a12", ha="right",
                 arrowprops=dict(arrowstyle="->", color="#7d1a12", lw=0.9))
    axM.set_ylabel(r"$-\log_{10}(p)$  LFMM (bio1)", fontsize=9)
    ytop = max(9.3, gea.nlp.max() + 1.6)
    axM.set_ylim(-0.4, ytop)
    gene_track(axM, ytop - 0.55, 0.34)                   # gene track above all data
    for sp in ("top", "right"):
        axM.spines[sp].set_visible(False)
    axM.tick_params(labelbottom=False)
    # dp colourbar
    smdp = plt.cm.ScalarMappable(cmap=DPCMAP, norm=diverging([0], 0, dpmax))
    cb = fig.colorbar(smdp, ax=axM, fraction=0.02, pad=0.005)
    cb.set_label(r"$\Delta p$ climate (hot $-$ cold, mean gen 1–3)", fontsize=7)
    cb.ax.tick_params(labelsize=6.5)
    # marker legend
    from matplotlib.lines import Line2D
    h = [Line2D([], [], marker=CLSMARK[c][0], ls="", mfc="grey", mec="white",
                ms=6, label=CLSNAME[c]) for c in ("snp", "smallindel", "sv")]
    axM.legend(handles=h, loc="upper left", frameon=False, fontsize=7,
               title="variant", title_fontsize=7.5)

    # ================= Panel B: founder haplotype ticks =================
    # left origin strip (home bio1)
    bnorm = diverging(b1, b1mid)
    stripx0, stripw = LO - 2600, 1500
    for i in range(nF):
        y = yrow[i]
        c = BLUERED(bnorm(b1[i])) if not np.isnan(b1[i]) else "#cccccc"
        axH.add_patch(Rectangle((stripx0, y - 0.5), stripw, 1.0, facecolor=c,
                                edgecolor="none", zorder=2, clip_on=False))
    # ALT ticks
    pos_arr = V.pos.values
    cls_arr = V.cls.values
    for cls, (mk, base) in (("snp", (".", 1.2)), ("smallindel", ("s", 4)),
                            ("sv", ("D", 9))):
        cols = np.where(cls_arr == cls)[0]
        xs, ys = [], []
        for j in cols:
            car = np.where(G[:, j] == 1)[0]
            xs.extend([pos_arr[j]] * len(car))
            ys.extend([yrow[i] for i in car])
        if cls == "snp":
            axH.scatter(xs, ys, s=0.6, marker=".", c="#7a7a7a", alpha=0.5,
                        linewidths=0, zorder=3)
        else:
            col = "#E67E22" if cls == "smallindel" else "#C0392B"
            axH.scatter(xs, ys, s=base, marker=mk, c=col, alpha=0.85,
                        edgecolors="none", zorder=4)
    axH.set_ylim(-1, nF)
    axH.set_ylabel(f"founders (n={nF}), cold → hot home bio1", fontsize=9)
    axH.set_yticks([])
    for sp in ("top", "right", "left"):
        axH.spines[sp].set_visible(False)
    # genomic position ruler under the founder panel (shared x)
    ticks = np.arange(13_090_000, HI + 1, 5_000)
    axH.set_xticks(ticks)
    axH.set_xticklabels([f"{t:,}" for t in ticks], fontsize=7, rotation=0)
    axH.tick_params(labelbottom=True, length=3)
    axH.set_xlabel(f"{CHROM} position (bp)", fontsize=9)
    axH.text(stripx0 + stripw / 2, nF + 1, "home\nbio1", ha="center", va="bottom",
             fontsize=6.5, color="#555", clip_on=False)
    # bio1 strip colourbar
    smb = plt.cm.ScalarMappable(cmap=BLUERED, norm=bnorm)
    cbb = fig.colorbar(smb, ax=axH, fraction=0.02, pad=0.005)
    cbb.set_label("mean home bio1 (C)", fontsize=7.5)
    cbb.ax.tick_params(labelsize=6.5)

    # ================= Panel C: founder LD triangle =================
    n = len(ldpos)
    sx = (ldpos.max() - ldpos.min()) / (n - 1)
    ii, jj = np.triu_indices(n, k=1)
    xc = ldpos.min() + ((ii + jj) / 2) / (n - 1) * (ldpos.max() - ldpos.min())
    yc = -(jj - ii) * sx / 2
    r2 = R2[ii, jj]
    ok = ~np.isnan(r2)
    axL.scatter(xc[ok], yc[ok], c=r2[ok], cmap=LDCMAP, vmin=0, vmax=1,
                marker="s", s=3, edgecolors="none")
    # connector from lead SV straight down
    axL.axvline(LEAD_POS, ls="--", color="#b2182b", lw=0.7, alpha=0.8)
    axL.text(LO, 0, "founder LD (r$^2$, MAF≥0.05 + lead SV)", ha="left",
             va="top", fontsize=7, color="#555")
    axL.annotate("lead SV sits in tiling block Chr2_5306 but is merged into the CARK\n"
                 "block Chr2_5291 (sparse-SV merge → 40 kb WZA unit). Founder r²≤0.05\n"
                 "to the CARK LD island, r²=0.88 to its own local block (AT2G30800/810)",
                 xy=(LEAD_POS, yc.min() * 0.28), xytext=(LO + 1500, yc.min() * 0.60),
                 fontsize=7, color="#7d1a12", ha="left",
                 arrowprops=dict(arrowstyle="->", color="#7d1a12", lw=0.9))
    axL.set_ylim(yc.min() * 1.02, sx)
    axL.axis("off")
    smL = plt.cm.ScalarMappable(cmap=LDCMAP, norm=mcolors.Normalize(0, 1))
    cbl = fig.colorbar(smL, ax=axL, fraction=0.02, pad=0.005)
    cbl.set_label(r"$r^2$", fontsize=8)
    cbl.ax.tick_params(labelsize=6.5)

    axL.set_xlim(LO - 3200, HI + 200)                    # position ruler lives on axH

    out = f"{HERE}/cark_combined_manhattan_hap"
    fig.savefig(out + ".png", dpi=170, bbox_inches="tight")
    fig.savefig(out + ".pdf", bbox_inches="tight")
    print(f"wrote {out}.png / .pdf")

    # also dump the GEA variant table for the record
    gea.sort_values("nlp", ascending=False).to_csv(f"{HERE}/cark_gea_variants_{AXIS}.csv",
                                                    index=False)
    print(f"wrote cark_gea_variants_{AXIS}.csv")


if __name__ == "__main__":
    main()
