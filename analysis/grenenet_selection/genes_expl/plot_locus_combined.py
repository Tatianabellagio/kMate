#!/usr/bin/env python
"""Generalized cam5/CARK-style combined locus figure for one candidate gene.

Same 3 stacked panels sharing a genomic x-axis:
  A  LFMM GEA Manhattan (the candidate's axis) — wza_in cohort variants (SNP/indel/
     SV), -log10 p, colour = climate Δp (hot-cold, mean gen 1-3); per-class genome-
     wide Bonferroni lines; lead variant emphasised; gene-model track (candidate hi).
  B  Founder panel (231) ordered cold→hot by home bio1, origin strip, ALT ticks.
  C  Founder-LD triangle (r², MAF≥0.05 + aggregated lead variant) + connector.
Plus a printed LD-confirm: founder r² of the lead variant vs the candidate gene's
own variants (is the SV/indel genuinely tied to the gene, or CARK-style detached?).

Usage (one candidate; JSON config as argv[1]):
  PY plot_locus_combined.py '{"gene":"AT4G11600","sym":"GPX6","chrom":"Chr4",
     "gstart":7009769,"gend":7011375,"vpos":7011705,"ref_len":1225,"alt_len":61,
     "axis":"pc1","pad":8000}'
env: kmate.  Compute node.
"""
from __future__ import annotations
import os, sys, json, subprocess
import numpy as np
import pandas as pd
import pysam
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(HERE, "..")))
import lib

WZAIN = f"{lib.GEA}/phase1_replication/results/multiaxis/wza_in_clq09_tile"
GROUP_MEANS = f"{lib.GEA}/common/results/group_means.npz"
BIO_CSV = ("/global/scratch/users/tbellg/gea_grene-net/key_files/"
           "1001g_regmap_grenet_ecotype_info_corrected_bioclim_2024May16.csv")
BCF = "/global/home/users/tbellg/miniforge3/envs/kmate/bin/bcftools"
BONF = {"snp": 7.606, "smallindel": 7.131, "sv": 5.742}
CLSMARK = {"snp": ("o", 22), "smallindel": ("s", 26), "sv": ("D", 40)}
CLSNAME = {"snp": "SNP", "smallindel": "indel", "sv": "SV"}
DPCMAP = mcolors.LinearSegmentedColormap.from_list("dp", ["#2166ac", "#f7f7f7", "#b2182b"])
BLUERED = mcolors.LinearSegmentedColormap.from_list("br", ["#2166ac", "#b3b3b3", "#b2182b"])
LDCMAP = mcolors.LinearSegmentedColormap.from_list(
    "ld", ["#ffffff", "#e0e0f0", "#b3a2c8", "#8073ac", "#542788", "#2d004b"])


def vclass(rl, al):
    d = abs(al - rl)
    return "snp" if d == 0 else ("smallindel" if d < 50 else "sv")


def diverging(mid, vmax):
    return mcolors.TwoSlopeNorm(vmin=mid - vmax, vcenter=mid, vmax=mid + vmax)


def region_vcf(cfg, lo, hi):
    ch = cfg["chrom"]; ci = ch.replace("Chr", "")
    panel = f"{lib.PROJ}/panel/arch3/chr{ci}/merged_231_chr{ci}_final.vcf.gz"
    out = f"{HERE}/regions/{cfg['sym']}_{ch}_{lo}_{hi}.vcf.gz"
    os.makedirs(f"{HERE}/regions", exist_ok=True)
    if not os.path.exists(out):
        subprocess.run([BCF, "view", "-r", f"{ch}:{lo}-{hi}", "-e", "AC=0",
                        "-Oz", "-o", out, panel], check=True)
        subprocess.run([BCF, "index", "-t", out], check=True)
    return out


def load_gea(cfg, lo, hi):
    frames = []
    for cls in ("snp", "smallindel", "sv"):
        f = f"{WZAIN}/lfmm_{cls}_gen9_{cfg['axis']}.csv"
        if not os.path.exists(f):
            continue
        d = pd.read_csv(f)
        d = d[(d.chrom == cfg["chrom"]) & (d.pos >= lo) & (d.pos <= hi) & (d.MAF > 0.05)].copy()
        d["cls"] = cls
        frames.append(d)
    g = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
    if not len(g):
        return g
    g["nlp"] = -np.log10(g.pval.clip(lower=1e-300))
    z = np.load(GROUP_MEANS, allow_pickle=False)
    div = (z["hot_g1"] - z["cold_g1"]) + (z["hot_g2"] - z["cold_g2"]) + (z["hot_g3"] - z["cold_g3"])
    key = pd.DataFrame({"chrom": z["chrom"], "pos": z["pos"], "ref_len": z["ref_len"],
                        "alt_len": z["alt_len"], "dp": div / 3.0})
    key = key[(key.chrom == cfg["chrom"]) & (key.pos >= lo) & (key.pos <= hi)]
    return g.merge(key, on=["chrom", "pos", "ref_len", "alt_len"], how="left")


def load_panel(vcf, ch, lo, hi):
    v = pysam.VariantFile(vcf)
    samp = [int(s) for s in v.header.samples]
    recs, G = [], []
    for r in v.fetch(ch, lo, hi):
        rl, al = len(r.ref), len(r.alts[0])
        recs.append((r.pos, rl, al, vclass(rl, al)))
        G.append([1 if r.samples[s].get("GT", (None,))[0] == 1 else 0 for s in v.header.samples])
    G = np.array(G, np.int8).T
    V = pd.DataFrame(recs, columns=["pos", "ref_len", "alt_len", "cls"])
    bio1 = pd.read_csv(BIO_CSV).set_index("ecotypeid")["bio1"].to_dict()
    b1 = np.array([bio1.get(e, np.nan) for e in samp])
    return np.array(samp), b1, V, G


def agg_lead(V, G, cfg):
    """Aggregate ALT carriers of the lead variant at its position (multiallelic-safe)."""
    size = abs(cfg["alt_len"] - cfg["ref_len"])
    at = V.pos.values == cfg["vpos"]
    if size >= 50:
        cols = np.where(at & (np.abs(V.alt_len.values - V.ref_len.values) >= max(50, size // 2)))[0]
    else:
        cols = np.where(at & (V.ref_len.values == cfg["ref_len"]) &
                        (V.alt_len.values == cfg["alt_len"]))[0]
    if not len(cols):
        cols = np.where(at)[0]
    return (G[:, cols].sum(1) > 0).astype(float) if len(cols) else None


def r2(a, b):
    a = a - a.mean(); b = b - b.mean()
    if a.std() == 0 or b.std() == 0:
        return np.nan
    return float((np.mean(a * b) / (a.std() * b.std())) ** 2)


def ld_confirm(V, G, cfg):
    lead = agg_lead(V, G, cfg)
    if lead is None or lead.sum() == 0:
        return None
    gm = G.mean(0); gm = np.minimum(gm, 1 - gm)
    ingene = np.where((V.pos.values >= cfg["gstart"]) & (V.pos.values <= cfg["gend"]) & (gm >= 0.05))[0]
    r2g = [r2(lead, G[:, j].astype(float)) for j in ingene]
    lo, hi = cfg["_lo"], cfg["_hi"]
    loc = np.where((gm >= 0.05) & (V.pos.values != cfg["vpos"]))[0]
    r2l = [r2(lead, G[:, j].astype(float)) for j in loc]
    return dict(n_carr=int(lead.sum()),
                max_gene=float(np.nanmax(r2g)) if r2g else np.nan,
                med_gene=float(np.nanmedian(r2g)) if r2g else np.nan,
                n_gene=len(ingene),
                max_local=float(np.nanmax(r2l)) if r2l else np.nan)


def gene_track(ax, cfg, y0, dy, lo, hi):
    genes = lib.load_genes()
    gs = genes[(genes.chrom == cfg["chrom"]) & (genes.end >= lo) & (genes.start <= hi)]
    for _, g in gs.iterrows():
        hl = g.gene == cfg["gene"]
        col = "#b2182b" if hl else "#555"
        ax.add_patch(Rectangle((g.start, y0 - dy / 2), g.end - g.start, dy, facecolor=col,
                     edgecolor="none", alpha=0.95 if hl else 0.5, zorder=3, clip_on=False))
        lab = cfg["sym"] if hl else g.gene
        ax.text((g.start + g.end) / 2, y0 + dy * 0.9, lab, ha="center", va="bottom",
                fontsize=6.5 if hl else 4.8, style="italic",
                color=col, fontweight="bold" if hl else "normal", clip_on=False)


def main():
    cfg = json.loads(sys.argv[1])
    pad = cfg.get("pad", 8000)
    ch = cfg["chrom"]
    lo = min(cfg["gstart"], cfg["vpos"]) - pad
    hi = max(cfg["gend"], cfg["vpos"] + cfg["ref_len"]) + pad
    cfg["_lo"], cfg["_hi"] = lo, hi
    vcf = region_vcf(cfg, lo, hi)
    gea = load_gea(cfg, lo, hi)
    samp, b1, V, G = load_panel(vcf, ch, lo, hi)
    ld = ld_confirm(V, G, cfg)
    print(f"[{cfg['sym']}] window {ch}:{lo}-{hi} | GEA vars {len(gea)} | panel {G.shape}")
    if ld:
        print(f"  LD-confirm: lead carriers={ld['n_carr']} | r²(lead↔{cfg['sym']} gene, "
              f"n={ld['n_gene']}): max={ld['max_gene']:.3f} med={ld['med_gene']:.3f} | "
              f"max r² local={ld['max_local']:.3f}")

    dpmax = np.nanpercentile(np.abs(gea.dp.dropna()), 98) if gea.dp.notna().any() else 0.05
    dnorm = diverging(0, dpmax or 0.05)
    b1mid = float(np.nanmean(b1))
    bnorm = diverging(b1mid, np.nanmax(np.abs(b1 - b1mid)))
    order = np.argsort(np.nan_to_num(b1, nan=1e9)); nF = len(samp)
    yrow = {int(i): nF - 1 - k for k, i in enumerate(order)}
    lead_row = gea[(gea.pos == cfg["vpos"])]

    fig = plt.figure(figsize=(12, 15))
    gs = fig.add_gridspec(3, 1, height_ratios=[3.0, 6.0, 2.0], hspace=0.06)
    axM = fig.add_subplot(gs[0]); axH = fig.add_subplot(gs[1], sharex=axM)
    axL = fig.add_subplot(gs[2], sharex=axM)

    for ax in (axM, axH):
        ax.axvspan(cfg["gstart"], cfg["gend"], color="#cdd1d5", alpha=0.4, zorder=0, lw=0)
        ax.axvline(cfg["vpos"], ls="--", color="#b2182b", lw=0.6, alpha=0.6, zorder=1)

    # A: Manhattan
    for cls, sub in gea.groupby("cls"):
        mk, sz = CLSMARK[cls]
        axM.scatter(sub.pos, sub.nlp, c=sub.dp.fillna(0), cmap=DPCMAP, norm=dnorm,
                    marker=mk, s=sz, edgecolors="white", linewidths=0.3, zorder=3, label=CLSNAME[cls])
    for cls in ("sv", "snp"):
        axM.axhline(BONF[cls], ls="--", color="#888", lw=0.8)
        axM.text(hi, BONF[cls], f"{CLSNAME[cls]} Bonf", ha="right", va="bottom", fontsize=6.5, color="#666")
    if len(lead_row):
        lr = lead_row.sort_values("nlp").iloc[-1]
        axM.scatter([cfg["vpos"]], [lr.nlp], marker=CLSMARK[lr.cls][0], s=110,
                    facecolor=DPCMAP(dnorm(lr.dp if lr.dp == lr.dp else 0)),
                    edgecolor="#b2182b", linewidths=1.5, zorder=5)
        sz = abs(cfg["alt_len"] - cfg["ref_len"])
        axM.annotate(f"{cfg['sym']} lead {CLSNAME[lr.cls]}\n({sz} bp, {cfg['region']}, {cfg['axis']})",
                     xy=(cfg["vpos"], lr.nlp), xytext=(cfg["vpos"], lr.nlp + 1.2),
                     fontsize=7, color="#7d1a12", ha="center",
                     arrowprops=dict(arrowstyle="->", color="#7d1a12", lw=0.8))
    ytop = max(9.0, (gea.nlp.max() if len(gea) else 8) + 1.8)
    axM.set_ylim(-0.4, ytop); axM.set_ylabel(rf"$-\log_{{10}}p$ LFMM ({cfg['axis']})", fontsize=9)
    gene_track(axM, cfg, ytop - 0.5, 0.34, lo, hi)
    axM.legend(handles=[Line2D([], [], marker=CLSMARK[c][0], ls="", mfc="grey", mec="white",
               ms=6, label=CLSNAME[c]) for c in ("snp", "smallindel", "sv")],
               loc="upper left", frameon=False, fontsize=7)
    for sp in ("top", "right"):
        axM.spines[sp].set_visible(False)
    axM.tick_params(labelbottom=False)
    cb = fig.colorbar(plt.cm.ScalarMappable(cmap=DPCMAP, norm=dnorm), ax=axM, fraction=0.02, pad=0.005)
    cb.set_label(r"$\Delta p$ (hot$-$cold, gen1-3)", fontsize=7); cb.ax.tick_params(labelsize=6.5)

    # B: founder panel
    stripx0, stripw = lo - (hi - lo) * 0.045, (hi - lo) * 0.025
    for i in range(nF):
        c = BLUERED(bnorm(b1[i])) if not np.isnan(b1[i]) else "#ccc"
        axH.add_patch(Rectangle((stripx0, yrow[i] - 0.5), stripw, 1.0, facecolor=c,
                      edgecolor="none", zorder=2, clip_on=False))
    for cls, (mk, base, col) in (("snp", (".", 0.6, "#7a7a7a")), ("smallindel", ("s", 4, "#E67E22")),
                                 ("sv", ("D", 9, "#C0392B"))):
        cols = np.where(V.cls.values == cls)[0]
        xs, ys = [], []
        for j in cols:
            car = np.where(G[:, j] == 1)[0]
            xs.extend([V.pos.values[j]] * len(car)); ys.extend([yrow[i] for i in car])
        axH.scatter(xs, ys, s=base, marker=mk, c=col, alpha=0.5 if cls == "snp" else 0.85,
                    edgecolors="none", zorder=3 if cls == "snp" else 4)
    axH.set_ylim(-1, nF); axH.set_yticks([])
    axH.set_ylabel(f"founders (n={nF}), cold → hot home bio1", fontsize=9)
    ticks = np.linspace(lo, hi, 6).astype(int)
    axH.set_xticks(ticks); axH.set_xticklabels([f"{t:,}" for t in ticks], fontsize=7)
    axH.set_xlabel(f"{ch} position (bp)", fontsize=9)
    for sp in ("top", "right", "left"):
        axH.spines[sp].set_visible(False)
    cbb = fig.colorbar(plt.cm.ScalarMappable(cmap=BLUERED, norm=bnorm), ax=axH, fraction=0.02, pad=0.005)
    cbb.set_label("mean home bio1 (C)", fontsize=7); cbb.ax.tick_params(labelsize=6.5)

    # C: LD triangle
    gm = G.mean(0); gm = np.minimum(gm, 1 - gm)
    keep = np.where((gm >= 0.05))[0]
    pos = V.pos.values[keep].astype(float); M = G[:, keep].astype(float)
    o = np.argsort(pos); pos, M = pos[o], M[:, o]
    n = len(pos)
    if n > 2:
        Mc = M - M.mean(0); sd = Mc.std(0); sd[sd == 0] = np.nan
        R2 = (Mc.T @ Mc) / (len(M) * np.outer(sd, sd)); R2 = R2 ** 2
        sx = (pos.max() - pos.min()) / (n - 1)
        ii, jj = np.triu_indices(n, 1)
        xc = pos.min() + ((ii + jj) / 2) / (n - 1) * (pos.max() - pos.min())
        yc = -(jj - ii) * sx / 2; rr = R2[ii, jj]; ok = ~np.isnan(rr)
        axL.scatter(xc[ok], yc[ok], c=rr[ok], cmap=LDCMAP, vmin=0, vmax=1, marker="s", s=3, edgecolors="none")
        axL.axvline(cfg["vpos"], ls="--", color="#b2182b", lw=0.6, alpha=0.7)
        if ld:
            axL.text(lo, 0, f"founder LD (r², MAF≥0.05)  |  r²(lead↔{cfg['sym']} gene) "
                     f"max={ld['max_gene']:.2f}", ha="left", va="top", fontsize=7, color="#555")
        axL.set_ylim(yc.min() * 1.02, sx)
        fig.colorbar(plt.cm.ScalarMappable(cmap=LDCMAP, norm=mcolors.Normalize(0, 1)),
                     ax=axL, fraction=0.02, pad=0.005).set_label(r"$r^2$", fontsize=8)
    axL.set_xlim(lo - (hi - lo) * 0.06, hi); axL.axis("off")

    out = f"{HERE}/loci/{cfg['sym']}_combined"
    os.makedirs(f"{HERE}/loci", exist_ok=True)
    fig.savefig(out + ".png", dpi=165, bbox_inches="tight")
    fig.savefig(out + ".pdf", bbox_inches="tight")
    print(f"  wrote {out}.png/.pdf")


if __name__ == "__main__":
    main()
