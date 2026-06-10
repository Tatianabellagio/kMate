#!/usr/bin/env python
"""Micro-pangenome of AT3G50830 (COR413-PM2) + the intronic insertion bubble.

Top:  TAIR10 gene model (exons/CDS/UTR, strand) with the SV locus marked.
Bottom: the Minigraph-Cactus bubble at Chr3:18,895,050 as a sequence graph — the
        reference (no-insertion) path + the alternate insertion alleles, each edge
        width ∝ number of founders taking that path (from the 231-founder panel).

  PY=/global/home/users/tbellg/miniforge3/envs/basic/bin/python   # matplotlib env
  $PY analysis/grenenet_gea/plot_micropangenome.py
(founder genotypes are precomputed into the script via bcftools beforehand.)
"""
import os, subprocess, sys
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Rectangle
matplotlib.rcParams["pdf.fonttype"] = 42; matplotlib.rcParams["ps.fonttype"] = 42

GFF = "/global/home/users/tbellg/ara_key_files/TAIR10_GFF3_genes_transposons.gff"
VCF = "/global/scratch/users/tbellg/kmate/panel/arch3/chr3/merged_231_chr3_final.vcf.gz"
BCF = "/global/home/users/tbellg/miniforge3/envs/kmate/bin/bcftools"
FIGDIR = "/global/scratch/users/tbellg/kmate/results/grenenet_gea/gea/figures"
GENE = "AT3G50830"; ISO = "AT3G50830.1"; CHROM = "Chr3"; SV_POS = 18895050


def gene_model():
    feats = {"CDS": [], "exon": [], "five_prime_UTR": [], "three_prime_UTR": []}
    span = None; strand = "-"
    for ln in open(GFF):
        if GENE not in ln:
            continue
        c = ln.rstrip("\n").split("\t")
        if len(c) < 9:
            continue
        ft, s, e, st, attr = c[2], int(c[3]), int(c[4]), c[6], c[8]
        if ft == "gene" and f"ID={GENE}" in attr:
            span = (s, e); strand = st
        elif ft in feats and ISO in attr:
            feats[ft].append((s, e))
    return span, strand, feats


def founder_paths():
    """Per-founder allele at the bubble -> counts per path (ref / each insertion)."""
    out = subprocess.run([BCF, "query", "-r", f"{CHROM}:{SV_POS}",
                          "-f", "%ID\t%REF\t%ALT[\t%GT]\n", VCF],
                         capture_output=True, text=True).stdout.strip("\n").split("\n")
    recs = []
    for ln in out:
        f = ln.split("\t"); rid, ref, alt = f[0], f[1], f[2]
        gts = f[3:]
        ins_len = len(alt) - len(ref)
        recs.append(dict(id=rid, ins_len=ins_len, gts=gts))
    n = len(recs[0]["gts"])
    counts = {}; nmiss = 0; nref = 0
    for i in range(n):
        carried = [r for r in recs if "1" in r["gts"][i]]
        allmiss = all("." in r["gts"][i] and not any(ch.isdigit() for ch in r["gts"][i]) for r in recs)
        if carried:
            key = f"+{carried[0]['ins_len']} bp"
            counts[key] = counts.get(key, 0) + 1
        elif allmiss:
            nmiss += 1
        else:
            nref += 1
    return n, nref, nmiss, counts, recs


def main():
    span, strand, feats = gene_model()
    gs, ge = span
    n, nref, nmiss, counts, recs = founder_paths()
    ncall = n - nmiss
    # order insertion alleles by founder count (desc)
    ins = sorted(counts.items(), key=lambda kv: -kv[1])
    print(f"founders={n}, called={ncall}, reference(no-ins)={nref}, missing={nmiss}, "
          f"insertion alleles={ins}")

    fig = plt.figure(figsize=(8.2, 5.2))
    gsx = fig.add_gridspec(2, 1, height_ratios=[1, 1.35], hspace=0.55)
    axg = fig.add_subplot(gsx[0]); axb = fig.add_subplot(gsx[1])

    # ---- gene model ----
    pad = 250
    axg.plot([gs, ge], [0.5, 0.5], "-", color="0.4", lw=1.2, zorder=1)   # intron backbone
    for s, e in feats["exon"]:                                            # UTR-height exon
        axg.add_patch(Rectangle((s, 0.40), e - s, 0.20, fc="#cfd8e3", ec="k", lw=.5, zorder=2))
    for s, e in feats["CDS"]:                                            # tall CDS
        axg.add_patch(Rectangle((s, 0.32), e - s, 0.36, fc="#5b7aa6", ec="k", lw=.5, zorder=3))
    # strand arrow (minus = 5'->3' is right->left)
    x0, x1 = (ge, gs) if strand == "-" else (gs, ge)
    axg.annotate("", xy=(x1, 0.85), xytext=(x0, 0.85),
                 arrowprops=dict(arrowstyle="->", color="0.5", lw=1))
    axg.text((gs + ge) / 2, 0.95, f"{GENE}  (COR413-PM2)  ·  {strand} strand", ha="center", va="bottom", fontsize=9)
    # SV marker in the intron
    axg.plot([SV_POS], [0.5], "v", color="crimson", ms=9, zorder=5)
    axg.annotate("insertion bubble\n(intron 1)", xy=(SV_POS, 0.5), xytext=(SV_POS, 0.05),
                 ha="center", va="bottom", fontsize=7.5, color="crimson",
                 arrowprops=dict(arrowstyle="-", color="crimson", lw=.6))
    axg.set_xlim(gs - pad, ge + pad); axg.set_ylim(0, 1.1)
    axg.set_yticks([]); axg.set_xlabel(f"{CHROM} position (bp)", fontsize=8)
    axg.tick_params(labelsize=7)
    for sp in ("top", "left", "right"):
        axg.spines[sp].set_visible(False)
    axg.set_title("Gene model", fontsize=9, loc="left")

    # ---- pangenome bubble ----
    axb.set_xlim(0, 1); axb.set_ylim(0, 1); axb.axis("off")
    axb.set_title(f"Pangenome bubble at {CHROM}:{SV_POS:,}  —  {n} founders", fontsize=9, loc="left")
    lx, rx = 0.14, 0.86
    for x in (lx, rx):                                                   # flanking ref anchors
        axb.add_patch(Rectangle((x - 0.035, 0.45), 0.07, 0.10, fc="0.6", ec="k", lw=.6))
    axb.text(lx, 0.36, "intron", ha="center", fontsize=7, color="0.4")
    axb.text(rx, 0.36, "intron", ha="center", fontsize=7, color="0.4")

    def wlin(cnt):                                                       # edge width ∝ founders
        return 1.2 + 6.5 * cnt / max(ncall, 1)
    # reference path (straight)
    axb.plot([lx + .035, rx - .035], [0.5, 0.5], "-", color="0.45", lw=wlin(nref),
             solid_capstyle="round", zorder=2)
    axb.text((lx + rx) / 2, 0.5 - 0.075, f"reference — no insertion   ({nref} founders)",
             ha="center", va="top", fontsize=7.5, color="0.3")
    # insertion alleles as arcs above, width ∝ count; highlight the lead (181 bp)
    ys = [0.74, 0.9, 1.02]
    for k, (lab, cnt) in enumerate(ins):
        y = ys[min(k, len(ys) - 1)]
        lead = (cnt == max(counts.values()))
        col = "crimson" if lead else "#8a8a8a"
        arc = FancyArrowPatch((lx + .035, 0.55), (rx - .035, 0.55),
                              connectionstyle=f"arc3,rad={-0.35 - 0.12 * k}",
                              arrowstyle="-", lw=wlin(cnt), color=col, alpha=.9, zorder=3)
        axb.add_patch(arc)
        af = cnt / ncall
        tag = f"{lab} insertion — {cnt} founders (AF {af:.2f})" + ("   ◄ candidate SV" if lead else "")
        axb.text((lx + rx) / 2, y + 0.02, tag, ha="center", va="bottom",
                 fontsize=8 if lead else 7, color=col, fontweight="bold" if lead else "normal")
    axb.text(0.5, 0.06,
             f"{ncall} founders genotyped · edge width ∝ #founders on each path",
             ha="center", fontsize=7, color="0.45")

    fig.suptitle("Micro-pangenome: COR413-PM2 and its intronic structural-variant bubble",
                 fontsize=11, y=0.99)
    os.makedirs(FIGDIR, exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(f"{FIGDIR}/micropangenome_COR413.{ext}", bbox_inches="tight", dpi=150)
    print(f"wrote {FIGDIR}/micropangenome_COR413.pdf (+.png)")


if __name__ == "__main__":
    main()
