#!/usr/bin/env python
"""Generalized cam5/CARK-style combined locus figure for one candidate gene.

Same 3 stacked panels sharing a genomic x-axis:
  A  Manhattan of the upstream scan — see `scan` below.
  B  Founder panel (231) ordered cold→hot by home bio1, origin strip, ALT ticks.
  C  Founder-LD triangle (r², MAF≥0.05 + aggregated lead variant) + connector.
Plus a printed LD-confirm: founder r² of the lead variant vs the candidate gene's
own variants (is the SV/indel genuinely tied to the gene, or CARK-style detached?).

`scan` selects which upstream scan panel A shows. A gene found only by one scan must
be plotted against THAT scan — a GWAS-only candidate drawn over the LFMM GEA panel
shows an empty Manhattan and reads as a null (this is what CML50 looked like before
the gwas mode existed).

  "gea"  (default)  LFMM per-record p on the candidate's climate axis, from
                    `wza_in_clq09_tile`; colour = climate Δp (hot−cold, mean gen 1-3);
                    Bonferroni lines per variant class.
  "gwas"            GEMMA per-garden p on the founder selection coefficient `s`, from
                    `r3_persite_gwas/results/gemma_gwas/persite_gwas_{snp,nonsnp,sv}.npz`,
                    for ONE garden (`garden`, default = the lead's median-strength
                    significant garden); colour = signed GEMMA Z (effect direction on
                    `s`). Adds panel D: the lead variant's p across all 30 gardens
                    ordered by garden bio1, so the recurrence count can be read
                    against the thresholds rather than taken on trust.

  ⚠ In gwas mode THREE Bonferroni lines are drawn (snp 7.55, nonsnp 7.02, sv 5.41)
  because the per-class scans carry very different marker counts (1.75M / 525k /
  12.8k). The SAME marker and the SAME p can clear the sv line and miss the nonsnp
  line, so a garden-recurrence count is only meaningful with its class named.

Usage (one candidate; JSON config as argv[1]):
  PY plot_locus_combined.py '{"gene":"AT4G11600","sym":"GPX6","chrom":"Chr4",
     "gstart":7009769,"gend":7011375,"vpos":7011705,"ref_len":1225,"alt_len":61,
     "axis":"pc1","pad":8000}'
  PY plot_locus_combined.py '{"gene":"AT5G04170","sym":"CML50","chrom":"Chr5",
     "gstart":1145431,"gend":1147781,"vpos":1146373,"ref_len":95,"alt_len":1,
     "scan":"gwas","pad":8000}'
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
sys.path.insert(0, os.path.abspath(os.path.join(HERE, "..", "..")))
import lib
sys.path.insert(0, os.path.join(HERE, "..", "convergence"))
import plot_theme as TH                                   # noqa: E402

TH.apply()

WZAIN = f"{lib.GEA}/r2_gea_nonsnp/phase1_replication/results/multiaxis/wza_in_clq09_tile"
GROUP_MEANS = f"{lib.GEA}/common/results/group_means.npz"
BIO_CSV = ("/global/scratch/users/tbellg/gea_grene-net/key_files/"
           "1001g_regmap_grenet_ecotype_info_corrected_bioclim_2024May16.csv")
BCF = "/global/home/users/tbellg/miniforge3/envs/kmate/bin/bcftools"
GWASD = f"{lib.GEA}/r3_persite_gwas/results/gemma_gwas"
PANEL = "/global/scratch/projects/fc_moilab/tbellg/kmate/panel/arch3"
BONF = {"snp": 7.606, "smallindel": 7.131, "sv": 5.742}
# gwas mode: per-class genome-wide Bonferroni from each scan's OWN marker count
# (snp 1,752,846 / nonsnp 525,043 / sv 12,789) -- they differ by ~140x, which is why
# one marker's "n gardens" depends on which class scan it is counted in.
BONF_GWAS = {"snp": 7.545, "nonsnp": 7.021, "sv": 5.408}
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


def _panel_sizes(ch, want):
    """chrom/pos -> |alt_len-ref_len|, LARGEST record at the position.

    The GEMMA result arrays carry only chrom/pos (no allele), so variant class for the
    Manhattan marker shape has to come back from the panel meta. Multiallelic positions
    take the largest record -- the same guess `build_gwas_pool.py` makes, and the same
    caveat applies (see memory `panel-multiallelic-pos-key-trap`). This only sets a
    marker glyph here; nothing quantitative depends on it.
    """
    f = f"{PANEL}/{ch.lower()}/var_pa_231_arch3_{ch.lower()}.meta.npz"
    d = np.load(f, allow_pickle=True)
    sel = np.isin(d["pos"], np.asarray(want))
    s = pd.DataFrame({"pos": d["pos"][sel],
                      "size": np.abs(d["alt_len"][sel].astype(int)
                                     - d["ref_len"][sel].astype(int))})
    return s.groupby("pos")["size"].max().to_dict()


def lead_gwas_profile(cfg):
    """Lead variant's -log10 p across all 30 gardens, per class scan.

    Returns (DataFrame[garden, bio1, nlp_<cls>...], sites) or None. This is what makes
    a recurrence count auditable: the p-values decay smoothly across gardens, so
    'n gardens' is wherever a class threshold happens to cut that gradient.
    """
    out, sites, bio1 = {}, None, None
    for cls in ("snp", "nonsnp", "sv"):
        f = f"{GWASD}/persite_gwas_{cls}.npz"
        if not os.path.exists(f):
            continue
        z = np.load(f, allow_pickle=True)
        i = np.where((z["chrom"] == cfg["chrom"].lower()) & (z["pos"] == cfg["vpos"]))[0]
        if not len(i):
            continue
        sites, bio1 = z["sites"], z["bio1"]
        out[cls] = -np.log10(np.clip(z["P"][int(i[0])], 1e-300, None))
    if not out:
        return None
    T = pd.DataFrame({"garden": sites, "bio1": bio1})
    for cls, v in out.items():
        T[cls] = v
    return T


def load_gwas(cfg, lo, hi, garden):
    """GEMMA per-garden p for every marker in the window, for ONE garden.

    snp + nonsnp are the two disjoint scans covering all markers; the sv scan is a
    SUBSET of nonsnp (re-thresholded on its own marker count), so plotting it too
    would draw SVs twice. Class for the marker glyph comes from the panel meta.
    """
    frames = []
    for cls in ("snp", "nonsnp"):
        f = f"{GWASD}/persite_gwas_{cls}.npz"
        if not os.path.exists(f):
            continue
        z = np.load(f, allow_pickle=True)
        gi = np.where(z["sites"] == garden)[0]
        if not len(gi):
            continue
        gi = int(gi[0])
        m = (z["chrom"] == cfg["chrom"].lower()) & (z["pos"] >= lo) & (z["pos"] <= hi)
        if not m.any():
            continue
        d = pd.DataFrame({"pos": z["pos"][m],
                          "nlp": -np.log10(np.clip(z["P"][m, gi], 1e-300, None)),
                          "Z": z["Z"][m, gi]})
        d["scan"] = cls
        frames.append(d)
    if not frames:
        return pd.DataFrame()
    g = pd.concat(frames, ignore_index=True)
    sz = _panel_sizes(cfg["chrom"], g.pos.unique())
    g["cls"] = [vclass(1, 1 + sz.get(int(p), 0)) for p in g.pos]
    # snp scan markers stay 'snp' even if the position also carries a longer record
    g.loc[g.scan == "snp", "cls"] = "snp"
    return g


def load_panel(vcf, ch, lo, hi):
    """Founder genotypes in the window, plus a NO-CALL mask.

    The mask is not cosmetic. This panel used to score `GT == 1` and call everything else
    0, which renders a no-call identically to a reference call. For a well-genotyped SNP
    that is harmless; for an SV it is a misstatement -- GPX6's 1.16 kb deletion is called
    in 34% of the 231 founders, so two thirds of the column was being drawn as REF when it
    is unknown (memory `sv-panel-support-asymmetry`). G stays 0/1 so the LD triangle and
    r2/D' are unchanged; MISS is carried alongside for display.
    """
    v = pysam.VariantFile(vcf)
    samp = [int(s) for s in v.header.samples]
    recs, G, MISS = [], [], []
    for r in v.fetch(ch, lo, hi):
        rl, al = len(r.ref), len(r.alts[0])
        recs.append((r.pos, rl, al, vclass(rl, al)))
        gt = [r.samples[s].get("GT", (None,))[0] for s in v.header.samples]
        G.append([1 if g == 1 else 0 for g in gt])
        MISS.append([g is None for g in gt])
    G = np.array(G, np.int8).T
    MISS = np.array(MISS, bool).T
    V = pd.DataFrame(recs, columns=["pos", "ref_len", "alt_len", "cls"])
    bio1 = pd.read_csv(BIO_CSV).set_index("ecotypeid")["bio1"].to_dict()
    b1 = np.array([bio1.get(e, np.nan) for e in samp])
    return np.array(samp), b1, V, G, MISS


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


def lead_cols(V, cfg):
    """Column indices of the focal variant (same rule as agg_lead, multiallelic-safe)."""
    size = abs(cfg["alt_len"] - cfg["ref_len"])
    at = V.pos.values == cfg["vpos"]
    if size >= 50:
        cols = np.where(at & (np.abs(V.alt_len.values - V.ref_len.values)
                              >= max(50, size // 2)))[0]
    else:
        cols = np.where(at & (V.ref_len.values == cfg["ref_len"]) &
                        (V.alt_len.values == cfg["alt_len"]))[0]
    return cols if len(cols) else np.where(at)[0]


def r2(a, b):
    a = a - a.mean(); b = b - b.mean()
    if a.std() == 0 or b.std() == 0:
        return np.nan
    return float((np.mean(a * b) / (a.std() * b.std())) ** 2)


def dprime(a, b):
    """|D'| — LD normalised by the max |D| the two allele frequencies permit.

    r² is bounded above by the frequency difference: a lead at 5% CANNOT exceed
    r²=0.06 against a 50% partner however perfectly nested it is. So r² alone
    systematically reads "unlinked" for every rare candidate, which is exactly the
    class this screen produces. D' is the frequency-robust companion. It has the
    opposite bias (it saturates at 1.0 for rare alleles on small samples), so the two
    are reported together and neither is read alone.
    """
    pA, pB = a.mean(), b.mean()
    D = float(np.mean(a * b) - pA * pB)
    dmax = (min(pA * (1 - pB), pB * (1 - pA)) if D > 0
            else min(pA * pB, (1 - pA) * (1 - pB)))
    return abs(D) / dmax if dmax > 0 else np.nan


def r2_ceiling(a, b):
    """Largest r² attainable between two variants at these allele frequencies."""
    p, q = sorted((float(a.mean()), float(b.mean())))
    return p * (1 - q) / (q * (1 - p)) if 0 < p and q < 1 else np.nan


def ld_confirm(V, G, cfg):
    lead = agg_lead(V, G, cfg)
    if lead is None or lead.sum() == 0:
        return None
    gm = G.mean(0); gm = np.minimum(gm, 1 - gm)
    # EXCLUDE the lead's own position. A lead inside the gene body (intron/UTR/CDS)
    # is otherwise compared against itself and max_gene reads a meaningless 1.000 --
    # `loc` below always excluded vpos, `ingene` did not, which is why max_gene could
    # come out ABOVE max_local. 40% of the candidate table has an in-gene lead.
    ingene = np.where((V.pos.values >= cfg["gstart"]) & (V.pos.values <= cfg["gend"])
                      & (gm >= 0.05) & (V.pos.values != cfg["vpos"]))[0]
    cols = [G[:, j].astype(float) for j in ingene]
    r2g = [r2(lead, c) for c in cols]
    dpg = [dprime(lead, c) for c in cols]
    ceil = [r2_ceiling(lead, c) for c in cols]
    lo, hi = cfg["_lo"], cfg["_hi"]
    loc = np.where((gm >= 0.05) & (V.pos.values != cfg["vpos"]))[0]
    r2l = [r2(lead, G[:, j].astype(float)) for j in loc]
    k = int(np.nanargmax(r2g)) if r2g and not np.all(np.isnan(r2g)) else None
    return dict(n_carr=int(lead.sum()), freq=float(lead.mean()),
                max_gene=float(np.nanmax(r2g)) if r2g else np.nan,
                med_gene=float(np.nanmedian(r2g)) if r2g else np.nan,
                n_gene=len(ingene),
                max_local=float(np.nanmax(r2l)) if r2l else np.nan,
                # how much of the achievable r² the best in-gene partner actually gets
                frac_ceiling=float(r2g[k] / ceil[k]) if k is not None and ceil[k] else np.nan,
                max_dprime=float(np.nanmax(dpg)) if dpg else np.nan,
                med_dprime=float(np.nanmedian(dpg)) if dpg else np.nan)


def load_gene_models(chrom, lo, hi):
    """Exon/CDS/UTR structure per gene in a window, from the transposons GFF.

    `lib.load_genes()` reads TAIR10_GFF3_genes.gff, which contains ONLY `gene` rows --
    no exons -- so a gene track built from it can only be a solid block and cannot say
    whether a variant is exonic or intronic. The sibling
    TAIR10_GFF3_genes_transposons.gff carries mRNA/exon/CDS/five_prime_UTR/
    three_prime_UTR, which is what JBrowse renders. Coordinates are the same TAIR10
    assembly, so this is purely additive detail on the same models.

    Returns {gene: {strand, start, end, exons[], cds[], n_iso}} for the LONGEST isoform
    of each gene (Arabidopsis genes here have 1-2; drawing every isoform makes a
    10-gene window unreadable, and n_iso records when one was collapsed).
    """
    tx, want = {}, chrom
    with open(lib.TAIR10_GENES_TE) as fh:
        for line in fh:
            if line[0] == "#":
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[0] != want:
                continue
            kind, s, e = f[2], int(f[3]), int(f[4])
            if e < lo or s > hi or kind not in ("mRNA", "exon", "CDS"):
                continue
            at = dict(kv.split("=", 1) for kv in f[8].rstrip(";").split(";") if "=" in kv)
            if kind == "mRNA":
                tid = at.get("ID", "")
                tx.setdefault(tid, dict(gene=at.get("Parent", ""), strand=f[6],
                                        start=s, end=e, exons=[], cds=[]))
            else:
                # a CDS Parent can list several ids ("AT5G04170.1,AT5G04170.1-Protein")
                for tid in at.get("Parent", "").split(","):
                    if tid in tx:
                        tx[tid]["exons" if kind == "exon" else "cds"].append((s, e))
    out = {}
    for tid, t in tx.items():
        g = t["gene"]
        prev = out.get(g)
        if prev is None or (t["end"] - t["start"]) > (prev["end"] - prev["start"]):
            out[g] = dict(t, n_iso=(prev or {}).get("n_iso", 0) + 1)
        else:
            prev["n_iso"] += 1
    for g in out:
        out[g]["exons"].sort()
        out[g]["cds"].sort()
    return out


def feature_at(model, start, end):
    """Which sub-genic feature a [start, end] interval hits, given one gene model."""
    def ov(iv):
        return [(a, b) for a, b in iv if a <= end and b >= start]
    if not model:
        return "no model"
    c, x = ov(model["cds"]), ov(model["exons"])
    if c:
        return "CDS"
    if x:
        return "exon (UTR)"
    ex = model["exons"]
    for (a1, b1), (a2, _) in zip(ex, ex[1:]):
        if start > b1 and end < a2:
            n = ex.index((a1, b1)) + 1
            n = n if model["strand"] == "+" else len(ex) - n
            return (f"intron {n} of {len(ex)-1} · {b1+1:,}-{a2-1:,} ({a2-b1-1:,} bp) · "
                    f"{start-b1-1:,} bp from the 5' donor")
    return "outside the transcript"


def gene_zoom_inset(ax, cfg, models, dstart, dend, featc):
    """Inset: the candidate gene alone, at gene scale, with the lead's footprint.

    In an 18 kb locus window a 2 kb gene is ~180 px, so exon/intron structure is
    visible but you cannot actually read whether the variant sits in an exon or an
    intron -- which is the first question asked of any coding-region candidate. This
    inset re-renders just the candidate transcript across the full inset width.
    """
    m = models.get(cfg["gene"])
    if m is None:
        return
    span = m["end"] - m["start"]
    gl, gr = m["start"] - span * 0.04, m["end"] + span * 0.04
    ins = ax.inset_axes([0.035, 0.46, 0.40, 0.30])
    ins.set_xlim(gl, gr); ins.set_ylim(-1.35, 1.5)
    col = "#b2182b"

    ins.plot([m["start"], m["end"]], [0, 0], color="#555", lw=1.0, zorder=2)
    step = max(span / 55, 1)
    mark = ">" if m["strand"] == "+" else "<"
    xs = np.arange(m["start"] + step / 2, m["end"], step)
    ins.plot(xs, np.zeros(len(xs)), marker=mark, ls="", ms=3.0, color="#555", zorder=3)
    for a, b in m["exons"]:
        ins.add_patch(Rectangle((a, -0.34), max(b - a, 1), 0.68, facecolor="#8a8a8a",
                      edgecolor="none", zorder=4))
    for a, b in m["cds"]:
        ins.add_patch(Rectangle((a, -0.62), max(b - a, 1), 1.24, facecolor="#4d4d4d",
                      edgecolor="none", zorder=5))
    ex = m["exons"] if m["strand"] == "+" else m["exons"][::-1]
    for k, (a, b) in enumerate(ex, 1):
        ins.text((a + b) / 2, 0.78, f"E{k}", ha="center", va="bottom", fontsize=5.5,
                 color="#444")
    # the lead's actual footprint, drawn to scale
    ins.axvspan(dstart, dend, color=col, alpha=0.30, lw=0, zorder=6)
    ins.plot([(dstart + dend) / 2], [-0.95], marker="^", ms=5, color=col, zorder=7)
    ins.text((dstart + dend) / 2, -1.32, f"{abs(cfg['alt_len']-cfg['ref_len'])} bp",
             ha="center", va="bottom", fontsize=5.8, color=col)


    ins.set_yticks([])
    ins.set_xticks([m["start"], m["end"]])
    ins.set_xticklabels([f"{m['start']:,}", f"{m['end']:,}"], fontsize=5.5)
    ins.tick_params(axis="x", length=2, pad=1)
    for sp in ("top", "right", "left"):
        ins.spines[sp].set_visible(False)
    ins.spines["bottom"].set_linewidth(0.5)
    # fully opaque: the parent's Bonferroni lines run through this y-range and
    # a semi-transparent inset lets them read as part of the gene model
    ins.patch.set_facecolor("white"); ins.patch.set_alpha(1.0)
    ins.add_patch(Rectangle((0, 0), 1, 1, transform=ins.transAxes, fill=False,
                            ec="#ccc", lw=0.6, zorder=1, clip_on=False))


def gene_track(ax, cfg, y0, dy, lo, hi, models):
    """JBrowse-style gene models: CDS thick, UTR thin, introns a line with chevrons."""
    genes = lib.load_genes()
    gs = genes[(genes.chrom == cfg["chrom"]) & (genes.end >= lo) & (genes.start <= hi)]
    for _, g in gs.iterrows():
        hl = g.gene == cfg["gene"]
        col = "#b2182b" if hl else "#666"
        al = 0.95 if hl else 0.55
        m = models.get(g.gene)
        if m is None:                                   # no transcript (e.g. a TE gene)
            ax.add_patch(Rectangle((g.start, y0 - dy * 0.25), g.end - g.start, dy * 0.5,
                         facecolor=col, edgecolor="none", alpha=al * 0.6, zorder=3,
                         clip_on=False))
        else:
            # intron backbone + strand chevrons
            ax.plot([m["start"], m["end"]], [y0, y0], color=col, lw=0.8, alpha=al,
                    zorder=3, clip_on=False, solid_capstyle="butt")
            step = max((hi - lo) / 90, 1)
            mark = ">" if m["strand"] == "+" else "<"
            xs = np.arange(m["start"] + step / 2, m["end"], step)
            if len(xs):
                ax.plot(xs, np.full(len(xs), y0), marker=mark, ls="", ms=2.3,
                        color=col, alpha=al * 0.85, zorder=4, clip_on=False)
            cds = m["cds"]
            for a, b in m["exons"]:                     # exon body = UTR height
                ax.add_patch(Rectangle((a, y0 - dy * 0.26), max(b - a, 1), dy * 0.52,
                             facecolor=col, edgecolor="none", alpha=al * 0.75,
                             zorder=5, clip_on=False))
            for a, b in cds:                            # CDS overlay = full height
                ax.add_patch(Rectangle((a, y0 - dy * 0.5), max(b - a, 1), dy,
                             facecolor=col, edgecolor="none", alpha=al, zorder=6,
                             clip_on=False))
        lab = cfg["sym"] if hl else g.gene
        ax.text((g.start + g.end) / 2, y0 + dy * 0.9, lab, ha="center", va="bottom",
                fontsize=6.5 if hl else 4.8, style="italic",
                color=col, clip_on=False)
    return models


def main():
    cfg = json.loads(sys.argv[1])
    pad = cfg.get("pad", 8000)
    ch = cfg["chrom"]
    lo = min(cfg["gstart"], cfg["vpos"]) - pad
    hi = max(cfg["gend"], cfg["vpos"] + cfg["ref_len"]) + pad
    cfg["_lo"], cfg["_hi"] = lo, hi
    scan = cfg.get("scan", "gea")
    vcf = region_vcf(cfg, lo, hi)
    prof = lead_gwas_profile(cfg) if scan == "gwas" else None
    if scan == "gwas":
        if prof is None:
            raise SystemExit(f"gwas mode: {ch}:{cfg['vpos']} is not in the GEMMA arrays")
        # default garden = the MEDIAN-strength significant garden in the lead's best
        # class scan, not the strongest. The strongest is the best case, not a typical
        # one, and picking it would overstate what a garden looks like.
        if "garden" in cfg:
            garden = int(cfg["garden"])
        else:
            cls = max((c for c in ("sv", "nonsnp", "snp") if c in prof),
                      key=lambda c: int((prof[c] >= BONF_GWAS[c]).sum()))
            s = prof[prof[cls] >= BONF_GWAS[cls]].sort_values(cls)
            s = s if len(s) else prof.sort_values(cls)
            garden = int(s.iloc[len(s) // 2].garden)
        cfg["_garden"] = garden
        gea = load_gwas(cfg, lo, hi, garden)
    else:
        gea = load_gea(cfg, lo, hi)
    samp, b1, V, G, MISS = load_panel(vcf, ch, lo, hi)
    # exon/CDS structure for the window, parsed once and reused by the gene track
    models = load_gene_models(ch, lo, hi)
    # the lead's footprint: for a deletion the REF spans pos..pos+ref_len-1 and the
    # first base is the VCF anchor, so the bases actually removed start at pos+1
    dstart = cfg["vpos"] + (1 if cfg["ref_len"] > cfg["alt_len"] else 0)
    dend = cfg["vpos"] + max(cfg["ref_len"] - 1, 0)
    featc = feature_at(models.get(cfg["gene"]), dstart, dend)
    print(f"  lead footprint {ch}:{dstart:,}-{dend:,} "
          f"({abs(cfg['alt_len']-cfg['ref_len'])} bp) lands in: {featc}")
    ld = ld_confirm(V, G, cfg)
    print(f"[{cfg['sym']}] window {ch}:{lo}-{hi} | GEA vars {len(gea)} | panel {G.shape}")
    if ld:
        print(f"  LD-confirm: lead carriers={ld['n_carr']} ({ld['freq']*100:.1f}%) | "
              f"r²(lead↔{cfg['sym']} gene, n={ld['n_gene']}): max={ld['max_gene']:.3f} "
              f"med={ld['med_gene']:.3f} | max r² local={ld['max_local']:.3f}")
        print(f"              frequency-aware: best r² is {ld['frac_ceiling']*100:.0f}% of "
              f"its ceiling | D' max={ld['max_dprime']:.3f} med={ld['med_dprime']:.3f}")

    if scan == "gwas":
        cvals = gea.Z if len(gea) else pd.Series([0.0])
        dnorm = diverging(0, float(np.nanpercentile(np.abs(cvals), 98)) or 1.0)
        cblab = r"GEMMA $Z$ on $s$ (garden %d)" % cfg["_garden"]
    else:
        dpmax = np.nanpercentile(np.abs(gea.dp.dropna()), 98) if gea.dp.notna().any() else 0.05
        dnorm = diverging(0, dpmax or 0.05)
        cblab = r"$\Delta p$ (hot$-$cold, gen1-3)"
    b1mid = float(np.nanmean(b1))
    bnorm = diverging(b1mid, np.nanmax(np.abs(b1 - b1mid)))
    order = np.argsort(np.nan_to_num(b1, nan=1e9)); nF = len(samp)
    yrow = {int(i): nF - 1 - k for k, i in enumerate(order)}
    lead_row = gea[(gea.pos == cfg["vpos"])]

    if scan == "gwas":
        fig = plt.figure(figsize=(12, 17.5))
        gs = fig.add_gridspec(4, 1, height_ratios=[3.0, 6.0, 2.0, 2.2], hspace=0.06)
        axG = fig.add_subplot(gs[3])          # own x-axis (gardens, not bp)
    else:
        fig = plt.figure(figsize=(12, 15))
        gs = fig.add_gridspec(3, 1, height_ratios=[3.0, 6.0, 2.0], hspace=0.06)
        axG = None
    axM = fig.add_subplot(gs[0]); axH = fig.add_subplot(gs[1], sharex=axM)
    axL = fig.add_subplot(gs[2], sharex=axM)

    for ax in (axM, axH):
        ax.axvspan(cfg["gstart"], cfg["gend"], color="#cdd1d5", alpha=0.4, zorder=0, lw=0)
        ax.axvline(cfg["vpos"], ls="--", color="#b2182b", lw=0.6, alpha=0.6, zorder=1)

    # A: Manhattan
    cvar = "Z" if scan == "gwas" else "dp"
    for cls, sub in gea.groupby("cls"):
        mk, sz = CLSMARK[cls]
        axM.scatter(sub.pos, sub.nlp, c=sub[cvar].fillna(0), cmap=DPCMAP, norm=dnorm,
                    marker=mk, s=sz, edgecolors="white", linewidths=0.3, zorder=3, label=CLSNAME[cls])
    bonf = BONF_GWAS if scan == "gwas" else BONF
    # only the SV-scan threshold is drawn: it is the one this candidate is called on
    lines = ("sv",) if scan == "gwas" else ("sv", "snp")
    for cls in lines:
        axM.axhline(bonf[cls], ls="--", color="#888", lw=0.8)
        axM.text(hi, bonf[cls], f"{CLSNAME[cls]} Bonf", ha="right", va="bottom",
                 fontsize=6.5, color="#666")
    if len(lead_row):
        lr = lead_row.sort_values("nlp").iloc[-1]
        axM.scatter([cfg["vpos"]], [lr.nlp], marker=CLSMARK[lr.cls][0], s=110,
                    facecolor=DPCMAP(dnorm(lr[cvar] if lr[cvar] == lr[cvar] else 0)),
                    edgecolor="#b2182b", linewidths=1.5, zorder=5)
    ytop = max(9.0, (gea.nlp.max() if len(gea) else 8) + 1.8)
    axM.set_ylim(-0.4, ytop)
    axM.set_ylabel(rf"$-\log_{{10}}p$ GEMMA on $s$ (garden {cfg['_garden']})"
                   if scan == "gwas" else
                   rf"$-\log_{{10}}p$ LFMM ({cfg['axis']})", fontsize=9)
    gene_track(axM, cfg, ytop - 0.5, 0.34, lo, hi, models)
    gene_zoom_inset(axM, cfg, models, dstart, dend, featc)
    axM.legend(handles=[Line2D([], [], marker=CLSMARK[c][0], ls="", mfc="grey", mec="white",
               ms=6, label=CLSNAME[c]) for c in ("snp", "smallindel", "sv")],
               loc="upper left", frameon=False, fontsize=7)
    for sp in ("top", "right"):
        axM.spines[sp].set_visible(False)
    axM.tick_params(labelbottom=False)
    cb = fig.colorbar(plt.cm.ScalarMappable(cmap=DPCMAP, norm=dnorm), ax=axM, fraction=0.02, pad=0.005)
    cb.set_label(cblab, fontsize=7); cb.ax.tick_params(labelsize=6.5)

    # B: founder panel
    # Two strips on the left: home temperature, then the FOCAL variant's own genotype in
    # three states. The background variation is pushed right down in alpha -- at ~600
    # variants x 231 founders it is texture, not data, and at full opacity it buried the
    # one column the figure exists to show.
    stripw = (hi - lo) * 0.020
    gap = (hi - lo) * 0.006
    x_gt = lo - (hi - lo) * 0.030                 # focal-genotype strip
    x_t = x_gt - stripw - gap                     # home-temperature strip
    for i in range(nF):
        c = BLUERED(bnorm(b1[i])) if not np.isnan(b1[i]) else "#CCCCCC"
        axH.add_patch(Rectangle((x_t, yrow[i] - 0.5), stripw, 1.0, facecolor=c,
                      edgecolor="none", zorder=2, clip_on=False))

    lc = lead_cols(V, cfg)
    if len(lc):
        alt = (G[:, lc] == 1).any(1)
        miss = MISS[:, lc].all(1)                 # no call in ANY record of the variant
        GT_C = {"alt": "#B2182B", "ref": "#E8E8E8", "miss": "#FFFFFF"}
        for i in range(nF):
            st = "alt" if alt[i] else ("miss" if miss[i] else "ref")
            axH.add_patch(Rectangle((x_gt, yrow[i] - 0.5), stripw, 1.0,
                          facecolor=GT_C[st], edgecolor="#BDBDBD" if st == "miss" else "none",
                          linewidth=0.2, zorder=2, clip_on=False))
        n_alt, n_miss = int(alt.sum()), int(miss.sum())
        n_ref = nF - n_alt - n_miss
        axH.annotate(f"carriers {n_alt}  ·  ref {n_ref}  ·  no call {n_miss}",
                     xy=(x_t, nF + 1.5), xycoords=("data", "data"), fontsize=7,
                     color="#5E5E5E", ha="left", va="bottom", annotation_clip=False)
        # carriers also marked in the body, so their haplotype context is readable
        for i in np.where(alt)[0]:
            axH.plot([cfg["vpos"]], [yrow[i]], marker="D", ms=3.2, color="#B2182B",
                     zorder=6, clip_on=False)

    for cls, (mk, base, col, al) in (
            ("snp", (".", 0.5, "#9E9E9E", 0.16)),
            ("smallindel", ("s", 2.4, "#E67E22", 0.22)),
            ("sv", ("D", 5.0, "#C0392B", 0.35))):
        cols = np.where(V.cls.values == cls)[0]
        xs, ys = [], []
        for j in cols:
            car = np.where(G[:, j] == 1)[0]
            xs.extend([V.pos.values[j]] * len(car)); ys.extend([yrow[i] for i in car])
        axH.scatter(xs, ys, s=base, marker=mk, c=col, alpha=al, edgecolors="none",
                    zorder=3)
    axH.set_ylim(-1, nF); axH.set_yticks([])
    axH.set_ylabel(f"founder (n={nF}) / home temperature", fontsize=9)
    ticks = np.linspace(lo, hi, 6).astype(int)
    axH.set_xticks(ticks); axH.set_xticklabels([f"{t:,}" for t in ticks], fontsize=7)
    axH.set_xlabel(f"{ch} position (bp)", fontsize=9)
    cbb = fig.colorbar(plt.cm.ScalarMappable(cmap=BLUERED, norm=bnorm), ax=axH,
                       fraction=0.02, pad=0.005)
    cbb.set_label("home temperature (°C)", fontsize=7); cbb.ax.tick_params(labelsize=6.5)

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
        axL.set_ylim(yc.min() * 1.02, sx)
        fig.colorbar(plt.cm.ScalarMappable(cmap=LDCMAP, norm=mcolors.Normalize(0, 1)),
                     ax=axL, fraction=0.02, pad=0.005).set_label(r"$r^2$", fontsize=8)
    axL.set_xlim(lo - (hi - lo) * 0.06, hi); axL.axis("off")

    # D (gwas only): the lead across all 30 gardens, cold -> hot.
    # The point of this panel is that p decays SMOOTHLY across gardens, so the
    # recurrence count is wherever a class threshold cuts the gradient -- not a
    # count of independent replications.
    if axG is not None and prof is not None:
        T = prof.sort_values("bio1").reset_index(drop=True)
        x = np.arange(len(T))
        cls_pref = [c for c in ("sv", "nonsnp", "snp") if c in T.columns]
        main = max(cls_pref, key=lambda c: int((T[c] >= BONF_GWAS[c]).sum()))
        axG.bar(x, T[main], width=0.72, zorder=3,
                color=[BLUERED(bnorm(v)) for v in T.bio1], edgecolor="none")
        # only the SV-scan threshold, matching panel A
        axG.axhline(BONF_GWAS[main], ls="--", lw=0.8, color="#b2182b", zorder=2)
        axG.text(-0.65, BONF_GWAS[main], f"{CLSNAME.get(main, main)} Bonf", ha="left",
                 va="bottom", fontsize=6.5, zorder=5, color="#b2182b",
                 bbox=dict(fc="white", ec="none", alpha=0.85, pad=0.8))
        nsig = {c: int((T[c] >= BONF_GWAS[c]).sum()) for c in cls_pref}
        axG.set_xticks(x)
        # garden number over its temperature, so the axis carries both without a caption
        axG.set_xticklabels([f"{int(g)}\n{b:.1f}" for g, b in zip(T.garden, T.bio1)],
                            fontsize=6)
        axG.set_xlabel("garden / temperature (°C)", fontsize=8.5)
        axG.set_ylabel(rf"$-\log_{{10}}p$, {cfg['sym']} lead", fontsize=9)
        axG.set_xlim(-0.8, len(T) - 0.2)
        axG.set_ylim(0, max(T[main].max() * 1.12, BONF_GWAS[main] * 1.15))
        for sp in ("top", "right"):
            axG.spines[sp].set_visible(False)
        axG.tick_params(axis="y", labelsize=7)
        print(f"  garden profile ({main} scan is the widest): " +
              ", ".join(f"{c}={nsig[c]}/{len(T)}" for c in cls_pref))

    # gwas figures get their own suffix so a gene plotted against both scans keeps
    # both, and so a gwas-mode run never silently overwrites a gea-mode figure
    out = f"{HERE}/loci/{cfg['sym']}_combined" + ("_gwas" if scan == "gwas" else "")
    os.makedirs(f"{HERE}/loci", exist_ok=True)
    fig.savefig(out + ".png", dpi=165, bbox_inches="tight")
    fig.savefig(out + ".pdf", bbox_inches="tight")
    print(f"  wrote {out}.png/.pdf")


if __name__ == "__main__":
    main()
