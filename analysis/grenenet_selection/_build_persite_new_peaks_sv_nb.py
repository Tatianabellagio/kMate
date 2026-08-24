#!/usr/bin/env python
"""Build+execute the SV-focused per-site new-peaks notebook — the SV twin of
persite_new_peaks.ipynb. Same 3-section structure but the STRICT-SV class throughout
(class_gwas_sv.npz vs class_gwas_snp.npz); the pooled non-SNP layer is dropped.

Sections:
  1. Recurrence dot-plot -- ALL SV Bonferroni peaks x gardens, gardens ORDERED and COLOURED by
     bio1 (annual mean T, cold->hot). Ring = SNP-blind (SV-sig here but the SNP scan has NO
     significant lead in the same clq0.9 block). Un-ringed big dots = SV peaks a SNP already tags.
  2. Per-garden mirror Manhattan (SNP up / SV down) for the top gardens by SV new-peak count.
  3. Zoom loci -- top recurrent SV peaks in their strongest garden.

Self-contained (computes inline from class_gwas_{snp,sv}.npz), figures INLINE (no loose PNGs),
runs in `basic` env, no chart titles per repo convention. Writes varexp/persite_sv_{new,all}_peaks.csv.
"""
import os
import nbformat as nbf
from nbformat.v4 import new_notebook, new_markdown_cell, new_code_cell
from nbconvert.preprocessors import ExecutePreprocessor

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
OUT = f"{ROOT}/analysis/grenenet_selection/notebooks/persite_new_peaks_sv.ipynb"
os.makedirs(os.path.dirname(OUT), exist_ok=True)

md_title = r"""# Per-site SV peaks — where SVs beat SNPs, garden by garden (SV twin)

The **SV twin** of `persite_new_peaks.ipynb`: identical structure, but the focal class is the
**strict-SV** scan (`class_gwas_sv.npz`, ~12.8k SV markers) rather than the pooled non-SNP layer.
The multi-trait (JOINT) and variance-partition analyses find SVs **redundant** with SNPs at the
polygenic level (passengers). **Per-garden (single-site) is the right lens** — and there, SVs
surface real peaks the SNP scan misses.

"**New peak**" = a clq0.9 block **significant in the SV scan whose SNPs were NOT significant in the
SNP-only scan** — the plain class-unique-hit definition. Section 1 shows **every SV Bonferroni
peak** (not just the SNP-blind ones): a **ring** marks the cells where the SV block is significant
but the SNP scan has no significant lead in that same LD block (SNP-blind / novel), while a large
un-ringed dot = an SV peak a SNP already tags (redundant). Gardens are **ordered and coloured by
bio1** (annual mean temperature, cold→hot). Honest framing: candidate local-adaptation /
incomplete-SNP-tagging loci, not a confirmed adaptive layer. No chart titles per repo convention."""

code_setup = r'''
import os, sys
import numpy as np, pandas as pd
from scipy import stats
import matplotlib as mpl, matplotlib.pyplot as plt
os.chdir("__ROOT__")
sys.path.insert(0, "analysis/grenenet_selection")
import lib
DIR = "analysis/grenenet_selection/varexp"
CHROM_ORDER = ["Chr1", "Chr2", "Chr3", "Chr4", "Chr5"]
# SNP grey (chrom-alternating); SV terracotta/red (matches the raw_manhattan SV twin)
SNP_GREY = {"Chr1": "#868E96", "Chr2": "#CED4DA", "Chr3": "#868E96", "Chr4": "#CED4DA", "Chr5": "#868E96"}
SNPG = "#B7C0CC"
SV_DK, SV_LT = "#C0392B", "#E8A58C"
SV_RED = {"Chr1": SV_DK, "Chr2": SV_LT, "Chr3": SV_DK, "Chr4": SV_LT, "Chr5": SV_DK}
BIO_CMAP = plt.cm.coolwarm       # cold = blue, hot = red

def load(name):
    d = np.load(f"{DIR}/class_gwas_{name}.npz", allow_pickle=True)
    chrom = np.array([c.replace("chr", "Chr") for c in d["chrom"]]); pos = d["pos"].astype(np.int64)
    blk = lib.assign_clq_blocks(chrom, pos, r2=0.9)
    return dict(chrom=chrom, pos=pos, Z=d["Z"], sites=list(d["sites"]), blk=blk)

SNP, SV = load("snp"), load("sv")      # SNP block-assign is the slow step (~1 min); SV is fast
sites = [int(s) for s in SNP["sites"]]
GENES = lib.load_genes()
def gene_at(c, p):
    g = GENES[GENES.chrom == c]; ov = g[(g.start <= p) & (g.end >= p)]
    r = ov.iloc[0] if len(ov) else g.iloc[(g.start - p).abs().argsort().iloc[0]]
    return r["name"] if isinstance(r["name"], str) and r["name"] else r["gene"]

# bio1 per garden -> ordering (cold->hot) + colour
bio1 = lib.load_climate().reindex(sites)["bio1"].astype(float)
site_bio1 = {int(s): float(bio1.loc[s]) for s in sites}
site_order = sorted(sites, key=lambda s: site_bio1[s])            # cold -> hot
xpos = {s: i for i, s in enumerate(site_order)}
BNORM = mpl.colors.Normalize(min(site_bio1.values()), max(site_bio1.values()))
def bio_color(s): return BIO_CMAP(BNORM(site_bio1[s]))

# genome x-coordinate (for the mirror Manhattans)
chrom_len = {c: int(SNP["pos"][SNP["chrom"] == c].max()) for c in CHROM_ORDER}
off, cum = {}, 0
for c in CHROM_ORDER:
    off[c] = cum; cum += chrom_len[c] + int(2e6)
ticks = [off[c] + chrom_len[c] / 2 for c in CHROM_ORDER]
def gx(chrom, pos): return np.array([off[c] for c in chrom]) + np.asarray(pos)
BONF = {"snp": -np.log10(0.05 / (len(set(SNP["blk"])) - 1)),
        "sv":  -np.log10(0.05 / (len(set(SV["blk"])) - 1))}

def nlp_site(D, si): return -np.log10(np.clip(2 * stats.norm.sf(np.abs(D["Z"][:, si])), 1e-300, 1))
def site_leads(D, si):
    df = pd.DataFrame({"block": D["blk"], "chrom": D["chrom"], "pos": D["pos"], "nlp": nlp_site(D, si)})
    df = df[df.block != ""]
    lead = df.loc[df.groupby("block")["nlp"].idxmax()].reset_index(drop=True)
    p = 10.0 ** (-lead["nlp"]); lead["fdr"] = lib.bh(p.to_numpy()) < 0.05; lead["bonf"] = p < 0.05 / len(p)
    return lead

# Two per (garden) block-level collections vs the SNP scan:
#   P (new)  -> SV block significant, SNPs NOT significant (both FDR & Bonferroni)
#   A (all)  -> every SV Bonferroni-significant block, flagged by whether SNPs are also sig (snp_bonf)
rows, arows = [], []
for si, sid in enumerate(sites):
    snp = site_leads(SNP, si).rename(columns={"nlp": "snp_nlp", "fdr": "snp_fdr", "bonf": "snp_bonf"})
    m = site_leads(SV, si).merge(snp[["block", "snp_nlp", "snp_fdr", "snp_bonf"]], on="block")
    g = m.nlp - m.snp_nlp
    fdr_new = m.fdr & (~m.snp_fdr)
    bonf_new = m.bonf & (~m.snp_bonf)
    for i in m.index[fdr_new | bonf_new]:
        r = m.loc[i]
        rows.append(dict(site=int(sid), block=r.block, chrom=r.chrom, pos=int(r.pos),
                         sv_nlp=round(r.nlp, 2), snp_nlp=round(r.snp_nlp, 2), gap=round(g[i], 2),
                         fdr_new=bool(fdr_new[i]), bonf_new=bool(bonf_new[i]), bio1=site_bio1[int(sid)]))
    for i in m.index[m.bonf]:
        r = m.loc[i]
        arows.append(dict(site=int(sid), block=r.block, chrom=r.chrom, pos=int(r.pos),
                          sv_nlp=round(r.nlp, 2), snp_nlp=round(r.snp_nlp, 2),
                          snp_bonf=bool(r.snp_bonf), bio1=site_bio1[int(sid)]))
P = pd.DataFrame(rows); P["gene"] = [gene_at(c, p) for c, p in zip(P.chrom, P.pos)]
A = pd.DataFrame(arows); A["gene"] = [gene_at(c, p) for c, p in zip(A.chrom, A.pos)]
P.to_csv(f"{DIR}/persite_sv_new_peaks.csv", index=False)
A.to_csv(f"{DIR}/persite_sv_all_peaks.csv", index=False)
for tier in ["fdr_new", "bonf_new"]:
    t = P[P[tier]]
    print(f"SV {tier.split('_')[0].upper():4s} new: {len(t)} instances / {t.site.nunique()} gardens / {t.block.nunique()} loci")
nov = A[~A.snp_bonf]
print(f"SV BONF all: {len(A)} instances / {A.block.nunique()} loci "
      f"({nov.block.nunique()} with a SNP-blind (novel) instance, {A.block.nunique() - nov.block.nunique()} always SNP-tagged)")
'''
code_setup = code_setup.replace("__ROOT__", ROOT)

md_recur = r"""## 1. Recurrence dot-plot (gardens ordered & coloured by bio1) — all SV Bonferroni peaks

**Rows = every SV peak** — each clq0.9 block Bonferroni-significant in the **SV-only** scan in ≥1
garden. Columns = the 30 gardens **ordered cold→hot by bio1** and **coloured by bio1** (blue=cold,
red=hot; colourbar at right). Dot size ∝ −log10p in that garden. **A black ring marks cells where
the SV block is significant but the SNP scan has NO significant lead in that same LD block** — the
SNP-blind / novel peaks. So: **ringed dot** = SV peak here that SNPs miss (novel); **large un-ringed
dot** = SV peak here that a SNP already tags (redundant / passenger). A row ringed across many
gardens is a recurrent SNP-blind locus; whether rings sit on the blue (cold) or red (hot) side
shows any climate patterning."""

code_recur = r'''
def recurrence():
    s = A                                            # every Bonferroni-significant SV block
    blocks = (s.groupby("block").agg(gene=("gene", "first"), chrom=("chrom", "first"),
              pos=("pos", "first"), n=("site", "nunique")).sort_values("n", ascending=False).reset_index())
    ring = P[P.bonf_new]                             # SNP-blind cells (SV-sig, SNPs non-sig)
    genuine = set(zip(ring.block, ring.site))
    lead = {}
    for _, b in blocks.iterrows():
        mask = SV["blk"] == b.block
        lead[b.block] = {sid: (nlp_site(SV, si)[mask].max() if mask.any() else 0.0) for si, sid in enumerate(sites)}
    fig, ax = plt.subplots(figsize=(12, 0.42 * len(blocks) + 1.4))
    for yi, (_, b) in enumerate(blocks.iterrows()):
        for sid in site_order:
            xi = xpos[sid]; v = lead[b.block][sid]
            ax.scatter(xi, yi, s=8 + v * 22, color=bio_color(sid), linewidths=0, zorder=2)
            if (b.block, sid) in genuine:
                ax.scatter(xi, yi, s=8 + v * 22, facecolors="none", edgecolors="k", linewidths=1.3, zorder=3)
    ax.set_yticks(range(len(blocks)))
    ax.set_yticklabels([f"{b.gene}\n{b.chrom}:{b.pos/1e6:.2f}Mb" for _, b in blocks.iterrows()], fontsize=8)
    ax.set_xticks(range(len(site_order))); ax.set_xticklabels([str(s) for s in site_order], fontsize=7, rotation=90)
    ax.set_xlabel("garden (site), ordered cold → hot by bio1"); ax.set_ylim(-0.6, len(blocks) - 0.4)
    ax.annotate("SV — all Bonferroni peaks: dot size ∝ −log10p    ○ ringed = SV-sig here but SNP scan not sig in the same block",
                (0.0, 1.02), xycoords="axes fraction", fontsize=8, va="bottom")
    cb = fig.colorbar(mpl.cm.ScalarMappable(norm=BNORM, cmap=BIO_CMAP), ax=ax, pad=0.01, fraction=0.025)
    cb.set_label("bio1 (°C, annual mean T)", fontsize=8)
    for sp in ["top", "right"]:
        ax.spines[sp].set_visible(False)
    fig.tight_layout(); plt.show()

recurrence()
'''

md_mirror = r"""## 2. Per-garden mirror Manhattan — SNP↑ / SV↓

For the gardens carrying the most SV new peaks: SNP −log10p **up** (grey), SV **down** (terracotta),
shared genome axis. The SV track is naturally sparse (~12.8k markers) — which is the point: SVs are
rare but bring these specific peaks. Grey dashed verticals mark that garden's SV new peaks
(Bonferroni), gene-labelled. **At the grey lines the SNP (up) side is low** — the honest "SV peak
where SNPs are blind" view. Dashed horizontal lines = block-level Bonferroni per class."""

code_mirror = r'''
byg = P[P.bonf_new].groupby("site").block.nunique().sort_values(ascending=False)
if not len(byg):
    byg = P[P.fdr_new].groupby("site").block.nunique().sort_values(ascending=False)
GARDENS = list(byg.head(4).index)                    # top gardens by SV new-peak count
for g in GARDENS:
    si = sites.index(g); pk = P[(P.site == g) & (P.bonf_new)].copy()
    pk["gx"] = gx(pk.chrom.values, pk.pos.values)
    snl, vnl = nlp_site(SNP, si), nlp_site(SV, si)
    sgx, vgx = gx(SNP["chrom"], SNP["pos"]), gx(SV["chrom"], SV["pos"])
    fig, ax = plt.subplots(figsize=(13, 5.5))
    for _, r in pk.iterrows():
        ax.axvline(r.gx, color="#6c757d", lw=0.8, ls="--", alpha=0.75, zorder=0)
    ax.scatter(sgx, snl, s=3, c=[SNP_GREY[c] for c in SNP["chrom"]], rasterized=True, linewidths=0, zorder=2)
    ax.scatter(vgx, -vnl, s=6, c=[SV_RED[c] for c in SV["chrom"]], rasterized=True, linewidths=0, alpha=0.9, zorder=2)
    ax.axhline(BONF["snp"], ls="--", lw=.8, c="0.4"); ax.axhline(-BONF["sv"], ls="--", lw=.8, c=SV_DK)
    ax.axhline(0, lw=.8, c="k")
    for _, r in pk.iterrows():
        ax.annotate(r.gene, (r.gx, -r.sv_nlp), fontsize=7, ha="center", va="top",
                    xytext=(0, -3), textcoords="offset points", color=SV_DK, zorder=5)
    ax.set_xticks(ticks); ax.set_xticklabels(CHROM_ORDER)
    ax.set_ylabel("−log10 p   (SNP ↑    SV ↓)")
    yl = max(snl.max(), vnl.max()) * 1.1; ax.set_xlim(0, cum); ax.set_ylim(-yl, yl)
    ax.annotate(f"garden {g} (bio1 {site_bio1[g]:.1f}°C)  —  {pk.block.nunique()} SV new peaks "
                f"(Bonferroni); grey lines = SV peaks, SNPs non-sig there",
                (0.01, 0.02), xycoords="axes fraction", fontsize=9, color="#495057", va="bottom")
    for sp in ["top", "right"]:
        ax.spines[sp].set_visible(False)
    fig.tight_layout(); plt.show()
'''

md_zoom = r"""## 3. Zoom loci — the strongest recurrent SV peaks

Top recurrent SV new-peaks in their strongest garden: local ±60 kb window, SNP (grey) vs SV
(terracotta) −log10p, gene models underneath, the SV lead ringed. SNPs flat, SV spiking."""

code_zoom = r'''
sv = P[P.bonf_new]
if not len(sv):
    sv = P[P.fdr_new]
top = (sv.groupby("block").agg(gene=("gene", "first"), chrom=("chrom", "first"),
       n=("site", "nunique")).sort_values("n", ascending=False).head(4).reset_index())
fig, axes = plt.subplots(1, len(top), figsize=(4.2 * len(top), 4.3), squeeze=False); axes = axes[0]
PAD = 60_000
for ax, (_, b) in zip(axes, top.iterrows()):
    best = sv[sv.block == b.block].sort_values("sv_nlp", ascending=False).iloc[0]
    si = sites.index(int(best.site)); c = b.chrom; center = int(best.pos); lo, hi = center - PAD, center + PAD
    for D, col, sz, z in [(SNP, SNPG, 10, 1), (SV, SV_DK, 45, 3)]:
        w = (D["chrom"] == c) & (D["pos"] >= lo) & (D["pos"] <= hi)
        ax.scatter(D["pos"][w] / 1e6, nlp_site(D, si)[w], s=sz, c=col, linewidths=0, zorder=z,
                   label="SNP" if col == SNPG else "SV")
    ax.scatter([center / 1e6], [best.sv_nlp], s=110, facecolors="none", edgecolors="k", linewidths=1.2, zorder=4)
    for _, gg in GENES[(GENES.chrom == c) & (GENES.start <= hi) & (GENES.end >= lo)].iterrows():
        ax.plot([max(gg.start, lo) / 1e6, min(gg.end, hi) / 1e6], [-0.5, -0.5], lw=4, c="#495057", solid_capstyle="butt")
    ax.set_xlabel(f"{c} (Mb)"); ax.set_ylabel("−log10 p"); ax.set_ylim(-1, max(best.sv_nlp, 3) * 1.15)
    ax.annotate(f"{b.gene}\ngarden {int(best.site)} (bio1 {site_bio1[int(best.site)]:.1f}°C)\nSV {best.sv_nlp:.1f} vs SNP {best.snp_nlp:.1f}",
                (0.03, 0.97), xycoords="axes fraction", ha="left", va="top", fontsize=8.5, color=SV_DK, zorder=5)
    for sp in ["top", "right"]:
        ax.spines[sp].set_visible(False)
axes[-1].legend(loc="upper right", frameon=False, markerscale=1.5, fontsize=8)
fig.tight_layout(); plt.show()
'''

md_bottom = r"""## Takeaway

At the **per-garden** level SVs genuinely bring peaks the SNP scan misses — Bonferroni-level at
several loci and **recurrent across independent gardens** (not FDR noise). Invisible to the
aggregate SNP story because **garden-specific (local)** — consistent with the polygenic/passenger
result at the meta level. With gardens ordered by bio1, inspect whether the ringed (SNP-blind) peaks
concentrate at hot vs cold gardens (climate patterning) or spread evenly. Honest framing: candidate
**local-adaptation / incomplete-SNP-tagging** loci, not a genome-wide adaptive layer (no GO
enrichment). `AT5G04170`/CML50 (calcium sensor) is the standout recurrent SV peak. This is the SV
twin of `persite_new_peaks.ipynb`; the pooled non-SNP layer lives there."""

nb = new_notebook(cells=[
    new_markdown_cell(md_title), new_code_cell(code_setup),
    new_markdown_cell(md_recur), new_code_cell(code_recur),
    new_markdown_cell(md_mirror), new_code_cell(code_mirror),
    new_markdown_cell(md_zoom), new_code_cell(code_zoom),
    new_markdown_cell(md_bottom),
])
ep = ExecutePreprocessor(timeout=2400, kernel_name="basic", startup_timeout=180)
ep.preprocess(nb, {"metadata": {"path": os.path.dirname(OUT)}})
with open(OUT, "w") as f:
    nbf.write(nb, f)
print(f"[built+executed] {OUT}")
