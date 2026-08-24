#!/usr/bin/env python
"""Annotated Manhattan for the GIGANTEA (GI / AT1G22770) non-SNP-only hit: per-site garden-4
non-SNP scan, genome-wide + a Chr1 zoom around GI. Marks the FDR-significant clq0.9 block
(Chr1_4247) and GI's gene body (which sits ~1.5kb upstream, in the +/-2kb flank). Block-level
Bonferroni + BH-FDR threshold lines (matching the notebook tables). Env: basic.
"""
import os, sys
import numpy as np, pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"figure.dpi": 110, "font.size": 9})

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import lib

OUT = f"{lib.GEA}/r3_persite_gwas/results/varexp"
CHROMS = ["chr1", "chr2", "chr3", "chr4", "chr5"]
SITE = 4
GI_START, GI_END = 8061844, 8067716
BLOCK = "Chr1_4247"

d = np.load(f"{OUT}/class_gwas_nonsnp.npz", allow_pickle=True)
si = int(np.where(d["sites"] == SITE)[0][0])
z = d["Z"][:, si]
p = 2 * stats.norm.sf(np.abs(z))
chrom, pos = d["chrom"], d["pos"]
lam = float(np.nanmedian(z ** 2) / stats.chi2.ppf(0.5, 1))

# block-level thresholds (match tables)
chrom_cap = np.array([c.replace("chr", "Chr") for c in chrom])
blk = lib.assign_clq_blocks(chrom_cap, pos, r2=0.9)
keep = blk != ""
lead = lib.collapse_to_blocks(pd.DataFrame({"block": blk[keep],
        "nlp": -np.log10(np.clip(p[keep], 1e-300, 1))}), stat="nlp", block_col="block")
lp = np.sort(10.0 ** (-lead["nlp"].to_numpy())); m = len(lp)
bonf_p = 0.05 / m
crit = lp <= (np.arange(1, m + 1) * 0.05 / m)
fdr_p = float(lp[crit][-1]) if crit.any() else None

# genome x-coords
df = pd.DataFrame({"chrom": chrom, "pos": pos, "p": p}).sort_values(["chrom", "pos"]).reset_index(drop=True)
off, centers, x = 0.0, [], np.zeros(len(df))
for ch in CHROMS:
    mk = (df.chrom == ch).to_numpy()
    x[mk] = df.pos[mk].to_numpy() + off
    centers.append(off + df.pos[mk].max() / 2); off += df.pos[mk].max() * 1.02
df["x"] = x
gi_x = float(df[(df.chrom == "chr1")].pos.iloc[0]) * 0  # placeholder
# GI genome-x = GI midpoint on chr1 (chr1 offset is 0)
gi_mid_x = (GI_START + GI_END) / 2

fig, (axG, axZ) = plt.subplots(2, 1, figsize=(11, 7))

# ---- genome-wide ----
nlp = -np.log10(np.clip(df["p"], 1e-300, 1))
for i, ch in enumerate(CHROMS):
    mk = (df.chrom == ch).to_numpy()
    axG.scatter(df.x[mk], nlp[mk], s=4, c=["#D55E00", "#f2a679"][i % 2], alpha=.5, edgecolors="none", rasterized=True)
axG.axhline(-np.log10(bonf_p), color="firebrick", lw=.9, ls="--", label=f"Bonferroni 0.05")
if fdr_p:
    axG.axhline(-np.log10(fdr_p), color="steelblue", lw=.9, ls=":", label="BH FDR q<0.05")
axG.axvline(gi_mid_x, color="green", lw=1.2, alpha=.8)
axG.annotate("GI / GIGANTEA\n(AT1G22770)", (gi_mid_x, axG.get_ylim()[1] * 0.95),
             color="green", fontsize=8.5, ha="left", va="top")
axG.set_xticks(centers); axG.set_xticklabels([c.replace("chr", "Chr") for c in CHROMS])
axG.set_ylabel(r"$-\log_{10}p$")
axG.set_title(f"non-SNP GWAS, per-site garden {SITE}  (lambda={lam:.3f}) — genome-wide", loc="left", fontsize=10)
axG.spines[["top", "right"]].set_visible(False); axG.legend(frameon=False, fontsize=8, loc="upper right")

# ---- Chr1 zoom around GI ----
W = 120_000
lo, hi = GI_START - W, GI_END + W
mkz = (df.chrom == "chr1") & (df.pos >= lo) & (df.pos <= hi)
z2 = df[mkz]
axZ.scatter(z2.pos, -np.log10(np.clip(z2.p, 1e-300, 1)), s=14, c="#D55E00", alpha=.7, edgecolors="none")
# block span
bmap = pd.read_csv(f"{lib.CLQ_BLOCKS_DIR}/chr1_clq0.9_blocks_clq0.9.tsv", sep="\t").sort_values("start_pos").reset_index(drop=True)
bs, be = int(bmap.iloc[4247].start_pos), int(bmap.iloc[4247].end_pos)
axZ.axvspan(bs, be, color="steelblue", alpha=.18, label=f"sig. block {BLOCK} (FDR)")
axZ.axvspan(GI_START, GI_END, color="green", alpha=.25, label="GI gene body")
axZ.axvspan(GI_END, bs, color="grey", alpha=.12, label=f"flank gap ({bs-GI_END:,} bp)")
if fdr_p:
    axZ.axhline(-np.log10(fdr_p), color="steelblue", lw=.9, ls=":")
axZ.axhline(-np.log10(bonf_p), color="firebrick", lw=.9, ls="--")
# other genes in region
genes = lib.load_genes()
gr = genes[(genes.chrom == "Chr1") & (genes.start <= hi) & (genes.end >= lo)]
for _, gg in gr.iterrows():
    axZ.axvspan(gg.start, gg.end, ymin=0, ymax=0.03, color="black", alpha=.5)
axZ.set_xlim(lo, hi)
axZ.set_xlabel("Chr1 position (bp)"); axZ.set_ylabel(r"$-\log_{10}p$")
axZ.set_title(f"Chr1 zoom: GI is ~{bs-GI_END:,} bp upstream of the FDR-significant block (flank hit, not in-block)",
              loc="left", fontsize=9.5)
axZ.spines[["top", "right"]].set_visible(False); axZ.legend(frameon=False, fontsize=8, loc="upper right")

fig.tight_layout()
fig.savefig(f"{OUT}/plots/gwas_plots/GI_locus_site4_nonsnp.png", dpi=120, bbox_inches="tight")
print(f"lam={lam:.3f} bonf_-log10={-np.log10(bonf_p):.2f} fdr_-log10={-np.log10(fdr_p):.2f}")
print(f"block {BLOCK}: {bs:,}-{be:,}; GI: {GI_START:,}-{GI_END:,}; gap {bs-GI_END:,} bp")
zt = z2.copy(); zt["nlp"] = -np.log10(np.clip(zt.p, 1e-300, 1))
print("top markers in zoom window:")
print(zt.sort_values("nlp", ascending=False).head(5)[["pos", "p", "nlp"]].to_string(index=False))
print(f"wrote {OUT}/gwas_plots/GI_locus_site4_nonsnp.png")
