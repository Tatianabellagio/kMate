#!/usr/bin/env python
"""Build+execute the SNP-only vs non-SNP-only class-split GWAS Manhattan notebook (user 2026-07-03).
Load-only (class_gwas_{snp,nonsnp}.npz + summary json); runs in `basic` env (matplotlib hangs in
`plotting`).
"""
import os
import nbformat as nbf
from nbformat.v4 import new_notebook, new_markdown_cell, new_code_cell
from nbconvert.preprocessors import ExecutePreprocessor

ROOT = "/global/scratch/users/tbellg/kmate"
OUT = f"{ROOT}/analysis/grenenet_gea/notebooks/class_split_gwas.ipynb"
os.makedirs(os.path.dirname(OUT), exist_ok=True)

md_title = r"""# SNP-only vs non-SNP-only founder GWAS — same peaks?

Companion to the kinship-level variance-partition work (`varexp_selection.py` / `varexp_untagged.py`):
those asked "do SNP and non-SNP (indel+SV) genomic-relatedness matrices explain the same selection-
trait variance" (yes, even for the SNP-untagged non-SNP layer). This notebook asks the sharper,
per-marker question: run an **actual single-marker association scan** separately on SNP-only vs
non-SNP-only markers, with the *same* trait and the *same* shared kinship correction — do you get
the same **peaks**?

**Important — what kind of GWAS is this?** Each marker is tested **once per site** (31 gardens, each
its own single-site LOCO-EMMAX scan — the per-site GWAS), and those 31 per-site z-scores are then
combined with a **Bolormaa multi-trait meta-analysis** across sites into three genome-wide contrasts:

- **JOINT** (31 df, χ²) — selected at *any* site (an omnibus "was this locus under selection somewhere")
- **GLOBAL** (1 df) — generalist: same-direction effect across all 31 sites
- **CLIMATE** (1 df) — local adaptation: effect tracks site bio1 (mean annual temperature)

So the Manhattans below are the **multi-site (multi-trait) GWAS**, not 31 separate per-site Manhattans
— but since the per-site z-scores are saved (`Z[marker, site]`), the last section also shows what one
individual **per-site** scan looks like, for concreteness.

Design (full detail in `VAREXP_SELECTION_HANDOFF.md`): 212 analyzable founders, the locked selection
coefficient `s` trait (rank-normalized per site), ONE shared LOCO kinship correction per chromosome
built from ALL classes pooled (so results differ only by which markers are *tested*, not by a
different correction), test markers MAC≥5/call≥90% split SNP (1,750,565) vs non-SNP/indel+SV pooled
(536,432)."""

code_load = r"""
import os, json
import numpy as np, pandas as pd
from scipy import stats
import matplotlib.pyplot as plt
plt.rcParams.update({'figure.dpi': 120, 'font.size': 10})
os.chdir("/global/scratch/users/tbellg/kmate")
OUT = "results/grenenet_gea/varexp"
CHROMS = ["chr1", "chr2", "chr3", "chr4", "chr5"]
CLASS_COLORS = {"snp": ["#3b4cc0", "#7aa0c4"], "nonsnp": ["#D55E00", "#f2a679"]}
CLASS_LABEL = {"snp": "SNP-only (M=1,750,565)", "nonsnp": "non-SNP-only: indel+SV (M=536,432)"}
summary = json.load(open(f"{OUT}/class_gwas_summary.json"))

def load_class(name):
    d = np.load(f"{OUT}/class_gwas_{name}.npz", allow_pickle=True)
    df = pd.DataFrame({"chrom": d["chrom"], "pos": d["pos"],
                        "p_joint": d["p_joint"], "q_joint": d["q_joint"],
                        "p_global": d["p_global"], "q_global": d["q_global"],
                        "p_clim": d["p_clim"], "q_clim": d["q_clim"]})
    return df, d["Z"], d["sites"], d["bio1"], d["lam_persite"]

snp_df, snp_Z, sites, bio1, snp_lam = load_class("snp")
nonsnp_df, nonsnp_Z, _, _, nonsnp_lam = load_class("nonsnp")
print(f"{len(sites)} sites | SNP M={len(snp_df):,}  non-SNP M={len(nonsnp_df):,}")
print("SNP:    ", {k: summary["snp"][k] for k in ("lambda_JOINT","lambda_GLOBAL","lambda_CLIM","bonferroni_joint","fdr_joint")})
print("non-SNP:", {k: summary["nonsnp"][k] for k in ("lambda_JOINT","lambda_GLOBAL","lambda_CLIM","bonferroni_joint","fdr_joint")})
"""

code_helpers = r"""
def genome_x(df):
    df = df.sort_values(["chrom", "pos"]).reset_index(drop=True)
    off, centers, x = 0.0, [], np.zeros(len(df))
    for ch in CHROMS:
        mk = (df.chrom == ch).to_numpy()
        if not mk.any():
            continue
        x[mk] = df.pos[mk].to_numpy() + off
        centers.append(off + df.pos[mk].max() / 2)
        off += df.pos[mk].max() * 1.02
    df = df.copy(); df["x"] = x
    return df, centers

def manhattan3(df, centers, colors, title, bonf, s):
    fig, ax = plt.subplots(3, 1, figsize=(12, 9), sharex=True)
    panels = [("p_joint", "JOINT — any-site selection (31df)", "lambda_JOINT", "fdr_joint", "best_p_joint"),
              ("p_global", "GLOBAL — generalist (same direction everywhere)", "lambda_GLOBAL", "fdr_global", None),
              ("p_clim", "CLIMATE — local adaptation (tracks site bio1)", "lambda_CLIM", "fdr_clim", None)]
    for axi, (col, tit, lamk, fdrk, bestk) in zip(ax, panels):
        nlp = -np.log10(np.clip(df[col], 1e-300, 1))
        for i, ch in enumerate(CHROMS):
            mk = (df.chrom == ch).to_numpy()
            axi.scatter(df.x[mk], nlp[mk], s=6, c=colors[i % 2], alpha=.55, edgecolors="none")
        axi.axhline(-np.log10(bonf), color="firebrick", lw=1, ls="--")
        extra = f", best p={s[bestk]:.2g}" if bestk else ""
        axi.set_ylabel(r"$-\log_{10}p$")
        axi.set_title(f"{tit}   |   $\\lambda_{{GC}}$={s[lamk]:.3f}, FDR q<.05: {s[fdrk]}{extra}",
                      loc="left", fontsize=10.5)
        axi.spines[["top", "right"]].set_visible(False)
    ax[2].set_xticks(centers); ax[2].set_xticklabels([c.replace("chr", "Chr") for c in CHROMS])
    ax[2].set_xlabel("genome position")
    fig.suptitle(title, fontsize=13, y=0.995)
    fig.tight_layout(rect=[0, 0, 1, 0.98])
    return fig
"""

code_snp_manhattan = r"""
bonf_snp = 0.05 / len(snp_df)
snp_x, snp_centers = genome_x(snp_df)
fig = manhattan3(snp_x, snp_centers, CLASS_COLORS["snp"],
                 f"SNP-only GWAS ({CLASS_LABEL['snp']}) — multi-site JOINT/GLOBAL/CLIMATE meta, 31 sites",
                 bonf_snp, summary["snp"])
fig.savefig(f"{OUT}/manhattan_snp.png", dpi=150, bbox_inches="tight")
plt.show()
"""

code_nonsnp_manhattan = r"""
bonf_ns = 0.05 / len(nonsnp_df)
ns_x, ns_centers = genome_x(nonsnp_df)
fig = manhattan3(ns_x, ns_centers, CLASS_COLORS["nonsnp"],
                 f"non-SNP-only GWAS ({CLASS_LABEL['nonsnp']}) — multi-site JOINT/GLOBAL/CLIMATE meta, 31 sites",
                 bonf_ns, summary["nonsnp"])
fig.savefig(f"{OUT}/manhattan_nonsnp.png", dpi=150, bbox_inches="tight")
plt.show()
"""

md_overlay = """## Same genome, both classes overlaid — do the peaks line up?

Both JOINT scans on one axis (SNP navy, non-SNP orange). The strongest hit in the whole notebook is
the same 20kb window on Chr2 (~2.34Mb) in *both* independently-run scans."""

code_overlay = r"""
fig, ax = plt.subplots(figsize=(12, 4.2))
snp_xo, centers = genome_x(snp_df)
ns_xo, _ = genome_x(nonsnp_df)
ax.scatter(snp_xo.x, -np.log10(np.clip(snp_xo.p_joint, 1e-300, 1)), s=6, c=CLASS_COLORS["snp"][0],
          alpha=.4, edgecolors="none", label=CLASS_LABEL["snp"])
ax.scatter(ns_xo.x, -np.log10(np.clip(ns_xo.p_joint, 1e-300, 1)), s=6, c=CLASS_COLORS["nonsnp"][0],
          alpha=.5, edgecolors="none", label=CLASS_LABEL["nonsnp"])
ax.axhline(-np.log10(bonf_snp), color="firebrick", lw=1, ls="--", label="SNP Bonferroni")
ax.set_xticks(centers); ax.set_xticklabels([c.replace("chr", "Chr") for c in CHROMS])
ax.set_ylabel(r"$-\log_{10}p_{\rm JOINT}$"); ax.set_xlabel("genome position")
ax.set_title("JOINT (any-site selection) — SNP-only vs non-SNP-only, overlaid", loc="left")
ax.spines[["top", "right"]].set_visible(False); ax.legend(frameon=False, markerscale=3, loc="upper right")
fig.tight_layout()
fig.savefig(f"{OUT}/manhattan_overlay_joint.png", dpi=150, bbox_inches="tight")
plt.show()
print(f"shared 20kb windows={summary['shared_windows_n']}  "
      f"Spearman rho={summary['window_spearman_rho']:.3f} (p={summary['window_spearman_p']:.2g})")
print(f"top-hit window Jaccard: 0.1%={summary['top0.001_jaccard']:.3f}  "
      f"0.5%={summary['top0.005_jaccard']:.3f}  1%={summary['top0.01_jaccard']:.3f}")
"""

md_persite = """## What does ONE per-site scan look like?

The plots above are the 31-site *meta*-analysis (JOINT/GLOBAL/CLIMATE). For concreteness, here is a
single individual per-site LOCO-EMMAX scan (the site with the strongest single-site signal in either
class), plotted the same way for both marker classes — this is what feeds into the meta-analysis
above, one such panel per site, 31 times over."""

code_persite = r"""
best_site_i = int(np.nanargmax(np.maximum(np.nanmax(snp_Z**2, axis=0), np.nanmax(nonsnp_Z**2, axis=0))))
site_id = int(sites[best_site_i])
fig, ax = plt.subplots(2, 1, figsize=(12, 6.5), sharex=True)
for axi, (name, df, Zmat, color) in zip(ax, [("SNP", snp_df, snp_Z, CLASS_COLORS["snp"]),
                                             ("non-SNP", nonsnp_df, nonsnp_Z, CLASS_COLORS["nonsnp"])]):
    z = Zmat[:, best_site_i]
    p = 2 * stats.norm.sf(np.abs(z))
    dfx = df.copy(); dfx["p_site"] = p
    dfx, centers = genome_x(dfx)
    nlp = -np.log10(np.clip(dfx["p_site"], 1e-300, 1))
    for i, ch in enumerate(CHROMS):
        mk = (dfx.chrom == ch).to_numpy()
        axi.scatter(dfx.x[mk], nlp[mk], s=6, c=color[i % 2], alpha=.55, edgecolors="none")
    axi.set_ylabel(r"$-\log_{10}p$")
    axi.set_title(f"site {site_id} — {name} single-site scan (bio1={bio1[best_site_i]:.1f}C)", loc="left")
    axi.spines[["top", "right"]].set_visible(False)
ax[1].set_xticks(centers); ax[1].set_xticklabels([c.replace("chr", "Chr") for c in CHROMS])
ax[1].set_xlabel("genome position")
fig.suptitle(f"One PER-SITE scan (site {site_id}) — not the meta-analysis above", fontsize=12, y=0.995)
fig.tight_layout(rect=[0, 0, 1, 0.97])
fig.savefig(f"{OUT}/manhattan_persite_example.png", dpi=150, bbox_inches="tight")
plt.show()
"""

md_bottom = """## Bottom line

- These Manhattans are the **multi-site (multi-trait) JOINT/GLOBAL/CLIMATE meta-analysis** across
  all 31 gardens (Bolormaa combination of 31 per-site LOCO-EMMAX scans) — not 31 separate per-site
  plots. The per-site panel above shows what one of those 31 inputs looks like on its own.
- **Variance explained is essentially identical between classes**: mean per-site genomic-control
  λ = 1.040 for both SNP and non-SNP scans (independently corroborating the earlier GRM-h² result).
- **CLIMATE is a genuine null for both classes** (0 FDR hits each) — matches the established
  climate-gradient null for this system.
- **JOINT has real, largely-concordant hits**: the single strongest locus genome-wide is the *same*
  20kb window (Chr2 ~2.34Mb) in both independently-run scans; genome-wide peak strength correlates
  significantly between classes (Spearman ρ=0.545); overlap is strongest at the very top hits and
  decays somewhat further down the list (top-0.1% Jaccard 0.67 → top-1% 0.29).
- **Net**: running this GWAS with SNPs alone vs. SVs/indels alone tells essentially the same story —
  same inflation, same dominant locus, substantially (not perfectly) overlapping broader hit set."""

nb = new_notebook(cells=[
    new_markdown_cell(md_title),
    new_code_cell(code_load),
    new_code_cell(code_helpers),
    new_code_cell(code_snp_manhattan),
    new_code_cell(code_nonsnp_manhattan),
    new_markdown_cell(md_overlay),
    new_code_cell(code_overlay),
    new_markdown_cell(md_persite),
    new_code_cell(code_persite),
    new_markdown_cell(md_bottom),
])
ep = ExecutePreprocessor(timeout=900, kernel_name="basic", startup_timeout=180)
ep.preprocess(nb, {"metadata": {"path": os.path.dirname(OUT)}})
with open(OUT, "w") as f:
    nbf.write(nb, f)
print(f"[built+executed] {OUT}")
