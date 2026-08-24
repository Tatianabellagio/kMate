#!/usr/bin/env python
"""Build + execute notebooks/sv_snp_ld_tagging.ipynb (run in the `basic` env).

Genome-wide writeup of the arch3-panel SV-SNP LD tagging comparison computed by
`build_sv_snp_ld.py` -> analysis/grenenet_gea/sv_snp_ld/sv_snp_ld_{panel,shortread}_Chr{1-5}.npz.
No notebook previously existed for these numbers (they only lived as raw .npz + a
partial SLURM log). The only prior write-up of this question ran on the pre-arch3
merged panel and is now archived (`archive/preprocess_qc/notebooks/ld_*.ipynb`) --
this notebook is its arch3-era replacement, restricted to the SV-SNP tagging
question (not the SNP-SNP LD-decay comparison, which doesn't carry over cleanly
since arch3 changed the SNP catalog itself).

  BPY=/global/home/users/tbellg/miniforge3/envs/basic/bin/python
  $BPY analysis/grenenet_gea/_build_sv_snp_ld_nb.py
"""
import os
import nbformat as nbf
from nbconvert.preprocessors import ExecutePreprocessor

PROJ = "/global/scratch/users/tbellg/kmate"
NBDIR = f"{PROJ}/analysis/grenenet_gea/notebooks"
OUT = f"{NBDIR}/sv_snp_ld_tagging.ipynb"

md_intro = r"""# SV-SNP LD tagging — arch3 panel vs GrENE-Net short-read SNPs

**Question:** in the current production (arch3) 231-founder panel, how well are
structural variants (>50bp) tagged by the nearest SNP within &plusmn;50kb — and how
much of that tagging is actually coming from arch3's own denser SNP catalog, vs
what the old GrENE-Net SNP-only catalog (`greneNet_final_v1.1`, no SVs) could
already see?

Two r&sup2; scans, same SV set, same founders, same &plusmn;50kb window, computed by
`analysis/grenenet_gea/build_sv_snp_ld.py`:

- **panel** &mdash; each SV's max founder-genotype r&sup2; against SNPs from arch3's own
  merged catalog (`var_pa_231_arch3`)
- **shortread** &mdash; the same SVs against SNPs from the old GrENE-Net SNP-only VCF
  (`greneNet_final_v1.1`, founder-aligned)

r&sup2; = squared Pearson correlation of the SV-presence vector vs SNP alt-dosage
vector over the 231 founders. All 5 chromosomes, genome-wide.

**Note on the prior analysis:** an earlier version of this question was answered
in `preprocess_qc/notebooks/ld_comparison.ipynb` / `ld_full_comparison.ipynb` /
`ld_sv_snp_histograms.ipynb`, but on the **pre-arch3** merged panel
(`founders_231_chr.vcf.gz`). That analysis is now archived
(`archive/preprocess_qc/`) since arch3 replaced that panel build; this notebook
answers the same question on current production data.
"""

code_load = r"""
import numpy as np, pandas as pd
import matplotlib.pyplot as plt

PROJ = "/global/scratch/users/tbellg/kmate"
LD = f"{PROJ}/analysis/grenenet_gea/sv_snp_ld"
CHROMS = ["Chr1", "Chr2", "Chr3", "Chr4", "Chr5"]

# --- style: light-mode, quiet chrome (dataviz skill defaults) ---
ACCENT = "#2a78d6"   # arch3 panel SNPs (blue, palette slot 1) -- the "new" series
BASE   = "#7f8da6"   # GrENE-Net short-read SNPs (muted gray-blue) -- the "old" series
INK        = "#0b0b0b"
INK_MUTED  = "#898781"
GRID       = "#e1e0d9"
plt.rcParams.update({
    "figure.facecolor": "#fcfcfb", "axes.facecolor": "#fcfcfb",
    "axes.edgecolor": INK_MUTED, "axes.labelcolor": INK,
    "text.color": INK, "xtick.color": INK_MUTED, "ytick.color": INK_MUTED,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.8,
    "axes.spines.top": False, "axes.spines.right": False,
    "font.size": 11, "figure.dpi": 110,
})

def size_class(sv_size):
    sv_size = np.asarray(sv_size)
    cls = np.full(sv_size.shape, "large_sv", dtype=object)
    cls[sv_size < 5000] = "medium_sv"
    cls[sv_size < 500] = "small_sv"
    return cls

rows = []
for chrom in CHROMS:
    zp = np.load(f"{LD}/sv_snp_ld_panel_{chrom}.npz", allow_pickle=True)
    zs = np.load(f"{LD}/sv_snp_ld_shortread_{chrom}.npz", allow_pickle=True)
    assert np.array_equal(zp["pos"], zs["pos"]), "panel/shortread SV sets must align"
    rows.append(pd.DataFrame({
        "chrom": chrom, "pos": zp["pos"], "sv_size": zp["sv_size"],
        "sv_class": size_class(zp["sv_size"]),
        "r2_panel": zp["best_r2"], "r2_shortread": zs["best_r2"],
    }))
df = pd.concat(rows, ignore_index=True)
print(f"{len(df):,} SVs genome-wide (Chr1-5), arch3 production panel")
print(df.sv_class.value_counts().rename("n"))
"""

code_dist = r"""
# 1. Distribution of best r2 -- panel vs shortread
fig, ax = plt.subplots(figsize=(9, 4.8))
bins = np.linspace(0, 1, 41)
ax.hist(df.r2_panel.dropna(),     bins=bins, color=ACCENT, alpha=0.85, label="arch3 panel SNPs",        density=True)
ax.hist(df.r2_shortread.dropna(), bins=bins, color=BASE,   alpha=0.55, label="GrENE-Net short-read SNPs", density=True)
for t in (0.2, 0.5, 0.8, 0.95):
    ax.axvline(t, color=INK_MUTED, lw=0.8, ls=(0, (3, 3)))
ax.set_xlabel("best r² to nearest SNP (±50kb)"); ax.set_ylabel("density")
ax.legend(frameon=False, loc="upper left")
ax.grid(axis="x", visible=False)
plt.tight_layout(); plt.show()
"""

code_survival = r"""
# 2. Survival curve: fraction of SVs with r2 >= x -- continuous version of the
# headline threshold table (the discrete 4-bar version was redundant with this
# and was dropped; the same threshold values are now annotated directly here).
xs = np.linspace(0, 1, 201)
surv_panel     = np.array([(df.r2_panel     >= x).mean() for x in xs])
surv_shortread = np.array([(df.r2_shortread >= x).mean() for x in xs])

thresholds = [0.2, 0.5, 0.8, 0.95]
headline = pd.DataFrame({
    "threshold": thresholds,
    "panel_pct":     [100*surv_panel[np.argmin(np.abs(xs - t))] for t in thresholds],
    "shortread_pct": [100*surv_shortread[np.argmin(np.abs(xs - t))] for t in thresholds],
})
print(headline.round({"panel_pct": 1, "shortread_pct": 1}).to_string(index=False))

fig, ax = plt.subplots(figsize=(9, 4.8))
ax.plot(xs, surv_panel,     color=ACCENT, lw=2, label="arch3 panel SNPs")
ax.fill_between(xs, surv_panel, color=ACCENT, alpha=0.10)
ax.plot(xs, surv_shortread, color=BASE,   lw=2, label="GrENE-Net short-read SNPs")
ax.fill_between(xs, surv_shortread, color=BASE, alpha=0.10)
for _, row in headline.iterrows():
    t = row.threshold
    ax.axvline(t, color=INK_MUTED, lw=0.8, ls=(0, (3, 3)))
    ax.annotate(f"{row.panel_pct:.0f}%", (t, row.panel_pct/100), xytext=(4, 6),
                textcoords="offset points", fontsize=8.5, color=ACCENT, fontweight="bold")
    ax.annotate(f"{row.shortread_pct:.0f}%", (t, row.shortread_pct/100), xytext=(4, -10),
                textcoords="offset points", fontsize=8.5, color=BASE, fontweight="bold")
ax.set_xlabel("r² threshold x"); ax.set_ylabel("fraction of SVs with best r² ≥ x")
ax.legend(frameon=False, loc="upper right")
ax.grid(axis="x", visible=False)
plt.tight_layout(); plt.show()
"""

code_paired = r"""
# 3. Paired per-SV comparison: same SV, panel r2 vs shortread r2
fig, ax = plt.subplots(figsize=(6.2, 6))
hb = ax.hexbin(df.r2_shortread, df.r2_panel, gridsize=45, mincnt=1,
                cmap="Blues", bins="log", extent=(0, 1, 0, 1))
ax.plot([0, 1], [0, 1], color=INK_MUTED, lw=1, ls=(0, (3, 3)))
ax.set_xlabel("best r² vs GrENE-Net short-read SNPs")
ax.set_ylabel("best r² vs arch3 panel SNPs")
cb = fig.colorbar(hb, ax=ax, label="SVs (log count)", shrink=0.85)
ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.set_aspect("equal")
ax.grid(visible=False)
plt.tight_layout(); plt.show()

gain = df.r2_panel - df.r2_shortread
n_gain_bin = int(((df.r2_panel > 0.8) & (df.r2_shortread <= 0.8)).sum())
print(f"{n_gain_bin:,} SVs ({100*n_gain_bin/len(df):.1f}%) cross the r²>0.8 'well-tagged' "
      f"bar with arch3 SNPs but NOT with GrENE-Net short-read SNPs alone --")
print("these are exactly the SVs a SNP-only (GrENE-Net-style) GEA would under-tag.")
"""

code_untagged = r"""
# 4. Where are the untagged SVs (r2 <= 0.2, essentially untagged), genome-wide?
# Small multiples per chromosome, rug-style: is the tagging gap concentrated
# near centromeres/repeat deserts, or scattered? arch3's untagged set is much
# smaller than GrENE-Net's, so a shared row (not shared y-scale) per source
# keeps both legible.
CENTROMERE_MB = {"Chr1": 15.086, "Chr2": 3.607, "Chr3": 13.588, "Chr4": 3.956, "Chr5": 11.726}
UNTAG_T = 0.2

untagged_sr = df[df.r2_shortread <= UNTAG_T]
untagged_pn = df[df.r2_panel <= UNTAG_T]
print(f"untagged (r2<={UNTAG_T}) vs GrENE-Net short-read SNPs: {len(untagged_sr):,} SVs "
      f"({100*len(untagged_sr)/len(df):.1f}%)")
print(f"untagged (r2<={UNTAG_T}) vs arch3 panel SNPs:          {len(untagged_pn):,} SVs "
      f"({100*len(untagged_pn)/len(df):.1f}%)")

fig, axes = plt.subplots(5, 1, figsize=(11, 10), sharex=False)
for ax, chrom in zip(axes, CHROMS):
    sr = untagged_sr.loc[untagged_sr.chrom == chrom, "pos"] / 1e6
    pn = untagged_pn.loc[untagged_pn.chrom == chrom, "pos"] / 1e6
    ax.eventplot(sr, lineoffsets=1, linelengths=0.8, color=BASE,
                 label="untagged: GrENE-Net short-read SNPs" if chrom == CHROMS[0] else None)
    ax.eventplot(pn, lineoffsets=2, linelengths=0.8, color=ACCENT,
                 label="untagged: arch3 panel SNPs" if chrom == CHROMS[0] else None)
    ax.axvline(CENTROMERE_MB[chrom], color="black", lw=1.2, ls=(0, (4, 2)),
               label="centromere" if chrom == CHROMS[0] else None)
    ax.set_yticks([1, 2]); ax.set_yticklabels(["short-read", "arch3"], fontsize=8)
    ax.set_ylabel(chrom, rotation=0, ha="right", va="center", fontsize=10)
    ax.set_ylim(0.3, 2.7)
    ax.grid(axis="y", visible=False)
axes[-1].set_xlabel("position (Mb)")
axes[0].legend(frameon=False, loc="upper right", ncol=3, fontsize=8)
plt.tight_layout(); plt.show()
"""

md_close = r"""## Takeaways

- **Within the arch3 panel**, SVs are almost perfectly tagged by the panel's own
  SNPs (median r&sup2; &asymp; 1.0; 85.8% at r&sup2;&gt;0.8 genome-wide) &mdash; consistent with
  the small (231-founder), heavily-related, selfing GrENE-Net panel having much
  stronger genome-wide LD than natural diversity panels report (e.g. Kang 2023's
  ~17% tagged at r&sup2;&gt;0.6 in outbred *Arabidopsis* populations,
  `docs/GEA_SV_LITERATURE_REVIEW.md`).
- **Restricted to the old GrENE-Net SNP-only catalog**, tagging drops sharply:
  47.4% at r&sup2;&gt;0.8, vs 85.8% with arch3's own SNPs. Roughly **1 in 3 SVs**
  that arch3's denser SNP catalog would call well-tagged would NOT have been
  caught as well-tagged by a SNP-only (pre-kMate) analysis.
- This numeric result feeds the broader "are SVs GEA passengers or an
  independent selection currency" thread
  (`analysis/grenenet_gea/VAREXP_SELECTION_HANDOFF.md`, driver/passenger
  investigation) &mdash; that work uses a stricter within-panel indel+SV tagging scan
  (`build_nonsnp_tagging.py`) and reaches a compatible ~73-86% tagged-at-r&sup2;&gt;0.8
  number.
- **Caveat:** this is genotype-LD tagging (can a SNP predict the SV's presence
  across founders), not evidence about selection or causality &mdash; see the
  driver/passenger and variance-partition docs for that question.
"""

cells = [
    nbf.v4.new_markdown_cell(md_intro),
    nbf.v4.new_code_cell(code_load.strip()),
    nbf.v4.new_markdown_cell("## 1. Distribution of best r²"),
    nbf.v4.new_code_cell(code_dist.strip()),
    nbf.v4.new_markdown_cell("## 2. Survival curve (headline tagging rate, continuous)"),
    nbf.v4.new_code_cell(code_survival.strip()),
    nbf.v4.new_markdown_cell("## 3. Paired per-SV comparison"),
    nbf.v4.new_code_cell(code_paired.strip()),
    nbf.v4.new_markdown_cell("## 4. Where are the untagged SVs?"),
    nbf.v4.new_code_cell(code_untagged.strip()),
    nbf.v4.new_markdown_cell(md_close),
]
nb = nbf.v4.new_notebook(cells=cells, metadata={"kernelspec":
        {"name": "python3", "display_name": "Python 3"}})

if __name__ == "__main__":
    os.makedirs(NBDIR, exist_ok=True)
    ep = ExecutePreprocessor(timeout=900, kernel_name="python3", startup_timeout=180)
    ep.preprocess(nb, {"metadata": {"path": NBDIR}})
    with open(OUT, "w") as f:
        nbf.write(nb, f)
    print("wrote + executed", OUT)
