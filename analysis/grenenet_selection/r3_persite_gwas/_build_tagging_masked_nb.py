#!/usr/bin/env python
"""Build + execute notebooks/sv_indel_tagging_masked.ipynb (run in the `basic` env).

Corrected, expanded replacement for `sv_snp_ld_tagging.ipynb` / `build_sv_snp_ld.py`.
Two bugs were found in the original SV-only tagging computation (see
`analysis/grenenet_selection/r3_persite_gwas/build_tagging_masked.py` docstring for the full story):
  1. Missingness silently coded as REF (no var_called masking).
  2. No MAC floor -- singleton anchors (53.5% of SVs!) could register spurious
     r2=1 "tags" from small-N coincidence.
This notebook uses the corrected data (`analysis/grenenet_selection/r3_persite_gwas/results/sv_snp_ld_v2/`,
built by `build_tagging_masked.py --mac-floor 2`), extended to indels, and
adds missingness-filter and MAF-filter sensitivity plots -- the combinatorial
matrix of {no filter, missingness filter, MAF filter} x {SV, indel} x
{panel, shortread}.

  BPY=/global/home/users/tbellg/miniforge3/envs/basic/bin/python
  $BPY analysis/grenenet_selection/r3_persite_gwas/_build_tagging_masked_nb.py
"""
import os
import nbformat as nbf
from nbconvert.preprocessors import ExecutePreprocessor

PROJ = "/global/scratch/users/tbellg/kmate"
NBDIR = f"{PROJ}/analysis/grenenet_selection/notebooks"
OUT = f"{NBDIR}/sv_indel_tagging_masked.ipynb"

md_intro = r"""# SV & indel-SNP tagging (corrected) — arch3 panel vs GrENE-Net short-read SNPs

**Replaces** `sv_snp_ld_tagging.ipynb` / `build_sv_snp_ld.py`, which had two bugs
found while reviewing its plots:

1. **Missingness silently coded as REF.** The original loader only used
   `var_pa` (ALT-carrier), never `var_called` (the called mask) -- so a
   missing genotype and a true REF call were indistinguishable, for both the
   anchor (SV/indel) and the candidate SNP. Panel SV records average ~60%
   missing (`panel_stats_arch3.ipynb` &sect;5), and missingness is strongly
   cohort-correlated (long-read vs short-read), so jointly-missing founders
   could inflate apparent agreement.
2. **No MAC floor.** 53.5% of SVs are singletons (MAC=1) -- a singleton
   anchor can register a spurious r&sup2;=1 "tag" against any other marker
   private to the exact same one founder, a small-N coincidence, not real LD.

Confirmed empirically that the GrENE-Net short-read VCF has **zero missing
genotypes** (sampled ~2M GT calls), so only the panel (arch3) side needed
masking. Fixed in `analysis/grenenet_selection/r3_persite_gwas/build_tagging_masked.py`: proper
pairwise-complete masked r&sup2; (vectorized via matrix-vector products, not a
per-SNP loop) + a configurable MAC floor (default &ge;2, i.e. drop singletons)
applied to both the anchor and every candidate SNP. Extended to **indels**
(not just SVs), and every anchor's own `ac`/`an`/`maf` is stored, so the
missingness-filter and MAF-filter sensitivity plots below are computed
post-hoc from one set of `--mac-floor 2` runs -- no separate compute per filter.
"""

code_load = r"""
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

PROJ = "/global/scratch/users/tbellg/kmate"
V2 = f"{PROJ}/analysis/grenenet_selection/r3_persite_gwas/results/sv_snp_ld_v2"
CHROMS = ["Chr1", "Chr2", "Chr3", "Chr4", "Chr5"]
NF = 231

sns.set_theme(style="whitegrid")
plt.rcParams.update({"figure.dpi": 110})

CLASS_ORDER = ["SV", "indel"]
ACCENT = "#2a78d6"   # arch3 panel SNPs
BASE   = "#7f8da6"   # GrENE-Net short-read SNPs
INK_MUTED = "#898781"

frames = []
for cls_key, cls_label in [("sv", "SV"), ("indel", "indel")]:
    for chrom in CHROMS:
        zp = np.load(f"{V2}/tagging_{cls_key}_panel_{chrom}.npz", allow_pickle=True)
        zs = np.load(f"{V2}/tagging_{cls_key}_shortread_{chrom}.npz", allow_pickle=True)
        assert np.array_equal(zp["pos"], zs["pos"])
        frames.append(pd.DataFrame({
            "cls": cls_label, "chrom": chrom, "pos": zp["pos"], "size": zp["size"],
            "ac": zp["ac"], "an": zp["an"], "maf": zp["maf"],
            "r2_panel": zp["best_r2"], "r2_shortread": zs["best_r2"],
        }))
df = pd.concat(frames, ignore_index=True)
df["miss_pct"] = 100 * (1 - df.an / NF)

print(f"{len(df):,} anchors genome-wide (MAC>=2, both classes)")
print(df.cls.value_counts().rename("n"))
print(f"\nmedian anchor missingness: SV {df.loc[df.cls=='SV','miss_pct'].median():.1f}%, "
      f"indel {df.loc[df.cls=='indel','miss_pct'].median():.1f}%")
"""

code_headline = r"""
# 1. Corrected headline tagging rate, SV vs indel, panel vs shortread
# (compare to the OLD buggy numbers, computed without masking/MAC-floor: SV panel/shortread
#  85.8%/47.4% at r2>0.8 -- see sv_snp_ld_tagging.ipynb, now archived/superseded)
rows = []
for cl in CLASS_ORDER:
    sub = df[df.cls == cl]
    for mode, col in [("panel", "r2_panel"), ("shortread", "r2_shortread")]:
        for t in (0.2, 0.5, 0.8, 0.95):
            rows.append({"cls": cl, "snp_source": mode, "threshold": t, "pct": 100*(sub[col] > t).mean()})
headline = pd.DataFrame(rows)
print(headline.pivot_table(index=["cls", "snp_source"], columns="threshold", values="pct").round(1))

fig, axes = plt.subplots(1, 2, figsize=(13, 4.8), sharey=True)
xs = np.linspace(0, 1, 201)
for ax, cl in zip(axes, CLASS_ORDER):
    sub = df[df.cls == cl]
    surv_panel = np.array([(sub.r2_panel >= x).mean() for x in xs])
    surv_short = np.array([(sub.r2_shortread >= x).mean() for x in xs])
    ax.plot(xs, surv_panel, color=ACCENT, lw=2, label="arch3 panel SNPs")
    ax.fill_between(xs, surv_panel, color=ACCENT, alpha=0.10)
    ax.plot(xs, surv_short, color=BASE, lw=2, label="GrENE-Net short-read SNPs")
    ax.fill_between(xs, surv_short, color=BASE, alpha=0.10)
    for t in (0.2, 0.5, 0.8, 0.95):
        ax.axvline(t, color=INK_MUTED, lw=0.8, ls=(0, (3, 3)))
        jp = np.argmin(np.abs(xs - t)); js = jp
        ax.annotate(f"{100*surv_panel[jp]:.0f}%", (t, surv_panel[jp]), xytext=(4, 6),
                    textcoords="offset points", fontsize=8, color=ACCENT, fontweight="bold")
        ax.annotate(f"{100*surv_short[js]:.0f}%", (t, surv_short[js]), xytext=(4, -10),
                    textcoords="offset points", fontsize=8, color=BASE, fontweight="bold")
    ax.annotate(cl, xy=(0.05, 0.08), xycoords="axes fraction", fontweight="bold")
    ax.set_xlabel("r² threshold x")
    ax.grid(axis="x", visible=False)
axes[0].set_ylabel("fraction of anchors with best r² ≥ x")
axes[0].legend(frameon=False, loc="upper right", fontsize=9)
plt.tight_layout(); plt.show()
"""

code_missfilter = r"""
# 2. Missingness-filter sensitivity: tagging rate (r2>0.8) as we require
# increasingly complete anchors. All anchors here already pass MAC>=2; this
# asks whether restricting further to well-called anchors changes the picture.
miss_thresholds = [100, 75, 50, 25, 10]   # max allowed % missing (100 = no extra filter)
rows = []
for cl in CLASS_ORDER:
    for mode, col in [("panel", "r2_panel"), ("shortread", "r2_shortread")]:
        sub = df[df.cls == cl]
        for mt in miss_thresholds:
            keep = sub[sub.miss_pct <= mt]
            rows.append({"cls": cl, "snp_source": mode, "miss_thresh": mt,
                         "pct_tagged": 100*(keep[col] > 0.8).mean(), "n": len(keep)})
missdf = pd.DataFrame(rows)
print(missdf.pivot_table(index=["cls", "snp_source"], columns="miss_thresh", values="pct_tagged").round(1))

fig, axes = plt.subplots(1, 2, figsize=(13, 4.8), sharey=True)
for ax, cl in zip(axes, CLASS_ORDER):
    for mode, color in [("panel", ACCENT), ("shortread", BASE)]:
        sub = missdf[(missdf.cls == cl) & (missdf.snp_source == mode)].sort_values("miss_thresh")
        ax.plot(sub.miss_thresh, sub.pct_tagged, "-o", color=color, lw=2, ms=5,
                label=f"{mode}" if cl == CLASS_ORDER[0] else None)
    ax.set_xlabel("max allowed % missing in anchor (100 = no filter)")
    ax.invert_xaxis()
    ax.annotate(cl, xy=(0.05, 0.92), xycoords="axes fraction", fontweight="bold")
    ax.grid(axis="x", visible=False)
axes[0].set_ylabel("% of anchors with r² > 0.8")
axes[0].legend(frameon=False, loc="lower left", fontsize=9)
plt.tight_layout(); plt.show()
"""

code_maffilter = r"""
# 3. MAF/MAC-filter sensitivity: tagging rate (r2>0.8) as we require
# increasingly common anchors (all already MAC>=2 at compute time).
mac_thresholds = [2, 5, 10, 20]
rows = []
for cl in CLASS_ORDER:
    for mode, col in [("panel", "r2_panel"), ("shortread", "r2_shortread")]:
        sub = df[df.cls == cl]
        mac = np.minimum(sub.ac, sub.an - sub.ac)
        for mt in mac_thresholds:
            keep = sub[mac >= mt]
            rows.append({"cls": cl, "snp_source": mode, "mac_thresh": mt,
                         "pct_tagged": 100*(keep[col] > 0.8).mean(), "n": len(keep)})
macdf = pd.DataFrame(rows)
print(macdf.pivot_table(index=["cls", "snp_source"], columns="mac_thresh", values="pct_tagged").round(1))

fig, axes = plt.subplots(1, 2, figsize=(13, 4.8), sharey=True)
for ax, cl in zip(axes, CLASS_ORDER):
    for mode, color in [("panel", ACCENT), ("shortread", BASE)]:
        sub = macdf[(macdf.cls == cl) & (macdf.snp_source == mode)].sort_values("mac_thresh")
        ax.plot(sub.mac_thresh, sub.pct_tagged, "-o", color=color, lw=2, ms=5,
                label=f"{mode}" if cl == CLASS_ORDER[0] else None)
    ax.set_xlabel("MAC ≥ (anchor)")
    ax.annotate(cl, xy=(0.05, 0.92), xycoords="axes fraction", fontweight="bold")
    ax.grid(axis="x", visible=False)
axes[0].set_ylabel("% of anchors with r² > 0.8")
axes[0].legend(frameon=False, loc="lower right", fontsize=9)
plt.tight_layout(); plt.show()
"""

code_dist = r"""
# 4. Distribution of best r2, SV vs indel, panel vs shortread
fig, axes = plt.subplots(1, 2, figsize=(13, 4.8), sharey=True)
bins = np.linspace(0, 1, 41)
for ax, cl in zip(axes, CLASS_ORDER):
    sub = df[df.cls == cl]
    ax.hist(sub.r2_panel.dropna(), bins=bins, color=ACCENT, alpha=0.85, density=True, label="arch3 panel SNPs")
    ax.hist(sub.r2_shortread.dropna(), bins=bins, color=BASE, alpha=0.55, density=True, label="GrENE-Net short-read SNPs")
    ax.annotate(cl, xy=(0.92, 0.92), xycoords="axes fraction", ha="right", fontweight="bold")
    ax.set_xlabel("best r² to nearest SNP (±50kb)")
    ax.grid(axis="x", visible=False)
axes[0].set_ylabel("density")
axes[0].legend(frameon=False, loc="upper right", fontsize=9)
plt.tight_layout(); plt.show()
"""

code_paired = r"""
# 5. Paired per-anchor comparison: same anchor, panel r2 vs shortread r2
fig, axes = plt.subplots(1, 2, figsize=(11.5, 5.6))
for ax, cl in zip(axes, CLASS_ORDER):
    sub = df[df.cls == cl]
    hb = ax.hexbin(sub.r2_shortread, sub.r2_panel, gridsize=45, mincnt=1,
                    cmap="Blues", bins="log", extent=(0, 1, 0, 1))
    ax.plot([0, 1], [0, 1], color=INK_MUTED, lw=1, ls=(0, (3, 3)))
    ax.set_xlabel("best r² vs GrENE-Net short-read SNPs")
    ax.annotate(cl, xy=(0.05, 0.92), xycoords="axes fraction", fontweight="bold")
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.set_aspect("equal")
    ax.grid(visible=False)
    fig.colorbar(hb, ax=ax, label="anchors (log count)", shrink=0.8)
axes[0].set_ylabel("best r² vs arch3 panel SNPs")
plt.tight_layout(); plt.show()

for cl in CLASS_ORDER:
    sub = df[df.cls == cl]
    n_gain = int(((sub.r2_panel > 0.8) & (sub.r2_shortread <= 0.8)).sum())
    print(f"{cl}: {n_gain:,} of {len(sub):,} ({100*n_gain/len(sub):.1f}%) cross r²>0.8 with arch3 SNPs "
          f"but NOT with GrENE-Net short-read SNPs alone")
"""

md_close = r"""## Takeaways

- **Genome-wide corrected numbers (r&sup2;&gt;0.8, MAC&ge;2):** SV 69.7% (panel) vs
  47.8% (short-read); indel 68.9% (panel) vs 49.7% (short-read). Both bugs
  mattered and pushed in different directions on the Chr1-only preview we
  checked first (panel dropped from the old buggy 85.8%, short-read rose
  from 47.4%), but the genome-wide corrected picture is unambiguous either
  way: **arch3's own SNP catalog tags SVs/indels roughly 20pp better than
  GrENE-Net's short-read-only catalog does**, even after removing every
  artifact we could find.
- **Indels track SVs closely** at the headline threshold (68.9% vs 69.7% for
  panel; 49.7% vs 47.8% for short-read) &mdash; the SV-specific "SNP-only misses
  it" story from the original notebook generalizes to indels too, not just SVs.
- **Missingness filtering matters a lot, and asymmetrically.** Restricting to
  well-called anchors (&le;10% missing) raises panel tagging modestly (SV
  70%&rarr;80%, indel 69%&rarr;82%) but raises short-read tagging much more
  (SV 48%&rarr;63%, indel 50%&rarr;67%). Makes sense: short-read SNPs are
  always fully called (confirmed, zero missing genotypes), so *all* of its
  degradation traces back to anchor missingness &mdash; cleaning that up helps
  it proportionally more than it helps panel mode, which had other sources
  of noise too.
- **MAC filtering beyond the base &ge;2 floor barely matters** &mdash; all four
  lines (SV/indel &times; panel/short-read) are essentially flat from MAC&ge;2
  through MAC&ge;20. The singleton-removal already captured the dominant
  small-N artifact; further rarity filtering doesn't change the picture.
- **Caveat, unchanged from before:** this is genotype-LD tagging (can a SNP
  predict the anchor's presence across founders), not evidence about
  selection or causality.
"""

cells = [
    nbf.v4.new_markdown_cell(md_intro),
    nbf.v4.new_code_cell(code_load.strip()),
    nbf.v4.new_markdown_cell("## 1. Corrected headline tagging rate (survival curve)"),
    nbf.v4.new_code_cell(code_headline.strip()),
    nbf.v4.new_markdown_cell("## 2. Missingness-filter sensitivity"),
    nbf.v4.new_code_cell(code_missfilter.strip()),
    nbf.v4.new_markdown_cell("## 3. MAF/MAC-filter sensitivity"),
    nbf.v4.new_code_cell(code_maffilter.strip()),
    nbf.v4.new_markdown_cell("## 4. Distribution of best r²"),
    nbf.v4.new_code_cell(code_dist.strip()),
    nbf.v4.new_markdown_cell("## 5. Paired per-anchor comparison"),
    nbf.v4.new_code_cell(code_paired.strip()),
    nbf.v4.new_markdown_cell(md_close),
]
nb = nbf.v4.new_notebook(cells=cells, metadata={"kernelspec":
        {"name": "python3", "display_name": "Python 3"}})

if __name__ == "__main__":
    os.makedirs(NBDIR, exist_ok=True)
    ep = ExecutePreprocessor(timeout=1800, kernel_name="python3", startup_timeout=180)
    ep.preprocess(nb, {"metadata": {"path": NBDIR}})
    with open(OUT, "w") as f:
        nbf.write(nb, f)
    print("wrote + executed", OUT)
