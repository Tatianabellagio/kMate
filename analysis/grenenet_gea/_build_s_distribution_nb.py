#!/usr/bin/env python
"""Build + execute the per-site s-distribution notebook (user 2026-07-02): one panel per site,
median selection coefficient s (plot-replicate logit-slope) vs INITIAL-FREQUENCY stratum, for
SNP / indel / SV. No threshold, no null — the distributional Option A. Load-only (reads the
precomputed s_dist_by_stratum.csv/_sitemeta.csv/.npz) so it runs in the `basic` env.
"""
import os
import nbformat as nbf
from nbformat.v4 import new_notebook, new_markdown_cell, new_code_cell
from nbconvert.preprocessors import ExecutePreprocessor

ROOT = "/global/scratch/users/tbellg/kmate"
OUT = f"{ROOT}/analysis/grenenet_gea/notebooks/s_distribution_by_site.ipynb"
os.makedirs(os.path.dirname(OUT), exist_ok=True)

md_title = """# Selection coefficient by class, SNP vs indel vs SV — DE-TRENDED, one panel per site

**Question (distributional, no null):** starting from the same initial frequency, are non-SNP
variants' selection coefficients shifted relative to the other variants there?

**s** = per-variant **plot-replicate** selection coefficient at a site = mean over the site's ~10–12
replicate plots of the OLS slope of $\\mathrm{logit}(p)$ on generation (gen 0 = shared founding
frequency). Plots-as-replicates separate real selection (same direction across independent plots)
from drift. Negative s = declining (purged); positive s = rising (favoured).

**De-trended by initial frequency.** Raw median s vs $p_0$ has a strong NEG→POS slope that is an
**artifact of the logit statistic at the 0/1 boundaries** (a Jensen effect — a zero-selection drift
simulation reproduces it; median *linear* $\\Delta p$ stays ~0), amplified because rare variants sit
on few founders (small effective $N_e$). That trend is shared by all classes, so we **subtract the
per-$p_0$-bin median of ALL variants** (pooled — not SNP-defined). Each panel then shows, per class,
$s_{\\text{class}}-s_{\\text{all}}$ at matched $p_0$: **flat 0 = behaves like same-frequency variants;
below 0 = more purifying selection.** SNP and indel should hug 0; SV's dip (if any) is the signal.

*Caveat:* pool AF is a global-mode founder projection, so a variant's s is a projection of founder-h
slopes; an SV shares its trajectory with SNPs on the same founders (co-occurrence, not proven
SV-specific). Sites ordered cold→hot (bio1)."""

code_load = r"""
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
plt.rcParams.update({'figure.dpi':110, 'font.size':8, 'axes.linewidth':0.6})
G = "/global/scratch/users/tbellg/kmate/results/grenenet_gea/sv_adaptive"
long = pd.read_csv(f"{G}/s_dist_by_stratum.csv")
long["resid"] = long["median"] - long["base"]     # de-trend: class median - ALL-class median per p0 bin
meta = pd.read_csv(f"{G}/s_dist_by_stratum_sitemeta.csv").sort_values("bio1").reset_index(drop=True)
npz  = np.load(f"{G}/s_dist_by_stratum.npz")
p0q = npz["p0q"]; mids = 0.5*(p0q[:-1]+p0q[1:])
# Okabe-Ito colourblind-safe: SNP grey, indel blue, SV vermillion
COL = {"SNP":"#888888", "indel":"#0072B2", "SV":"#D55E00"}
print(f"{len(meta)} sites, {long.stratum.nunique()} initial-frequency strata")
print(meta[['site','bio1','n_plots','shift_sv','shift_ind','wil_p_sv']].to_string(index=False))
"""

code_grid = r"""
sites = meta.site.tolist()
ncol, nrow = 5, int(np.ceil(len(sites)/5))
fig, axes = plt.subplots(nrow, ncol, figsize=(15, 2.5*nrow), sharex=True, sharey=True)
axes = axes.ravel()
for ax in axes[len(sites):]:
    ax.axis("off")
for i, s in enumerate(sites):
    ax = axes[i]; d = long[long.site==s]
    for cls in ["SNP","indel","SV"]:
        dc = d[d.cls==cls].sort_values("stratum")
        if dc.empty: continue
        x = mids[dc.stratum.values]
        ax.plot(x, dc["resid"], "-", color=COL[cls], lw=1.6 if cls=="SV" else 1.0,
                marker="o", ms=3 if cls=="SV" else 2.2, label=cls, zorder=3 if cls=="SV" else 2,
                alpha=0.95 if cls=="SV" else 0.8)
    ax.axhline(0, color="k", lw=0.7, ls=":")           # 0 = behaves like same-frequency variants
    mrow = meta[meta.site==s].iloc[0]
    star = "*" if (mrow.wil_p_sv==mrow.wil_p_sv and mrow.wil_p_sv<0.05) else ""
    ax.set_title(f"site {s}  (bio1 {mrow.bio1:.0f})  ΔSV={mrow.shift_sv:+.2f}{star}",
                 fontsize=7.5)
    ax.set_xscale("log"); ax.set_xlim(mids[0]*0.85, mids[-1]*1.15)
    ax.set_ylim(-0.22, 0.12)
    ax.tick_params(labelsize=6)
fig.supxlabel("initial frequency $p_0$ (log)", y=0.005, fontsize=10)
fig.supylabel("s minus same-frequency baseline   (below 0 = more purged than typical variant)",
              x=0.005, fontsize=10)
fig.legend(handles=[Line2D([0],[0],color=COL[c],lw=2,marker="o",ms=4,label=c)
                    for c in ["SNP","indel","SV"]],
           loc="upper right", ncol=3, fontsize=9, frameon=False, bbox_to_anchor=(0.99,1.005))
fig.suptitle("DE-TRENDED selection coefficient (class − same-$p_0$ baseline) — one panel per site "
             "(cold→hot).  ΔSV = mean SV−baseline shift; * = Wilcoxon p<0.05",
             fontsize=11, y=1.005)
fig.tight_layout(rect=[0.02,0.02,1,0.99])
fig.savefig(f"{G}/s_distribution_by_site.png", dpi=130, bbox_inches="tight")
plt.show()
print("saved s_distribution_by_site.png")
"""

md_read = """### How to read
- **y = 0 (dotted line)** means the class moves like the *typical variant that started at the same
  frequency* — the artifactual $p_0$ trend has been subtracted, so flat-at-0 is the neutral
  expectation.
- **Below 0 = more purifying (more purged); above 0 = more favoured**, relative to same-frequency
  variants.
- **SNP (grey) and indel (blue) hug 0** in essentially every panel — the non-SNP category behaves
  like the bulk. **SV (orange)** is the line to watch; a dip below 0 (mostly at low $p_0$, warmer
  sites) is the modest SV purifying excess. `ΔSV` in each title = SV−baseline averaged over strata."""

code_summary = r"""
fig, ax = plt.subplots(1, 2, figsize=(11, 4))
# (1) pooled de-trended: median-over-sites of (class - baseline) per stratum
piv = (long.groupby(["stratum","cls"])["resid"].median().unstack())
for cls in ["SNP","indel","SV"]:
    ax[0].plot(mids, piv[cls].values, "-o", color=COL[cls], lw=1.8 if cls=="SV" else 1.2,
               ms=4, label=cls)
ax[0].axhline(0, color="k", lw=0.6, ls=":"); ax[0].set_xscale("log")
ax[0].set_xlabel("initial frequency $p_0$")
ax[0].set_ylabel("s − same-$p_0$ baseline (pooled over sites)")
ax[0].set_title("Pooled: SNP/indel ≈ 0; SV dips below at low $p_0$"); ax[0].legend(frameon=False)
# (2) per-site SV excess-vs-baseline vs climate
m = meta
ax[1].axhline(0, color="k", lw=0.6, ls=":")
ax[1].scatter(m.bio1, m.shift_sv, c="#D55E00", s=28, label="SV", zorder=3)
ax[1].scatter(m.bio1, m.shift_ind, c="#0072B2", s=18, alpha=.6, label="indel")
ax[1].set_xlabel("site mean annual temp (bio1)")
ax[1].set_ylabel("mean s − baseline (freq-matched)")
nneg = int((m.shift_sv<0).sum())
ax[1].set_title(f"SV below baseline at {nneg}/{len(m)} sites")
ax[1].legend(frameon=False, title=None)
fig.tight_layout(); fig.savefig(f"{G}/s_distribution_summary.png", dpi=130, bbox_inches="tight")
plt.show()
print(f"SV median excess-vs-baseline across sites = {m.shift_sv.median():+.4f} (below at {nneg}/{len(m)})")
"""

md_take = """### Takeaway
- Once the logit/initial-frequency artifact is removed, **SNP and indel sit on the neutral line
  (0)** at every frequency — the non-SNP category (98% indels) behaves like the bulk.
- **SV sits slightly below 0** (median excess-vs-baseline small and negative, at ~22/31 sites) — a
  modest **excess of purifying selection**, concentrated in rare SVs at warmer sites; SV is not
  systematically *above* 0 anywhere (no positive-selection signal).
- From the threshold analyses this SV dip is **entirely insertion-driven** (svins purged, svdel not)
  — consistent with an ancestral/derived-insertion polarity or founder-calling effect rather than
  selection against SV length per se. Resolving that needs ancestral polarization or an independent
  (local-mode/vg) SV allele frequency, not the founder projection."""

nb = new_notebook(cells=[
    new_markdown_cell(md_title),
    new_code_cell(code_load),
    new_code_cell(code_grid),
    new_markdown_cell(md_read),
    new_code_cell(code_summary),
    new_markdown_cell(md_take),
])
ep = ExecutePreprocessor(timeout=600, kernel_name="basic", startup_timeout=180)
ep.preprocess(nb, {"metadata": {"path": os.path.dirname(OUT)}})
with open(OUT, "w") as f:
    nbf.write(nb, f)
print(f"[built+executed] {OUT}")
