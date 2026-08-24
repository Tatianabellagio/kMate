"""Recreate 07_genomic_density.png as SNP-only, for the lab-meeting 'what SNP-only
studies miss' comparison slide. Mirrors _build_panel_stats_nb.py section 7 exactly
(same binning, centromeres, palette, no titles), but draws only the SNP layer.

Also writes a 2-panel version: SNP-only (top) vs SNP+indel+SV (bottom) per chrom,
so the missed non-SNP layer is visible in one figure.
"""
import os
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

PROJ = "/global/scratch/projects/fc_moilab/tbellg/kmate"
CHROMS = [1, 2, 3, 4, 5]
PLOTDIR = f"{PROJ}/analysis/grenenet_selection/notebooks/plots"

sns.set_theme(style="whitegrid")
plt.rcParams.update({"figure.dpi": 110})
CLASS_ORDER = ["SNP", "indel", "SV"]
_pal = sns.color_palette("deep", 10)
CLASS_COLOR = {"SNP": _pal[0], "indel": _pal[2], "SV": _pal[3]}   # blue, green, red
CENTROMERE_MB = {"Chr1": 15.086, "Chr2": 3.607, "Chr3": 13.588, "Chr4": 3.956, "Chr5": 11.726}
BIN_BP = 200_000

# --- load only meta (pos, ref_len, alt_len) -> class per record ---
rec_chrom, rec_pos, rec_cls = [], [], []
for c in CHROMS:
    cl, Cl = f"chr{c}", f"Chr{c}"
    m = np.load(f"{PROJ}/panel/arch3/{cl}/var_pa_231_arch3_{cl}.meta.npz", allow_pickle=True)
    pos = m["pos"]; rl = m["ref_len"].astype(np.int64); al = m["alt_len"].astype(np.int64)
    ldiff = np.abs(al - rl)
    snp = (rl == 1) & (al == 1)
    sv = ldiff > 50
    indel = (~snp) & (~sv)
    cls = np.where(snp, "SNP", np.where(sv, "SV", "indel"))
    rec_chrom.append(np.full(len(pos), Cl)); rec_pos.append(pos); rec_cls.append(cls)
    print(f"  {Cl}: {len(pos):,} records")

rec = pd.DataFrame({"chrom": np.concatenate(rec_chrom),
                    "pos": np.concatenate(rec_pos),
                    "cls": np.concatenate(rec_cls)})

# --- headline counts for the slide ---
tot = rec.cls.value_counts()
n_snp = int(tot.get("SNP", 0)); n_indel = int(tot.get("indel", 0)); n_sv = int(tot.get("SV", 0))
n_all = len(rec); n_nonsnp = n_indel + n_sv
print("\n=== genome-wide record counts ===")
print(f"  SNP   {n_snp:>10,}")
print(f"  indel {n_indel:>10,}")
print(f"  SV    {n_sv:>10,}")
print(f"  ALL   {n_all:>10,}")
print(f"  non-SNP (missed by SNP-only): {n_nonsnp:,} = {100*n_nonsnp/n_all:.1f}% of the panel "
      f"({n_sv:,} SVs + {n_indel:,} indels)")

# ---------------------------------------------------------------------------
# Figure A: SNP-only density (pair with the original 07_genomic_density.png)
# ---------------------------------------------------------------------------
fig, axes = plt.subplots(len(CHROMS), 1, figsize=(11, 2.0*len(CHROMS)), sharex=False)
for ax, c in zip(axes, CHROMS):
    Cl = f"Chr{c}"
    sub = rec[rec.chrom == Cl]
    bins = np.arange(0, sub.pos.max() + BIN_BP, BIN_BP)
    h, _ = np.histogram(sub.loc[sub.cls == "SNP", "pos"], bins=bins)
    ax.bar(bins[:-1] / 1e6, h, width=BIN_BP/1e6, color=CLASS_COLOR["SNP"], alpha=0.6,
           label="SNP" if c == CHROMS[0] else None, align="edge", linewidth=0)
    ax.axvline(CENTROMERE_MB[Cl], color="black", lw=1.2, ls=(0, (4, 2)),
               label="centromere" if c == CHROMS[0] else None)
    ax.set_ylabel(Cl, rotation=0, ha="right", va="center", fontsize=10)
    ax.grid(axis="x", visible=False)
axes[-1].set_xlabel("position (Mb)")
axes[0].legend(frameon=False, loc="upper right", ncol=2, fontsize=9)
plt.tight_layout()
outA = f"{PLOTDIR}/07_genomic_density_snp_only.png"
fig.savefig(outA, dpi=150, bbox_inches="tight")
print(f"\nwrote {outA}")

# ---------------------------------------------------------------------------
# Figure B: side-by-side per chrom -- SNP-only (left) vs all classes (right),
# shared y per row so the extra indel+SV mass is directly visible.
# ---------------------------------------------------------------------------
fig, axes = plt.subplots(len(CHROMS), 2, figsize=(15, 2.0*len(CHROMS)), sharex="row")
for r, c in enumerate(CHROMS):
    Cl = f"Chr{c}"
    sub = rec[rec.chrom == Cl]
    bins = np.arange(0, sub.pos.max() + BIN_BP, BIN_BP)
    hists = {cls: np.histogram(sub.loc[sub.cls == cls, "pos"], bins=bins)[0] for cls in CLASS_ORDER}
    ymax = max(hists["SNP"].max(), (hists["SNP"] + hists["indel"] + hists["SV"]).max()) * 1.05

    axL, axR = axes[r]
    # left: SNP only
    axL.bar(bins[:-1]/1e6, hists["SNP"], width=BIN_BP/1e6, color=CLASS_COLOR["SNP"],
            alpha=0.85, align="edge", linewidth=0)
    # right: stacked SNP + indel + SV
    bottom = np.zeros_like(hists["SNP"])
    for cls in CLASS_ORDER:
        axR.bar(bins[:-1]/1e6, hists[cls], width=BIN_BP/1e6, bottom=bottom,
                color=CLASS_COLOR[cls], alpha=0.85, align="edge", linewidth=0,
                label=cls if r == 0 else None)
        bottom = bottom + hists[cls]
    for ax in (axL, axR):
        ax.axvline(CENTROMERE_MB[Cl], color="black", lw=1.1, ls=(0, (4, 2)))
        ax.set_ylim(0, ymax)
        ax.grid(axis="x", visible=False)
    axL.set_ylabel(Cl, rotation=0, ha="right", va="center", fontsize=10)
    # corner annotation instead of subplot titles (project convention)
    if r == 0:
        axL.annotate("SNP only", xy=(0.5, 1.02), xycoords="axes fraction",
                     ha="center", va="bottom", fontsize=12, fontweight="bold")
        axR.annotate("SNP + indel + SV  (kMate)", xy=(0.5, 1.02), xycoords="axes fraction",
                     ha="center", va="bottom", fontsize=12, fontweight="bold")
        axR.legend(frameon=False, loc="upper right", ncol=3, fontsize=9)
axes[-1, 0].set_xlabel("position (Mb)")
axes[-1, 1].set_xlabel("position (Mb)")
plt.tight_layout()
outB = f"{PLOTDIR}/07_genomic_density_snp_vs_all.png"
fig.savefig(outB, dpi=150, bbox_inches="tight")
print(f"wrote {outB}")
