#!/usr/bin/env python
"""Standalone single-panel versions of the founder-confound figures, for slides / the paper.

Same data and same house style as `sv_founder_confound.ipynb` -- these are the individual
panels pulled out at presentation size, not a different analysis. Style matches
`notebooks/plots/s_nofilter_climate_scatter.png`: light grid, no spines, plain grey rho/p at
bottom-left, frameless legend, no interpretation inside the plot.

Writes results/sv_adaptive/plots/fig_*.png
"""
import os, sys
import numpy as np, pandas as pd, matplotlib as mpl, matplotlib.pyplot as plt
from scipy import stats
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import lib

G = f"{lib.GEA}/r1_sv_negative_selection/results/sv_adaptive"
P = f"{G}/plots"
os.makedirs(P, exist_ok=True)
plt.rcParams.update({"figure.dpi": 130, "font.size": 9, "axes.linewidth": 0.6})

INS = "#D55E00"; DEL = "#0072B2"; BASE = "0.8"; CAC = "#009E73"; PGC = "#CC79A7"
FIT = "0.35"; ZERO = "0.6"


def style(ax, axis="both"):
    ax.set_axisbelow(True); ax.grid(color="0.88", lw=0.7, axis=axis)
    for sp in ax.spines.values():
        sp.set_visible(False)


def pfmt(p):
    return f"{p:.1e}" if p < 1e-4 else f"{p:.4f}"


def stat(ax, txt, loc=(0.03, 0.03)):
    ax.annotate(txt, xy=loc, xycoords="axes fraction", fontsize=9, color="0.3")


def fitline(ax, x, y, color=FIT):
    b, a0 = np.polyfit(x, y, 1); xs = np.array([x.min(), x.max()])
    ax.plot(xs, a0 + b * xs, color=color, lw=1.5, ls="--", zorder=2)


C = np.load(f"{G}/founder_climate_confound.npz", allow_pickle=True)
S = np.load(f"{G}/founder_sv_content.npz", allow_pickle=True)
GRp = np.load(f"{G}/founder_graph_representation.npz", allow_pickle=True)
TC = np.load(f"{G}/founder_assembly_technology.npz", allow_pickle=True)

keep = C["keep"].astype(bool)
home1 = C["home_bio1"][keep]; g1 = C["gamma_bio1"][keep]
is_cac = C["is_cactus"][keep].astype(bool)
M = {k: S[k][keep].astype(float) for k in S.files if k != "founders"}
KB = M["bp_ins"] / 1000.0
asm_mb = GRp["asm_mb"]; have = np.isfinite(asm_mb)
tech = TC["tech"].astype(str); okt = TC["ok"].astype(bool)
XL = "ecotype origin bio1 (mean annual temp, C)"
YL = "founder kb of inserted sequence"

# ---------------------------------------------------------------- 1. headline scatter
fig, ax = plt.subplots(figsize=(6.2, 5.2))
style(ax)
ax.scatter(home1, KB, c=home1, cmap="coolwarm", s=46, zorder=3, edgecolor="none")
fitline(ax, home1, KB)
ax.set_xlabel(XL); ax.set_ylabel(YL)
r, p = stats.spearmanr(home1, KB)
stat(ax, f"ρ = {r:+.2f}    p = {pfmt(p)}")
fig.tight_layout(); fig.savefig(f"{P}/fig_kb_vs_origin.png", bbox_inches="tight"); plt.close(fig)

# ---------------------------------------------------------------- 2. by panel half
fig, ax = plt.subplots(figsize=(6.2, 5.2))
style(ax)
for m, c, lab in [(is_cac, CAC, f"cactus (n={is_cac.sum()})"),
                  (~is_cac, PGC, f"PanGenie (n={(~is_cac).sum()})")]:
    ax.scatter(home1[m], KB[m], s=32, color=c, alpha=0.8, linewidths=0, label=lab, zorder=3)
    fitline(ax, home1[m], KB[m], color=c)
ax.set_xlabel(XL); ax.set_ylabel(YL)
rc = stats.spearmanr(KB[is_cac], home1[is_cac]).correlation
rp = stats.spearmanr(KB[~is_cac], home1[~is_cac]).correlation
stat(ax, f"cactus ρ = {rc:+.2f}    PanGenie ρ = {rp:+.2f}")
ax.legend(frameon=False, fontsize=8, loc="upper right")
fig.tight_layout(); fig.savefig(f"{P}/fig_kb_vs_origin_by_panel.png", bbox_inches="tight"); plt.close(fig)

# ---------------------------------------------------------------- 3. by platform
fig, ax = plt.subplots(figsize=(6.6, 5.2))
style(ax)
tl = [t for t in ["CLR", "HiFi", "ONT", "ONT_R10.4", "ONT_HiFi"] if (okt & (tech == t)).sum() >= 2]
TCOL = dict(zip(tl, [INS, DEL, CAC, PGC, "#E69F00"]))
for t in tl:
    m = okt & (tech == t)
    ax.scatter(home1[m], KB[m], s=32, color=TCOL[t], alpha=0.8, linewidths=0, label=t, zorder=3)
    if m.sum() >= 6:
        fitline(ax, home1[m], KB[m], color=TCOL[t])
ax.set_xlabel(XL); ax.set_ylabel(YL)
ax.legend(frameon=False, fontsize=8, loc="upper right")
rr = [f"{t} {stats.spearmanr(KB[okt & (tech == t)], home1[okt & (tech == t)]).correlation:+.2f}"
      for t in tl if (okt & (tech == t)).sum() >= 6]
stat(ax, "ρ within platform:  " + "   ".join(rr))
fig.tight_layout(); fig.savefig(f"{P}/fig_kb_vs_origin_by_platform.png", bbox_inches="tight"); plt.close(fig)

# ---------------------------------------------------------------- 4. assembly size
fig, ax = plt.subplots(figsize=(6.2, 5.2))
style(ax)
ax.scatter(home1[have], asm_mb[have], c=home1[have], cmap="coolwarm", s=46, zorder=3, edgecolor="none")
fitline(ax, home1[have], asm_mb[have])
ax.axhline(119.15, color=ZERO, lw=0.8, ls=":")
ax.set_xlabel(XL); ax.set_ylabel("assembly size, 5 chromosomes (Mb)")
r, p = stats.spearmanr(home1[have], asm_mb[have])
stat(ax, f"ρ = {r:+.2f}    p = {pfmt(p)}")
fig.tight_layout(); fig.savefig(f"{P}/fig_assembly_size_vs_origin.png", bbox_inches="tight"); plt.close(fig)

# ---------------------------------------------------------------- 5. measure bars
MEAS = [("n_ins", "ins, count", INS), ("kb_ins", "ins, kb", INS),
        ("n_ins_frac", "ins, count / load", INS), ("bp_ins_frac", "ins, bp per variant", INS),
        ("n_del", "del, count", DEL), ("kb_del", "del, kb", DEL),
        ("n_del_frac", "del, count / load", DEL), ("bp_del_frac", "del, bp per variant", DEL),
        ("n_sv", "SV, count", BASE), ("kb_sv", "SV, kb", BASE),
        ("n_sv_frac", "SV, count / load", BASE), ("bp_sv_frac", "SV, bp per variant", BASE)]
M["kb_ins"] = KB; M["kb_del"] = M["bp_del"] / 1000.0; M["kb_sv"] = M["bp_sv"] / 1000.0
vals = [stats.spearmanr(M[k], home1).correlation for k, _, _ in MEAS]
fig, ax = plt.subplots(figsize=(6.6, 5.2))
style(ax, axis="x")
yy = np.arange(len(MEAS))[::-1]
ax.axvline(0, color=ZERO, lw=0.8, ls=":")
ax.barh(yy, vals, color=[c for _, _, c in MEAS], alpha=0.85, height=0.6, linewidth=0)
ax.set_yticks(yy); ax.set_yticklabels([l for _, l, _ in MEAS], fontsize=8)
ax.set_xlabel("ρ (measure, ecotype origin bio1)")
fig.tight_layout(); fig.savefig(f"{P}/fig_measures_vs_origin_bars.png", bbox_inches="tight"); plt.close(fig)

# ------------------------------------------------- 6-8. geography + climate space
KEY = "/global/scratch/users/tbellg/gea_grene-net/key_files"
fid = C["founders"].astype("U6")[keep]
home18 = C["home_bio18"][keep]
eco = pd.read_csv(f"{KEY}/1001g_regmap_grenet_ecotype_info_corrected_bioclim_2024May16.csv")
eco["ecotypeid"] = eco["ecotypeid"].astype(str); em = eco.set_index("ecotypeid")
BIOS = [f"bio{i}" for i in range(1, 20)]
H = np.column_stack([[float(em.loc[f, b]) for f in fid] for b in BIOS])
Hz = (H - H.mean(0)) / H.std(0)
U, sv_, _ = np.linalg.svd(Hz, full_matrices=False); PCs = U * sv_
varexp = sv_ ** 2 / (sv_ ** 2).sum()

geo = pd.read_csv(f"{KEY}/1001g_regmap_grenet_ecotype_info_corrected_2024May16.csv")
geo["ecotype_id"] = geo["ecotype_id"].astype(str); gm = geo.set_index("ecotype_id")
lat = np.array([float(gm.loc[f, "Latitude_corrected"]) for f in fid])
lon = np.array([float(gm.loc[f, "Longitude_corrected"]) for f in fid])

vmin, vmax = np.quantile(KB, [0.02, 0.98])
norm = mpl.colors.Normalize(vmin, vmax); cmap = mpl.cm.YlGnBu
CBL = "kb of inserted sequence"

fig, ax = plt.subplots(figsize=(11, 5.6))
style(ax)
sc = ax.scatter(lon, lat, c=KB, cmap=cmap, norm=norm, s=50, linewidths=0.3,
                edgecolor="0.35", zorder=3)
ax.set_xlabel("longitude"); ax.set_ylabel("latitude")
ax.set_aspect(1 / np.cos(np.deg2rad(np.nanmean(lat))))
rl = stats.spearmanr(lat, KB)
stat(ax, f"vs latitude:  ρ = {rl.correlation:+.2f}    p = {pfmt(rl.pvalue)}")
cb = fig.colorbar(sc, ax=ax, label=CBL, fraction=0.022, pad=0.012); cb.outline.set_visible(False)
fig.tight_layout(); fig.savefig(f"{P}/fig_map_kb.png", bbox_inches="tight"); plt.close(fig)

fig, ax = plt.subplots(figsize=(6.2, 5.2))
style(ax)
sc = ax.scatter(PCs[:, 0], PCs[:, 1], c=KB, cmap=cmap, norm=norm, s=46,
                linewidths=0.3, edgecolor="0.35", zorder=3)
ax.set_xlabel(f"origin-climate PC1 ({varexp[0]*100:.1f}% var)")
ax.set_ylabel(f"origin-climate PC2 ({varexp[1]*100:.1f}% var)")
r1 = stats.spearmanr(PCs[:, 0], KB)
stat(ax, f"vs PC1:  ρ = {r1.correlation:+.2f}    p = {pfmt(r1.pvalue)}")
cb = fig.colorbar(sc, ax=ax, label=CBL, fraction=0.046, pad=0.02); cb.outline.set_visible(False)
fig.tight_layout(); fig.savefig(f"{P}/fig_climate_pc_kb.png", bbox_inches="tight"); plt.close(fig)

fig, ax = plt.subplots(figsize=(6.2, 5.2))
style(ax)
sc = ax.scatter(home1, home18, c=KB, cmap=cmap, norm=norm, s=46,
                linewidths=0.3, edgecolor="0.35", zorder=3)
ax.set_xlabel(XL); ax.set_ylabel("ecotype origin bio18 (precip warmest quarter, mm)")
rb1 = stats.spearmanr(home1, KB).correlation; rb18 = stats.spearmanr(home18, KB).correlation
stat(ax, f"vs bio1 ρ = {rb1:+.2f}    vs bio18 ρ = {rb18:+.2f}")
cb = fig.colorbar(sc, ax=ax, label=CBL, fraction=0.046, pad=0.02); cb.outline.set_visible(False)
fig.tight_layout(); fig.savefig(f"{P}/fig_origin_climate_kb.png", bbox_inches="tight"); plt.close(fig)

print("[wrote]")
for f in ("fig_kb_vs_origin", "fig_kb_vs_origin_by_panel", "fig_kb_vs_origin_by_platform",
          "fig_assembly_size_vs_origin", "fig_measures_vs_origin_bars",
          "fig_map_kb", "fig_climate_pc_kb", "fig_origin_climate_kb"):
    print(f"  {P}/{f}.png")
