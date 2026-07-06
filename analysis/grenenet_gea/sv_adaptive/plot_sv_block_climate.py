#!/usr/bin/env python3
"""Do the SV-bearing top-JOINT blocks light up gardens of a particular CLIMATE TYPE?

We've scored GEA assuming LINEAR climate (adaptive = up-in-cold/down-in-warm). But the
JOINT-not-GLOBAL / sparse-garden result suggests each block is selected at a small subset of
gardens. This asks whether that subset is climate-structured (a temperature band, a
precip/temp 'type', the two extremes) or scattered -- visually + quantitatively.

Per block = its lead-JOINT hap-cluster marker's per-garden z (kinship-corrected founder-
selection association; z>0 = alt-haplotype founders rose at that garden). Climate = REAL bio1
(temp) and bio12 (precip) per garden from lib.load_climate() (NOT the npz 'bio1', which is
bioPC1 for this tag). Diverging color = z (RdBu, 0=neutral, CVD-safe); |z| as size/outline.

Outputs (results/grenenet_gea/sv_adaptive/):
  svclim_1_heatmap.png      39 SV blocks x 30 gardens, gardens ordered by temp, rows clustered
  svclim_2_climate_space.png small-multiples: convincing blocks, gardens in temp x precip space
  svclim_3_aggregate.png     which gardens/climates get lit across all 39 blocks
Run in `basic` env (matplotlib works there; plotting env hangs). Deterministic.
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from scipy.cluster.hierarchy import linkage, leaves_list
import lib

OUT = Path("results/grenenet_gea/sv_adaptive")
ZTHR, TOPFRAC, ZCLIP = 2.5, 0.01, 4.0
CMAP = "RdBu_r"                                   # diverging, CVD-safe; red=+ (up), blue=- (down)
NORM = TwoSlopeNorm(vcenter=0, vmin=-ZCLIP, vmax=ZCLIP)

raw = lib.multisite_gwas_raw("clq90_pc1")
Z, sites = raw["Z"], raw["sites"]
lead = pd.DataFrame({"unit": raw["unit"], "chi2": raw["chi2_joint"], "row": np.arange(len(raw["unit"]))})
lead = lead.loc[lead.groupby("unit")["chi2"].idxmax()]
u2row, u2chi = dict(zip(lead.unit, lead.row)), dict(zip(lead.unit, lead.chi2))

L = pd.read_csv(OUT / "sv_landscape_clq0.9.csv")
L = L[L.block_id.isin(u2row)].copy()
L["chi2_joint"] = L.block_id.map(u2chi)
ntop = int(TOPFRAC * len(L))
sv = L[(L.chi2_joint.rank(ascending=False) <= ntop) & (L.has_sv == 1)].sort_values("chi2_joint", ascending=False)

# REAL climate per garden (bio1=temp, bio12=precip); npz 'bio1' is bioPC1 -> don't use it
clim = lib.load_climate().reindex(sites)
bio1, bio12 = clim["bio1"].to_numpy(), clim["bio12"].to_numpy()
genes = lib.load_genes()
def gspan(c, s, e):
    g = genes[(genes.chrom == c) & (genes.start <= e) & (genes.end >= s)]
    return ";".join(dict.fromkeys(g.name.fillna(g.gene)))[:22] if len(g) else "(intergenic)"

Zsv = np.vstack([Z[u2row[b]] for b in sv.block_id])          # (39, 30)
labels = [f"{gspan(r.chrom, r.start_pos, r.end_pos)}" for _, r in sv.iterrows()]
n_strong = (np.abs(Zsv) > ZTHR).sum(1)
print(f"{len(sv)} SV-bearing top-{TOPFRAC:.0%} JOINT blocks | {(n_strong>=3).sum()} with >=3 strong gardens")

# ---- quantitative: is each block's z ~ linear temp, quadratic (extremes), or scattered? ----
b1z = (bio1 - bio1.mean()) / bio1.std()
quad = b1z**2 - (b1z**2).mean()
rows = []
for i, (_, r) in enumerate(sv.iterrows()):
    z = Zsv[i]; strong = np.abs(z) > ZTHR
    r_lin = np.corrcoef(z, b1z)[0, 1]
    r_quad = np.corrcoef(z, quad)[0, 1]
    sg_t = bio1[strong]
    rows.append(dict(block=r.block_id, gene=labels[i], n_strong=int(strong.sum()),
                     r_lin_temp=round(r_lin, 2), r_quad_extremes=round(r_quad, 2),
                     strong_temp=f"{sg_t.min():.0f}-{sg_t.max():.0f}C" if strong.any() else "-"))
q = pd.DataFrame(rows)
q.to_csv(OUT / "sv_block_climate_shape.csv", index=False)
nl = (q.r_lin_temp.abs() > 0.4).sum(); nq = (q.r_quad_extremes.abs() > 0.4).sum()
print(f"climate shape of the 39 blocks: |r linear-temp|>0.4: {nl} | |r quadratic-extremes|>0.4: {nq} | "
      f"neither (scattered): {((q.r_lin_temp.abs()<=0.4)&(q.r_quad_extremes.abs()<=0.4)).sum()}")

# ================= FIG 1: heatmap, gardens ordered by temp, rows clustered =================
gi = np.argsort(bio1)                                        # gardens cold -> warm
ro = leaves_list(linkage(Zsv, method="average", metric="correlation")) if len(Zsv) > 2 else np.arange(len(Zsv))
M = Zsv[np.ix_(ro, gi)]
fig, ax = plt.subplots(figsize=(11, 9))
im = ax.imshow(M, cmap=CMAP, norm=NORM, aspect="auto")
# secondary encoding (not color-alone): ring the |z|>2.5 cells
sy, sx = np.where(np.abs(M) > ZTHR)
ax.scatter(sx, sy, s=14, facecolors="none", edgecolors="k", linewidths=0.7)
ax.set_xticks(range(len(gi))); ax.set_xticklabels([f"{int(sites[g])}" for g in gi], fontsize=6, rotation=90)
ax.set_yticks(range(len(ro))); ax.set_yticklabels([labels[r] for r in ro], fontsize=6)
ax.set_xlabel("garden (ordered cold → warm by bio1);  °C on top strip"); ax.set_ylabel("SV-bearing top-JOINT block (gene)")
# top strip: garden temperature
axt = ax.inset_axes([0, 1.01, 1, 0.03])
axt.imshow(bio1[gi][None, :], cmap="magma", aspect="auto"); axt.set_xticks([]); axt.set_yticks([])
axt.set_ylabel("bio1", rotation=0, ha="right", va="center", fontsize=7)
cb = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.02); cb.set_label("founder-selection z at garden  (red = alt-hap rose, blue = fell)")
ax.set_title("Do SV-block selection hits align with temperature?  (○ = |z|>2.5)\n"
             "reds/blues clustering at one end = linear climate; both ends = extremes; scattered = none",
             fontsize=10, loc="left")
fig.tight_layout(); fig.savefig(OUT / "svclim_1_heatmap.png", dpi=150, bbox_inches="tight"); plt.close(fig)

# ============ FIG 2: climate-space small multiples for the convincing (>=3-garden) blocks ====
conv = np.where(n_strong >= 3)[0]
ncol = 4; nrow = int(np.ceil(len(conv) / ncol))
fig, axs = plt.subplots(nrow, ncol, figsize=(3.1 * ncol, 2.9 * nrow), squeeze=False)
for k, i in enumerate(conv):
    a = axs[k // ncol][k % ncol]; z = Zsv[i]
    a.scatter(bio1, bio12, s=8, c="0.8", edgecolors="none", zorder=1)         # all gardens (context)
    st = np.abs(z) > ZTHR
    sc = a.scatter(bio1[st], bio12[st], c=z[st], cmap=CMAP, norm=NORM,
                   s=30 + 45 * (np.abs(z[st]) - ZTHR), edgecolors="k", linewidths=0.5, zorder=3)
    for g in np.where(st)[0]:
        a.annotate(int(sites[g]), (bio1[g], bio12[g]), fontsize=5, ha="center", va="center", zorder=4)
    a.set_title(labels[i], fontsize=7.5, loc="left")
    a.tick_params(labelsize=6); a.spines[["top", "right"]].set_visible(False)
for k in range(len(conv), nrow * ncol):
    axs[k // ncol][k % ncol].axis("off")
fig.supxlabel("bio1  (mean annual temp, °C)", fontsize=9); fig.supylabel("bio12  (annual precip, mm)", fontsize=9)
fig.suptitle("Where each SV block is selected, in climate space (temp × precip). "
             "Grey = all 30 gardens; colored = |z|>2.5 (red up / blue down).\n"
             "Question: do the lit gardens cluster in a climate region (a 'type'), or scatter?",
             fontsize=10)
cb = fig.colorbar(sc, ax=axs, fraction=0.015, pad=0.01); cb.set_label("founder-selection z")
fig.savefig(OUT / "svclim_2_climate_space.png", dpi=150, bbox_inches="tight"); plt.close(fig)

# ================= FIG 3: aggregate — which gardens/climates get lit across all 39 =========
up = (Zsv > ZTHR).sum(0); dn = (Zsv < -ZTHR).sum(0)          # per garden: #blocks up / down
fig, ax = plt.subplots(1, 2, figsize=(12.5, 4.6))
# (a) diverging bars vs temp order
xo = np.argsort(bio1); x = np.arange(len(sites))
ax[0].bar(x, up[xo], color="#b2182b", label="blocks selected UP (z>2.5)")
ax[0].bar(x, -dn[xo], color="#2166ac", label="blocks selected DOWN (z<-2.5)")
ax[0].axhline(0, color="k", lw=.8); ax[0].set_xticks(x); ax[0].set_xticklabels([int(sites[g]) for g in xo], fontsize=6, rotation=90)
ax[0].set_xlabel("garden (cold → warm)"); ax[0].set_ylabel("# of the 39 SV blocks lit here")
ax[0].set_title("(a) Are certain gardens repeatedly lit? (ordered by temp)", fontsize=9.5, loc="left")
ax[0].legend(fontsize=7, frameon=False); ax[0].spines[["top", "right"]].set_visible(False)
# (b) climate space, point size = total blocks lighting the garden, split up/down by color share
tot = up + dn
sc = ax[1].scatter(bio1, bio12, s=20 + 40 * tot, c=(up - dn), cmap=CMAP,
                   norm=TwoSlopeNorm(0, vmin=-max(1, tot.max()), vmax=max(1, tot.max())),
                   edgecolors="k", linewidths=0.5)
for g in range(len(sites)):
    if tot[g] >= 3: ax[1].annotate(int(sites[g]), (bio1[g], bio12[g]), fontsize=6, ha="center", va="center")
ax[1].set_xlabel("bio1 (temp, °C)"); ax[1].set_ylabel("bio12 (precip, mm)")
ax[1].set_title("(b) climate space — size = # blocks lit, color = net up(red)/down(blue)", fontsize=9.5, loc="left")
ax[1].spines[["top", "right"]].set_visible(False)
cb = fig.colorbar(sc, ax=ax[1], fraction=0.04, pad=0.02); cb.set_label("net (up − down) blocks")
fig.tight_layout(); fig.savefig(OUT / "svclim_3_aggregate.png", dpi=150, bbox_inches="tight"); plt.close(fig)

print("\nper-garden lit counts (top by total), temp/precip:")
agg = pd.DataFrame(dict(site=[int(s) for s in sites], bio1=bio1.round(1), bio12=bio12.round(0),
                        up=up, down=dn, total=tot)).sort_values("total", ascending=False)
print(agg.head(12).to_string(index=False))
print(f"\n[done] -> svclim_1_heatmap.png, svclim_2_climate_space.png, svclim_3_aggregate.png, sv_block_climate_shape.csv")
