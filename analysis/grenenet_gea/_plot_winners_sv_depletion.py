#!/usr/bin/env python
"""Q2 plot: is the 'winners are SV-depleted / short-read' tendency CLIMATE-graded (which could drive
the hot-site SV purging) or a uniform background? y = per-site winner-SV-depletion / winner-modality,
x = bio1; purge sites highlighted. Flat vs bio1 => cannot explain the climate-graded purging.
Env: basic. Reads winners_sv_depletion.csv. Saves notebooks/plots/winners_sv_depletion.png."""
import numpy as np, pandas as pd, matplotlib.pyplot as plt
from scipy import stats
G = "/global/scratch/users/tbellg/kmate/analysis/grenenet_gea/sv_adaptive"
PLOTS = "/global/scratch/users/tbellg/kmate/analysis/grenenet_gea/notebooks/plots"
df = pd.read_csv(f"{G}/winners_sv_depletion.csv")

panels = [("rho_s_sv", "winner SV-depletion   ρ(founder s, SV-count)\n(negative = winners SV-poor)"),
          ("win_pct_cactus", "fraction of winners that are cactus (long-read)\n(baseline 0.35; low = winners short-read)")]
fig, axes = plt.subplots(1, 2, figsize=(11, 4.6))
for ax, (col, lab) in zip(axes, panels):
    ax.set_axisbelow(True); ax.grid(color="0.9", lw=0.7)
    for sp in ax.spines.values(): sp.set_visible(False)
    base = 0.35 if col == "win_pct_cactus" else 0.0
    ax.axhline(base, color="0.6", lw=0.8, ls=":")
    o = df[~df.purge]; p = df[df.purge]
    ax.scatter(o.bio1, o[col], c="0.7", s=45, label="other sites", zorder=3)
    ax.scatter(p.bio1, p[col], c="#D55E00", s=70, label="purge sites", zorder=4, edgecolor="k", linewidth=0.5)
    b, a0 = np.polyfit(df.bio1, df[col], 1); x = np.array([df.bio1.min(), df.bio1.max()])
    ax.plot(x, a0 + b * x, color="0.35", lw=1.3, ls="--")
    r = stats.spearmanr(df.bio1, df[col])
    ax.annotate(f"vs bio1: ρ={r.statistic:+.2f}  p={r.pvalue:.3f}", xy=(0.03, 0.05),
                xycoords="axes fraction", fontsize=9, color="0.3")
    ax.set_xlabel("site mean annual temperature  bio1 (°C)"); ax.set_ylabel(lab, fontsize=8)
    ax.legend(frameon=False, fontsize=8, loc="upper right")
fig.tight_layout(); fig.savefig(f"{PLOTS}/winners_sv_depletion.png", dpi=130, bbox_inches="tight")
print(f"saved {PLOTS}/winners_sv_depletion.png")
print(f"corr(rho_s_sv, bio1): {stats.spearmanr(df.bio1, df.rho_s_sv).statistic:+.3f} p={stats.spearmanr(df.bio1, df.rho_s_sv).pvalue:.3f}")
print(f"corr(win_pct_cactus, bio1): {stats.spearmanr(df.bio1, df.win_pct_cactus).statistic:+.3f} p={stats.spearmanr(df.bio1, df.win_pct_cactus).pvalue:.3f}")
