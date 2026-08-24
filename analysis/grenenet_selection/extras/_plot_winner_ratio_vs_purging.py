#!/usr/bin/env python
"""All 30 sites: long-read(cactus) fraction of the top-decile winning ecotypes vs SV purging
(shift_sv = mean s - SNP baseline; more negative = more purged). Directly tests the artifact: do
sites with MORE SV purging have FEWER long-read winners? Env: basic. Saves
notebooks/plots/winner_cactus_fraction_vs_purging.png."""
import numpy as np, pandas as pd, matplotlib.pyplot as plt
from scipy import stats
from adjustText import adjust_text
G = "/global/scratch/users/tbellg/kmate/analysis/grenenet_selection/sv_adaptive"
PLOTS = "/global/scratch/users/tbellg/kmate/analysis/grenenet_selection/notebooks/plots"
df = pd.read_csv(f"{G}/winners_sv_depletion.csv")
sm = pd.read_csv(f"{G}/s_dist_by_stratum_sitemeta.csv")[["site", "shift_sv"]]
df = df.merge(sm, on="site", how="inner")
BASE = 80 / 231.0

fig, ax = plt.subplots(figsize=(8, 5))
ax.set_axisbelow(True); ax.grid(color="0.9", lw=0.7)
for sp in ax.spines.values(): sp.set_visible(False)
ax.axhline(BASE, color="0.5", lw=1.0, ls="--")
ax.annotate(f"panel baseline {BASE:.2f}  (1 cactus : 1.9 PG)", xy=(0.98, BASE),
            xycoords=("axes fraction", "data"), ha="right", va="bottom", fontsize=8, color="0.4")
ax.axvline(0, color="0.7", lw=0.8, ls=":")
o = df[~df.purge]; p = df[df.purge]
ax.scatter(o.shift_sv, o.win_pct_cactus, c="0.7", s=55, zorder=3, label="other sites")
ax.scatter(p.shift_sv, p.win_pct_cactus, c="#D55E00", s=80, zorder=4, edgecolor="k", linewidth=0.6,
           label="purge sites (4,32,43,60,26)")
b, a0 = np.polyfit(df.shift_sv, df.win_pct_cactus, 1); x = np.array([df.shift_sv.min(), df.shift_sv.max()])
ax.plot(x, a0 + b * x, color="0.35", lw=1.4, ls=":")
r = stats.spearmanr(df.shift_sv, df.win_pct_cactus)
ax.annotate(f"ρ = {r.statistic:+.2f}   p = {r.pvalue:.3f}", xy=(0.03, 0.04),
            xycoords="axes fraction", fontsize=9, color="0.35")
texts = [ax.text(rr.shift_sv, rr.win_pct_cactus, str(int(rr.site)), fontsize=6, color="0.45") for _, rr in df.iterrows()]
adjust_text(texts, ax=ax, expand=(1.4, 1.6), arrowprops=dict(arrowstyle="-", color="0.75", lw=0.4))
ax.set_xlabel("SV purging   mean (s − SNP baseline)     ← more purged")
ax.set_ylabel("long-read (cactus) fraction of winners")
ax.legend(frameon=False, fontsize=8, loc="upper left")
fig.tight_layout(); fig.savefig(f"{PLOTS}/winner_cactus_fraction_vs_purging.png", dpi=130, bbox_inches="tight")
print(f"saved {PLOTS}/winner_cactus_fraction_vs_purging.png")
print(f"rho(win_pct_cactus, shift_sv) = {r.statistic:+.3f} p={r.pvalue:.3f}")
