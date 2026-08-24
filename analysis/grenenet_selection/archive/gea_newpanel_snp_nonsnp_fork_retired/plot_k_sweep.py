#!/usr/bin/env python
"""Plot the K sweep: GIF-vs-K (calibration) and hits-vs-K (power) for the
latent-factor-adjusted partial-rank GEA, per axis, SNP vs non-SNP. Run in basic env.
"""
import os
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = "/global/scratch/users/tbellg/kmate/analysis/grenenet_gea/gea_newpanel/results/k_sweep"
d = pd.read_csv(f"{OUT}/k_sweep.csv")
AXES = ["bio1", "bio18", "pc1"]
CLR = {"snp": "#c1443c", "nonsnp": "#2e7d5b"}

fig, axes = plt.subplots(2, 3, figsize=(14, 8), sharex=True)
for j, ax in enumerate(AXES):
    aG, aH = axes[0, j], axes[1, j]
    for cls in ["snp", "nonsnp"]:
        s = d[(d.cls == cls) & (d.axis == ax)].sort_values("K")
        aG.plot(s.K, s.GIF, "-o", color=CLR[cls], ms=4, label=cls)
        aH.plot(s.K, s.n_bonf_raw, "-o", color=CLR[cls], ms=4, label=f"{cls} Bonf(raw)")
        aH.plot(s.K, s.n_bonf_gifcal, "--s", color=CLR[cls], ms=3, alpha=0.6,
                label=f"{cls} Bonf(GIF-cal)")
    aG.axhline(1.0, ls=":", c="k", lw=0.8); aG.axhspan(0.85, 1.25, color="0.9", zorder=0)
    aG.axvline(16, ls=":", c="0.5", lw=0.8)
    aH.axvline(16, ls=":", c="0.5", lw=0.8)
    aG.set_title(f"{ax}"); aG.set_ylabel("GIF" if j == 0 else "")
    aH.set_ylabel("# Bonferroni blocks/records" if j == 0 else ""); aH.set_xlabel("K latent factors")
    if j == 0:
        aG.legend(fontsize=8, frameon=False); aH.legend(fontsize=7, frameon=False)
fig.suptitle("K sweep — LF-adjusted partial-Spearman: GIF (top) vs hits (bottom). "
             "Dotted vline = K=16 default; grey band = GIF credible [0.85,1.25]", y=0.995)
fig.tight_layout()
fig.savefig(f"{OUT}/k_sweep_figure.png", dpi=140, bbox_inches="tight")
print("wrote", f"{OUT}/k_sweep_figure.png")
