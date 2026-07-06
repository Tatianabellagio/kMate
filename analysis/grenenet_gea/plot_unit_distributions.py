#!/usr/bin/env python
"""4-panel distribution figure for the final dynamic-LD unit map: variants/unit, panel
k-mers/unit, within-unit LD (mean r2), and PC1-VE (gen9 pools). Covered vs desert overlaid.
Run in the plotting/basic env."""
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

D = pd.read_csv("results/grenenet_gea/blocks_mcf90/final_units_dynld_K500_metrics.csv")
cov = D[D.covered]; des = D[~D.covered]
fig, ax = plt.subplots(2, 2, figsize=(13, 9))

# 1: variants/unit (log)
b = np.logspace(0, 3, 45)
ax[0, 0].hist(D.n_variants, bins=b, color="#C44E52", alpha=.85)
ax[0, 0].set_xscale("log"); ax[0, 0].axvline(D.n_variants.median(), ls="--", color="k", lw=1)
ax[0, 0].set_xlabel("# variants / unit"); ax[0, 0].set_ylabel("# units")
ax[0, 0].set_title(f"Variants per unit (median {int(D.n_variants.median())})")

# 2: panel k-mers/unit (log), threshold line
b2 = np.logspace(0, 4.5, 45)
ax[0, 1].hist(des.panel_kmers.clip(lower=1), bins=b2, color="#999999", alpha=.7, label="desert")
ax[0, 1].hist(cov.panel_kmers, bins=b2, color="#4C72B0", alpha=.75, label="covered")
ax[0, 1].set_xscale("log"); ax[0, 1].axvline(500, ls="--", color="red", lw=1.5, label="K=500")
ax[0, 1].set_xlabel("panel k-mers / unit"); ax[0, 1].set_ylabel("# units")
ax[0, 1].set_title(f"Panel k-mers per unit (median {int(D.panel_kmers.median())})"); ax[0, 1].legend()

# 3: within-unit LD r2
ax[1, 0].hist(D.mean_r2.dropna(), bins=40, color="#55A868", alpha=.85)
ax[1, 0].axvline(D.mean_r2.median(), ls="--", color="k", lw=1)
ax[1, 0].set_xlabel("mean within-unit LD (founder r²)"); ax[1, 0].set_ylabel("# units")
ax[1, 0].set_title(f"Within-unit LD (median r² {D.mean_r2.median():.2f})")

# 4: PC1-VE on gen9 pools
ax[1, 1].hist(des.pc1_ve.dropna(), bins=40, color="#999999", alpha=.7, label="desert")
ax[1, 1].hist(cov.pc1_ve.dropna(), bins=40, color="#DD8452", alpha=.8, label="covered")
ax[1, 1].axvline(0.7, ls="--", color="red", lw=1.2)
ax[1, 1].set_xlabel("PC1-VE (gen9 pools)"); ax[1, 1].set_ylabel("# units")
ax[1, 1].set_title(f"PC1 explained variance (median {D.pc1_ve.median():.2f})"); ax[1, 1].legend()

plt.suptitle(f"Final dynamic-LD units (K=500), n={len(D):,}  |  covered {int(D.covered.sum()):,} ({100*D.covered.mean():.0f}%)",
             fontsize=13)
plt.tight_layout(); plt.savefig("results/grenenet_gea/blocks_mcf90/unit_distributions.png", dpi=110)
print("wrote results/grenenet_gea/blocks_mcf90/unit_distributions.png")
for c in ["n_variants", "panel_kmers", "mean_r2", "pc1_ve"]:
    v = D[c].dropna()
    print(f"  {c:>12}: median {v.median():.3f}  q25 {v.quantile(.25):.3f}  q75 {v.quantile(.75):.3f}")
