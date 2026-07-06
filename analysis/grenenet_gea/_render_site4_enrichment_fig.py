#!/usr/bin/env python
"""Render the site-4 SV-enrichment pilot figure from precomputed plot-data.

Part A: per-variant median |s| vs founding-freq p0, one line per variant class.
Part B: fraction of blocks carrying an SV vs block selection-score decile.
Part C: size-matched null histogram of SV-bearing fraction with selected &
        non-selected vlines; title shows the fold and empirical p-value.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = "/global/scratch/users/tbellg/kmate/results/grenenet_gea/site_temporal"
IN = f"{BASE}/site4_sv_enrichment_plotdata.npz"
OUT = f"{BASE}/site4_sv_enrichment.png"

d = np.load(IN)

fig, (axA, axB, axC) = plt.subplots(1, 3, figsize=(15, 4.5))

# Part A --------------------------------------------------------------
axA.plot(d["A_binmid"], d["A_med_snp"], marker="o", label="SNP")
axA.plot(d["A_binmid"], d["A_med_smallindel"], marker="s", label="small indel")
axA.plot(d["A_binmid"], d["A_med_sv"], marker="^", label="SV")
axA.set_xscale("log")
axA.set_xlabel("founding frequency p0")
axA.set_ylabel("median |s|")
axA.set_title("A: p0-matched median |s| by class")
axA.legend()

# Part B --------------------------------------------------------------
axB.plot(d["B_decile"], d["B_sv_bearing"], marker="o", color="tab:purple")
axB.set_xlabel("block selection-score decile")
axB.set_ylabel("fraction of blocks bearing an SV")
axB.set_title("B: SV-bearing vs selection decile")

# Part C --------------------------------------------------------------
axC.hist(d["C_null_frac"], bins=40, color="0.7", edgecolor="none")
f_sel = float(d["C_f_sel"])
f_non = float(d["C_f_non"])
axC.axvline(f_sel, color="tab:red", lw=2, label=f"selected ({f_sel:.3f})")
axC.axvline(f_non, color="tab:blue", lw=2, ls="--", label=f"non-selected ({f_non:.3f})")
axC.set_xlabel("SV-bearing fraction (size-matched null)")
axC.set_ylabel("count")
axC.set_title(f"C: fold={float(d['fold']):.2f}, p_emp={float(d['p_emp']):.2f}")
axC.legend()

fig.tight_layout()
fig.savefig(OUT, dpi=130)
print("wrote", OUT)
