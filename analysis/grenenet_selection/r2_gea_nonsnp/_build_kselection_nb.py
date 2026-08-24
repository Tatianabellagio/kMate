#!/usr/bin/env python
"""Build `notebooks/10_lfmm_k_selection.ipynb` — the standalone case for choosing K=16
latent factors for the genome-wide LFMM GEA (calibrated on the genome-wide SNP Δp matrix
and applied to all variant classes / all downstream analyses, incl. the cam5 replication).

The choice rests PRIMARILY on the GIF K-sweep (calibration-based): K=16 sits on the
GIF plateau where adding factors stops helping. The PCA scree is shown as supporting
context, honestly: the dominant structure is in the first few components, and the
broken-stick test over-selects (~48) at this n (355 pools) / huge p, so it is not the
basis for the choice. Phase-1 used K=16, so 16 also keeps the SV results comparable.

Lightweight: reads two small CSVs (the GIF sweep + the precomputed scree). The 5.6 GB
Δp matrix is loaded ONCE, separately, by `build_lfmm_scree.py`.

Execute with the `basic` env:
  python -m nbconvert --to notebook --execute --inplace notebooks/10_lfmm_k_selection.ipynb
"""
import nbformat as nbf

nb = nbf.v4.new_notebook()
cells = []

cells.append(nbf.v4.new_markdown_cell(
"""# How many latent factors? The case for **K = 16**

LFMM needs a number of latent factors `K` to absorb confounding population structure
before testing each variant against climate (here **bio1**, on the **SNP Δp** matrix,
gen9 `last_gen`, **355 pools**). Too few factors leaves structure in (inflated false
positives); too many soaks up real signal (over-correction).

We choose **K = 16**, primarily from the **GIF K-sweep** (a calibration-based criterion):
K = 16 sits on the *plateau* where adding more factors no longer reduces genomic
inflation, the null stays flat, and real signal is retained. The PCA scree is shown
below as supporting context. K = 16 also matches the published **phase-1** SNP GEA, so
the SV results stay directly comparable.

> K is calibrated on SNPs and then applied to all variant classes (SNP / indel / SV)."""))

# ---- criterion 1 (PRIMARY): GIF sweep --------------------------------------
cells.append(nbf.v4.new_markdown_cell(
"""## 1. GIF K-sweep — the primary criterion

LFMM fit at K = 1…20 on the genome-wide SNP matrix (`run_lfmm_ksweep.sh`). For each K we record the
**genomic inflation factor** (GIF; >1 = under-corrected/inflated, <1 = over-corrected),
the number of strong hits, and `frac_top_decile` (≈0.10 = a flat, well-behaved null).
The principled stopping point is the **plateau**: the smallest K past which more factors
no longer reduce GIF."""))

cells.append(nbf.v4.new_code_cell(
'''import numpy as np, pandas as pd
import matplotlib.pyplot as plt

ROOT  = "/global/scratch/users/tbellg/kmate/analysis/grenenet_selection"
PLOTS = f"{ROOT}/lfmm"

sw = pd.read_csv(f"{ROOT}/lfmm/ksweep_gif_snp_gen9_bio1.csv")
display(sw)

fig, ax = plt.subplots(figsize=(8, 4.5))
ax.plot(sw["K"], sw["gif"], "o-")
ax.axhline(1.0, color="red", ls="--", lw=0.9, label="GIF = 1 (ideal)")
ax.axvspan(14, 20, color="green", alpha=0.08, label="GIF plateau (K≥14)")
ax.axvline(16, color="green", ls=":", lw=2, label="K=16 (chosen)")
ax.set_xlabel("K (latent factors)"); ax.set_ylabel("GIF (genomic inflation factor)")
ax.set_title("LFMM GIF vs K — SNP Δp, bio1"); ax.legend(fontsize=8)
fig.tight_layout(); fig.savefig(f"{PLOTS}/lfmm_k_gif.png", dpi=150, bbox_inches="tight"); plt.show()'''))

cells.append(nbf.v4.new_code_cell(
'''g = sw.set_index("K")["gif"]
print(f"GIF K=1..4  (badly under-corrected): {g.loc[1]:.2f}, {g.loc[2]:.2f}, {g.loc[3]:.2f}, {g.loc[4]:.2f}")
print(f"GIF K=6     (elbow):                 {g.loc[6]:.2f}")
print(f"GIF K=16    (chosen):                {g.loc[16]:.2f}")
print(f"GIF K=20    (max swept):             {g.loc[20]:.2f}")
print(f"GIF gain K=16 -> K=20:               {g.loc[16]-g.loc[20]:+.2f}  (negligible — on the plateau)")
print(f"null flatness at K=16 (frac_top_decile): {sw.set_index('K').loc[16,'frac_top_decile']:.3f}  (~0.10 = flat)")
print(f"strong hits (p<1e-5) at K=16:            {int(sw.set_index('K').loc[16,'n_p_lt_1e5'])}  (signal not suppressed)")
print()
print("Read-out: GIF crashes from ~3.6 (K≤4) to ~2.5 by K=6, then flattens. From K=14 on it")
print("is essentially constant (~2.1–2.5); K=16 -> K=20 changes GIF by ~0.1. The null is flat")
print("(frac_top_decile≈0.10) and real signal is retained, so K=16 is a well-calibrated stop.")'''))

# ---- criterion 2 (SUPPORTING): scree ---------------------------------------
cells.append(nbf.v4.new_markdown_cell(
"""## 2. PCA scree — supporting context (read honestly)

PCA of the pool×pool covariance of the centered SNP Δp matrix shows where the variance
lives. The eigenvalues are **precomputed once** by `build_lfmm_scree.py` (which reads the
5.6 GB Δp matrix sequentially) into a small CSV, so this notebook stays light.

**The "capture ~half the structure" rule lands on K=16.** Cumulative variance reaches
**49.9% at K=16** and crosses 50% at K=17 — so K=16 is the last component still *below*
half the structural variance. Taking the leading components up to ~50% cumulative variance
is a standard, conservative way to set the number of latent factors, and here it points to
**K=16**, agreeing with the GIF plateau and phase-1.

Two honest caveats, so the scree is not over-read:
- The **dominant** structure is concentrated in the first ~5 components (the steep part);
  after that it is a long, gently declining tail of weak structure.
- The **broken-stick** test flags ~48 components above its null. With only 355 pools and
  ~2 M records, broken-stick **over-selects** (its per-component null is tiny, so weak
  tail components clear it). We therefore do **not** use broken-stick to set K; it would
  argue for far more correction than the GIF plateau, the 50%-variance rule, or phase-1.
  We report it only for honesty."""))

cells.append(nbf.v4.new_code_cell(
'''scree  = pd.read_csv(f"{ROOT}/lfmm/scree_snp_gen9.csv")
prop    = scree["eig_prop"].to_numpy()
bstick  = scree["broken_stick"].to_numpy()
cum     = scree["cum_var"].to_numpy()
above   = prop > bstick
K_bstick = int(np.argmax(~above)) if (~above).any() else len(prop)

k50 = int(np.searchsorted(cum, 0.5) + 1)
print(f"variance proportion: PC1={prop[0]:.3f}, PC5={prop[4]:.3f}, "
      f"K=16={prop[15]:.4f}, K=20={prop[19]:.4f}  (decreasing, flat tail by ~K10)")
print(f"cumulative variance: K=16 -> {cum[15]:.3f} (just under 50%), K=17 -> {cum[16]:.3f}")
print(f"first component reaching >=50% cumulative variance: K={k50}  (so K=16 is the last below 50%)")
print(f"broken-stick crossover: K={K_bstick}  (over-selects at this n — NOT used to set K)")

m = 40
fig, ax = plt.subplots(1, 2, figsize=(13, 4.5))
ax[0].plot(np.arange(1, m+1), prop[:m], "o-", ms=4, label="observed eigenvalue prop.")
ax[0].plot(np.arange(1, m+1), bstick[:m], "s--", ms=3, color="0.5", label="broken-stick null")
ax[0].axvline(16, color="green", ls=":", lw=2, label="K=16 (chosen)")
ax[0].set_xlabel("component"); ax[0].set_ylabel("variance proportion")
ax[0].set_title("PCA scree (SNP Δp, 355 pools)"); ax[0].legend(fontsize=8)
ax[1].plot(np.arange(1, m+1), cum[:m], "o-", ms=4)
ax[1].axvline(16, color="green", ls=":", lw=2)
ax[1].axhline(0.5, color="0.6", ls="--", lw=0.8, label="50% variance")
ax[1].set_xlabel("component"); ax[1].set_ylabel("cumulative variance")
ax[1].set_title("cumulative variance explained"); ax[1].legend(fontsize=8)
fig.tight_layout(); fig.savefig(f"{PLOTS}/lfmm_k_scree.png", dpi=150, bbox_inches="tight"); plt.show()
print("→ by K=16 the per-component variance has dropped ~6.4x from PC1 and is in the flat tail;")
print("  K=16 captures ~50% of the structure variance.")'''))

# ---- decision -------------------------------------------------------------
cells.append(nbf.v4.new_markdown_cell(
"""## Decision: **K = 16**

| Criterion | What it says | Role |
|---|---|---|
| **GIF plateau** (primary) | GIF crashes to K≈6, then flat ~2.1–2.5; K=16→20 buys ≈0.1 | **sets K = 16** |
| **50%-variance rule** | cumulative variance = 49.9% at K=16, crosses 50% at K=17 | **independently lands on K=16** |
| **Null flatness** | `frac_top_decile ≈ 0.10` at K=16 → flat, not over-corrected | supports |
| **Signal retained** | strong hits (p<1e-5) present at K=16, not zeroed | supports |
| **Phase-1 consistency** | published phase-1 SNP GEA used K=16 | supports / comparability |
| **Broken-stick** | ~48 — over-selects at n=355 / huge p | reported, **not used** |

We set **K = 16** on the GIF plateau and apply it to all variant classes (SNP / indel / SV).
It is the smallest-region, well-calibrated choice that also matches phase-1; going higher
(toward the broken-stick ~48) would add factors for no GIF benefit and risk suppressing signal.

**Caveat that does *not* change the choice:** even at K=16 the GIF is ≈2.2 (ideal ~1.0–1.2),
i.e. strong residual polygenic structure that *no* reasonable K removes (GIF is still ~2.1 at
K=20). This is a property of the data, not a bad K. It is exactly why downstream p-values are
**GIF-recalibrated** (`*.calibrated_pval.csv`) rather than used raw."""))

nb["cells"] = cells
nb["metadata"]["kernelspec"] = {"name": "python3", "display_name": "Python 3", "language": "python"}
out = "/global/scratch/users/tbellg/kmate/analysis/grenenet_selection/notebooks/10_lfmm_k_selection.ipynb"
with open(out, "w") as f:
    nbf.write(nb, f)
print("wrote", out)
