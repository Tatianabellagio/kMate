#!/usr/bin/env python
"""Build the LFMM-K calibration notebook: PCA scree / broken-stick on the SNP Δp matrix
to propose a number of latent factors K for LFMM, and (once the GIF sweep job finishes)
read its results to pick the best-calibrated K. K is calibrated on SNPs and applied to
all variant classes.
"""
import nbformat as nbf

nb = nbf.v4.new_notebook()
cells = []

cells.append(nbf.v4.new_markdown_cell(
"""# LFMM K calibration — how many latent factors?

The phase-1 replication ran LFMM at **K=16**, a value inherited from the *old-panel* phase-1
calibration. This `last_gen` set is only **355 pools**, so 16 factors is ~4.5% of n and may
**over-correct** (suppressing real signal). Here we re-derive K from this run's data, on the
**SNP Δp** matrix, and apply the chosen K to all classes.

Two independent criteria:
1. **PCA scree / broken-stick** on the SNP Δp pool covariance (this notebook) — a structure-based K.
2. **GIF K-sweep** (`run_ksweep.sh` job): LFMM at K=1..20, pick the smallest K with genomic
   inflation factor **GIF ≈ 1** and a flat null p-value histogram. Loaded in the last cell."""))

cells.append(nbf.v4.new_code_cell(
'''import numpy as np, pandas as pd
import matplotlib.pyplot as plt

LFMM = "/global/scratch/users/tbellg/kmate/analysis/grenenet_gea/r2_gea_nonsnp/phase1_replication/results/lfmm"
PLOTS = "/global/home/users/tbellg/scratch/kmate/analysis/grenenet_gea/r2_gea_nonsnp/cam5_replication/plots"
STEM = f"{LFMM}/lfmm_snp_gen9"

np_pools, nr = np.loadtxt(f"{STEM}_dims.txt", dtype=int)
print("SNP Δp matrix:", np_pools, "pools x", f"{nr:,}", "records")

# memmap the C-order [pools x rec] float64 Δp; subsample records for a fast, faithful scree
Y = np.memmap(f"{STEM}_Y.f64", dtype=np.float64, mode="r", shape=(np_pools, nr))
rng = np.random.default_rng(0)
nsub = min(200_000, nr)
cols = np.sort(rng.choice(nr, nsub, replace=False))
Ys = np.asarray(Y[:, cols], dtype=np.float64)        # [pools x nsub]
print("subsampled", f"{nsub:,}", "records for PCA scree")

# center each record across pools (PCA of pool structure, as LFMM sees it)
Yc = Ys - Ys.mean(axis=0, keepdims=True)
C = (Yc @ Yc.T) / Yc.shape[1]                        # [pools x pools] pool covariance
evals = np.linalg.eigvalsh(C)[::-1]                  # descending
evals = np.clip(evals, 0, None)
prop = evals / evals.sum()
cum = np.cumsum(prop)''' ))

cells.append(nbf.v4.new_code_cell(
'''# broken-stick null: expected proportion for component k = (1/p) * sum_{i=k}^{p} 1/i
p = len(evals)
bstick = np.array([np.sum(1.0 / np.arange(k, p + 1)) for k in range(1, p + 1)]) / p
# structure-based K = number of leading components whose eigenvalue exceeds broken-stick
K_bstick = int(np.argmax(prop < bstick))   # first component falling below the null
print(f"broken-stick K = {K_bstick}  (components above the null)")
for crit in (0.5, 0.7, 0.8, 0.9):
    print(f"  components for {int(crit*100)}% variance: {int(np.searchsorted(cum, crit) + 1)}")

fig, ax = plt.subplots(1, 2, figsize=(13, 4.5))
m = 40
ax[0].plot(np.arange(1, m + 1), prop[:m], "o-", ms=4, label="observed eigenvalue prop.")
ax[0].plot(np.arange(1, m + 1), bstick[:m], "s--", ms=3, color="0.5", label="broken-stick null")
ax[0].axvline(K_bstick, color="red", ls=":", label=f"K={K_bstick} (broken-stick)")
ax[0].axvline(16, color="green", ls=":", label="K=16 (phase-1)")
ax[0].set_xlabel("component"); ax[0].set_ylabel("variance proportion"); ax[0].set_title("PCA scree (SNP Δp)")
ax[0].legend(fontsize=8)
ax[1].plot(np.arange(1, m + 1), cum[:m], "o-", ms=4)
ax[1].axvline(K_bstick, color="red", ls=":"); ax[1].axvline(16, color="green", ls=":")
ax[1].axhline(0.8, color="0.6", ls="--", lw=0.8)
ax[1].set_xlabel("component"); ax[1].set_ylabel("cumulative variance"); ax[1].set_title("cumulative variance explained")
fig.tight_layout()
fig.savefig(f"{PLOTS}/lfmm_k_scree.png", dpi=150, bbox_inches="tight")
print("saved", f"{PLOTS}/lfmm_k_scree.png")
plt.show()''' ))

cells.append(nbf.v4.new_markdown_cell(
"""## GIF K-sweep results

Run `sbatch run_ksweep.sh` first (LFMM at K=1..20 on the SNP matrix, recording GIF +
null-flatness). Then this cell reads `ksweep_gif_snp_gen9_bio1.csv`. **Pick the smallest K
with GIF ≈ 1** (GIF > 1 = under-corrected / inflated; GIF < 1 = over-corrected / signal
suppressed). `frac_top_decile` ≈ 0.10 indicates a flat (well-calibrated) null."""))

cells.append(nbf.v4.new_code_cell(
'''import os
KS = "/global/scratch/users/tbellg/kmate/analysis/grenenet_gea/r2_gea_nonsnp/cam5_replication/results/ksweep_gif_snp_gen9_bio1.csv"
if os.path.exists(KS):
    sw = pd.read_csv(KS)
    display(sw)
    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.plot(sw["K"], sw["gif"], "o-")
    ax.axhline(1.0, color="red", ls="--", lw=0.9, label="GIF = 1 (well-calibrated)")
    ax.axvline(16, color="green", ls=":", label="K=16 (phase-1)")
    ax.set_xlabel("K (latent factors)"); ax.set_ylabel("GIF (genomic inflation factor)")
    ax.set_title("LFMM GIF vs K — SNP Δp, bio1"); ax.legend(fontsize=8)
    # smallest K with GIF within 0.05 of 1.0
    ok = sw[(sw["gif"] - 1.0).abs() <= 0.05]
    if len(ok):
        print("smallest K with |GIF-1|<=0.05:", int(ok["K"].min()))
    fig.tight_layout(); fig.savefig(f"{PLOTS}/lfmm_k_gif.png", dpi=150, bbox_inches="tight"); plt.show()
else:
    print("GIF sweep not done yet — run:  sbatch run_ksweep.sh   then re-run this cell.")
    print("expected:", KS)''' ))

nb["cells"] = cells
nb["metadata"]["kernelspec"] = {"name": "python3", "display_name": "Python 3", "language": "python"}
out = "/global/home/users/tbellg/scratch/kmate/analysis/grenenet_gea/r2_gea_nonsnp/cam5_replication/notebooks/lfmm_k_calibration.ipynb"
with open(out, "w") as f:
    nbf.write(nb, f)
print("wrote", out)
