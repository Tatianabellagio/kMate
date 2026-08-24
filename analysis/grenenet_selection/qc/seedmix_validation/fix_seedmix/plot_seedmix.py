"""Visualize per-sample / per-chromosome seed-mix h: production vs fix, and
whether the chromosome-average approaches equimolar 1/231.
Run in the `basic` env (matplotlib). Reads seedmix_allsamples_h.npz.
"""
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = "/global/scratch/users/tbellg/kmate/analysis/grenenet_selection/qc/seedmix_validation/fix_seedmix"
d = np.load(f"{OUT}/seedmix_allsamples_h.npz", allow_pickle=True)
fo = d["founders"].astype(str); chroms = d["chroms"].astype(str)
F = len(fo); u = 1.0 / F
samples = sorted({k.split("__")[1] for k in d.files if k.startswith("agg_fix__")})
NS = len(samples)

# ---- Figure 1: per-chrom scatter + chrom-mean, prod vs fix, for each sample ----
fig, axes = plt.subplots(2, NS, figsize=(3.0 * NS, 6.2), sharey=True)
if NS == 1: axes = axes.reshape(2, 1)
for j, sid in enumerate(samples):
    order = np.argsort(-d[f"agg_fix__{sid}"])       # sort founders by fix aggregate
    for r, tag in enumerate(["prod", "fix"]):
        ax = axes[r, j]
        per = d[f"perchrom_{tag}__{sid}"][:, order]  # 5 x F
        agg = d[f"agg_{tag}__{sid}"][order]
        x = np.arange(F)
        for ci in range(per.shape[0]):
            ax.scatter(x, np.maximum(per[ci], 1e-16), s=2, alpha=0.25, color="#8899bb", rasterized=True)
        ax.plot(x, np.maximum(agg, 1e-16), color="#cc3311", lw=0.8, label="chrom-mean")
        ax.axhline(u, color="k", ls="--", lw=0.8, label="1/231")
        ax.set_yscale("log"); ax.set_ylim(1e-16, 0.1)
        nabs = int((agg < 1e-3).sum())
        ax.set_title(f"{sid} {'PROD' if tag=='prod' else 'FIX'}\nabsorbed={nabs}", fontsize=8)
        if j == 0: ax.set_ylabel(f"{'production' if tag=='prod' else 'fix'}\nh (log)", fontsize=8)
        if r == 1: ax.set_xlabel("founder (sorted)", fontsize=7)
        ax.tick_params(labelsize=6)
axes[0, 0].legend(fontsize=6, loc="lower left")
fig.suptitle("Seed-mix per-founder h: light=per-chromosome, red=chrom-mean, dashed=1/231", fontsize=10)
fig.tight_layout(rect=[0, 0, 1, 0.96])
fig.savefig(f"{OUT}/plots/seedmix_perchrom_prod_vs_fix.png", dpi=130)
print("wrote seedmix_perchrom_prod_vs_fix.png")

# ---- Figure 2: does the average approach 1/231? distribution + absorbed count ----
fig, ax = plt.subplots(1, 3, figsize=(15, 4.2))
# (a) violin of aggregate h per sample, prod vs fix, with 1/231
pos = np.arange(NS)
dataP = [d[f"agg_prod__{s}"] for s in samples]
dataF = [d[f"agg_fix__{s}"] for s in samples]
vp = ax[0].violinplot(dataP, positions=pos - 0.18, widths=0.32, showmedians=True)
vf = ax[0].violinplot(dataF, positions=pos + 0.18, widths=0.32, showmedians=True)
for b in vp["bodies"]: b.set_facecolor("#bbbbbb"); b.set_alpha(0.7)
for b in vf["bodies"]: b.set_facecolor("#33aa55"); b.set_alpha(0.7)
ax[0].axhline(u, color="k", ls="--", lw=1, label="1/231")
ax[0].set_yscale("log"); ax[0].set_ylim(1e-8, 0.05)
ax[0].set_xticks(pos); ax[0].set_xticklabels([s.replace("SEEDMIX_", "") for s in samples], fontsize=8)
ax[0].set_ylabel("aggregate h (log)"); ax[0].set_title("aggregate h per sample (gray=prod, green=fix)")
ax[0].legend(fontsize=8)
# (b) n_absorbed per sample
sm = json.load(open(f"{OUT}/seedmix_allsamples_summary.json"))
nap = [sm[f"{s}_prod"]["n_abs"] for s in samples]; naf = [sm[f"{s}_fix"]["n_abs"] for s in samples]
w = 0.38
ax[1].bar(pos - w/2, nap, w, color="#bbbbbb", label="production")
ax[1].bar(pos + w/2, naf, w, color="#33aa55", label="fix")
ax[1].set_xticks(pos); ax[1].set_xticklabels([s.replace("SEEDMIX_", "") for s in samples], fontsize=8)
ax[1].set_ylabel("# founders absorbed (h<1e-3)"); ax[1].set_title("missing founders per sample")
ax[1].legend(fontsize=8)
# (c) RMSE of aggregate vs 1/231
rp = [sm[f"{s}_prod"]["rmse_vs_uniform"] for s in samples]
rf = [sm[f"{s}_fix"]["rmse_vs_uniform"] for s in samples]
ax[2].bar(pos - w/2, rp, w, color="#bbbbbb", label="production")
ax[2].bar(pos + w/2, rf, w, color="#33aa55", label="fix")
ax[2].set_xticks(pos); ax[2].set_xticklabels([s.replace("SEEDMIX_", "") for s in samples], fontsize=8)
ax[2].set_ylabel("RMSE(aggregate h, 1/231)"); ax[2].set_title("distance of chrom-average from equimolar")
ax[2].legend(fontsize=8)
fig.tight_layout()
fig.savefig(f"{OUT}/plots/seedmix_average_vs_uniform.png", dpi=130)
print("wrote seedmix_average_vs_uniform.png")
