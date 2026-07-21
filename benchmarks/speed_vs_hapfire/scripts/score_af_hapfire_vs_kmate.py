#!/usr/bin/env python3
"""hapFIRE-vs-kMate per-variant (AF) accuracy, matched pools / native panels --
the AF-side counterpart to score_hapfire_vs_kmate.py's h-accuracy comparison.

Both tools get the SAME pool composition (identical founders-meta + seed, see
score_hapfire_vs_kmate.py's docstring) but score against a DIFFERENT, tool-native
truth: kMate against arch3's SNP records (already computed by
poolsize_depth/scripts/score_snp_vs_nonsnp.py -> snp_vs_nonsnp_table.tsv),
hapFIRE against the greneNet VCF's SNP genotypes projected through the SAME
pool_weights.tsv used for its own greneNet-fasta sim (extract_greneNet_gt_matrix.py
-> work/greneNet_gt_matrix.npz). This is NOT a shared per-site truth -- arch3 and
greneNet are different variant-calling lineages, so the two tools are scored
against different SNP sets of different sizes (see n_snps columns below). It IS
the fair, matched-pool comparison: same founder draw, each tool's own natural
truth, mirroring how the h-accuracy figure is already built.

kMate's snp_R2/snp_RMSE come from score_snp_vs_nonsnp.py's table, which now
reports the >=90%-called (n_called>=208) filtered values by default -- so kMate's
side of this comparison is scored on sites called in >=90% of the 231-founder
panel, the fair basis against hapFIRE's fully-imputed zero-missingness panel.

Reuses ALL already-computed outputs -- no new simulation or estimator runs are
triggered by this script itself: hapFIRE's snp_frequency.txt (from the
greneNet_fair grid, benchmarks/speed_vs_hapfire/results/greneNet_fair/), kMate's
snp_R2/snp_RMSE (from poolsize_depth/scripts/score_snp_vs_nonsnp.py's output),
and the greneNet truth-AF matrix (from extract_greneNet_gt_matrix.py, run once).

Run with the `basic` env:
    /global/home/users/tbellg/miniforge3/envs/basic/bin/python \
        benchmarks/speed_vs_hapfire/scripts/score_af_hapfire_vs_kmate.py
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = "/global/scratch/users/tbellg/kmate"
BASE = f"{ROOT}/benchmarks/speed_vs_hapfire"
HF_SIMS = f"{BASE}/sims_greneNet"
HF_RES = f"{BASE}/results/greneNet_fair"
GT_MATRIX = f"{BASE}/work/greneNet_gt_matrix.npz"
OUT = f"{BASE}/results"
os.makedirs(OUT, exist_ok=True)

POOL_SIZES = [2, 5, 20, 50, 150, 231]
DEPTHS = [1, 10]
SEEDS = [42, 43, 44, 45, 46]
REL_EPS = 1e-6


def r2_rmse(est, truth):
    d = est - truth
    ss_tot = np.sum((truth - truth.mean()) ** 2)
    degenerate = ss_tot < (REL_EPS * truth.mean()) ** 2 * len(truth)
    r2 = np.nan if degenerate else 1 - np.sum(d ** 2) / ss_tot
    rmse = float(np.sqrt(np.mean(d ** 2)))
    return float(r2), rmse


print("loading greneNet genotype matrix...", flush=True)
gt = np.load(GT_MATRIX, allow_pickle=True)
gt_chrom, gt_pos, gt_founders, gt_dosage = gt["chrom"], gt["pos"], gt["founders"].astype(str), gt["dosage"]
founder_row = {f: i for i, f in enumerate(gt_founders)}
print(f"  {gt_dosage.shape[0]:,} SNPs x {gt_dosage.shape[1]} founders", flush=True)

kmate_snp = pd.read_csv(f"{ROOT}/benchmarks/speed_vs_hapfire/results/snp_vs_nonsnp_table.tsv", sep="\t")

rows = []
for n in POOL_SIZES:
    for cov in DEPTHS:
        for seed in SEEDS:
            hf_dir = f"{HF_RES}/n{n}_cov{cov}_s{seed}"
            sf_path = f"{hf_dir}/greneNet_hapfire_n{n}_cov{cov}_s{seed}_snp_frequency.txt"
            pw_path = f"{HF_SIMS}/cov{cov}_n{n}_g0_s{seed}_greneNet_chr1/pool_weights.tsv"
            if not (os.path.exists(sf_path) and os.path.exists(pw_path)):
                print(f"  [skip hapfire] n={n} cov={cov} s={seed}: missing output")
                hf_r2 = hf_rmse = hf_n = None
            else:
                pw = pd.read_csv(pw_path, sep="\t")
                pw["founder"] = pw["founder"].astype(str)
                weight = np.zeros(len(gt_founders))
                for f, w in zip(pw["founder"], pw["weight"]):
                    if f in founder_row:
                        weight[founder_row[f]] = w
                truth_af = gt_dosage @ weight

                sf = pd.read_csv(sf_path, sep="\t", header=None, names=["chrom", "pos", "freq"],
                                dtype={"chrom": str})
                truth_df = pd.DataFrame({"chrom": gt_chrom, "pos": gt_pos, "truth_af": truth_af})
                m = sf.merge(truth_df, on=["chrom", "pos"], how="inner")
                hf_r2, hf_rmse = r2_rmse(m["freq"].values, m["truth_af"].values)
                hf_n = len(m)

            ks = kmate_snp[(kmate_snp.N == n) & (kmate_snp.depth == cov) & (kmate_snp.seed == seed)]
            if len(ks):
                km_r2, km_rmse, km_n = ks["snp_R2"].iloc[0], ks["snp_RMSE"].iloc[0], int(ks["snp_n"].iloc[0])
            else:
                print(f"  [skip kmate] n={n} cov={cov} s={seed}: missing snp_vs_nonsnp row")
                km_r2 = km_rmse = km_n = None

            if hf_r2 is None and km_r2 is None:
                continue
            rows.append(dict(N=n, depth=cov, seed=seed,
                             kmate_af_R2=km_r2, kmate_af_RMSE=km_rmse, kmate_n_snps=km_n,
                             hapfire_af_R2=hf_r2, hapfire_af_RMSE=hf_rmse, hapfire_n_snps=hf_n))
            print(f"n={n:>3} cov={cov:>2}x s={seed}: "
                  f"kMate(arch3) R2={km_r2} RMSE={km_rmse} (n={km_n}) | "
                  f"hapFIRE(greneNet) R2={hf_r2} RMSE={hf_rmse} (n={hf_n})")

df = pd.DataFrame(rows)
table_path = f"{OUT}/hapfire_vs_kmate_af_table.tsv"
df.to_csv(table_path, sep="\t", index=False)
print(f"\nwrote {table_path} ({len(df)} rows)")

# ---------------------------------------------------------------------------
GREY = "#4d4d4d"
plt.rcParams.update({
    "text.color": GREY, "axes.labelcolor": GREY,
    "xtick.color": GREY, "ytick.color": GREY, "font.family": "sans-serif",
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.spines.left": False, "axes.spines.bottom": False,
    "axes.grid": True, "grid.color": "#dddddd", "grid.linewidth": 0.8,
    "axes.axisbelow": True, "xtick.bottom": False, "ytick.left": False,
})
TOOL_COLOR = {"kmate": "#54a24b", "hapfire": "#e45756"}
TOOL_LABEL = {"kmate": "kMate (arch3 SNPs, ≥90% called)", "hapfire": "hapFIRE (greneNet SNPs)"}
JITTER_RNG = np.random.RandomState(0)
PRESENT_N = [n for n in POOL_SIZES if n in df.N.unique()]


def box_panel(ax, sub, kmate_col, hapfire_col):
    x = np.arange(len(PRESENT_N))
    w = 0.32
    for ti, (tool, col) in enumerate([("kmate", kmate_col), ("hapfire", hapfire_col)]):
        vals = [sub[sub.N == n][col].dropna().values for n in PRESENT_N]
        pos = x + (ti - 0.5) * w
        ax.boxplot(vals, positions=pos, widths=w * 0.9, patch_artist=True,
                   showfliers=False, medianprops=dict(color="0.3", lw=1.3),
                   boxprops=dict(facecolor="none", edgecolor="0.3"),
                   whiskerprops=dict(color="0.3"), capprops=dict(color="0.3"))
        for p, v in zip(pos, vals):
            if len(v) == 0: continue
            jitter = JITTER_RNG.uniform(-w * 0.25, w * 0.25, size=len(v))
            ax.scatter(np.full(len(v), p) + jitter, v, color=TOOL_COLOR[tool],
                       s=20, alpha=0.85, edgecolor="white", linewidth=0.4, zorder=3)
    ax.set_xticks(x)
    ax.set_xticklabels([str(n) for n in PRESENT_N])


fig, axes = plt.subplots(2, len(DEPTHS), figsize=(4.6 * len(DEPTHS), 6.6), squeeze=False, sharex="col", sharey="row")
for ci, depth in enumerate(DEPTHS):
    sub = df[df.depth == depth]
    box_panel(axes[0][ci], sub, "kmate_af_R2", "hapfire_af_R2")
    axes[0][ci].text(0.5, 1.04, f"{depth}×", transform=axes[0][ci].transAxes,
                     fontsize=10, color="#999999", ha="center", va="bottom")
    if ci == 0: axes[0][ci].set_ylabel("R² (AF, own-panel SNPs)", fontsize=10)
    box_panel(axes[1][ci], sub, "kmate_af_RMSE", "hapfire_af_RMSE")
    axes[1][ci].set_xlabel("number of pooled founders")
    if ci == 0: axes[1][ci].set_ylabel("RMSE (AF, own-panel SNPs)", fontsize=10)

from matplotlib.lines import Line2D
handles = [Line2D([0], [0], marker="o", color="none", markerfacecolor=TOOL_COLOR[t],
                  markeredgecolor="white", markersize=7, label=TOOL_LABEL[t])
          for t in ("kmate", "hapfire")]
fig.tight_layout()
fig.legend(handles=handles, loc="lower center", ncol=2, fontsize=9,
           bbox_to_anchor=(0.5, -0.03), frameon=False)
out_path = f"{OUT}/hapfire_vs_kmate_af_accuracy.png"
fig.savefig(out_path, dpi=140, bbox_inches="tight")
plt.close(fig)
print("saved", out_path)
