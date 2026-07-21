#!/usr/bin/env python3
"""kMate-only AF accuracy, SNP vs non-SNP (indel/SV), same p231 N x depth grid
as hapfire_vs_kmate_af_accuracy.png -- but hue = variant type instead of tool
(kMate has no indel/SV comparator in hapFIRE, so this is kMate-internal).
Reuses already-computed kMate tsv outputs + recomb_truth_raw.tsv.gz, no new
compute needed.

The PRIMARY snp_R2/snp_RMSE/nonsnp_* columns are computed on records called in
>=90% of the 231-founder panel (n_called>=208), for a fair comparison against
hapFIRE's fully-imputed zero-missingness greneNet panel; unfiltered values are
kept in *_allrecords columns. This table feeds score_af_hapfire_vs_kmate.py, so
the kMate side of the cross-tool AF comparison inherits the same filter.

Run with the `basic` env:
    /global/home/users/tbellg/miniforge3/envs/basic/bin/python \
        benchmarks/poolsize_depth/scripts/score_snp_vs_nonsnp.py
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = "/global/scratch/users/tbellg/kmate"
P231 = f"{ROOT}/benchmarks/p231"
OUT = f"{ROOT}/benchmarks/speed_vs_hapfire/results"
os.makedirs(OUT, exist_ok=True)

POOL_SIZES = [2, 5, 20, 50, 150, 231]
DEPTHS = [1, 10]
SEEDS = [42, 43, 44, 45, 46]
REL_EPS = 1e-6
# Missingness filter: keep kMate records called in >=90% of the 231-founder panel
# (n_called >= 208, i.e. <=10% missing). hapFIRE's greneNet panel is fully imputed
# (zero missingness), so scoring kMate on comparably-complete sites is the fair
# apples-to-apples basis; this is the PRIMARY reported number. All-records values
# are retained in the *_allrecords columns for provenance.
FULL_PANEL = 231
MISS_THRESH = int(np.ceil(0.9 * FULL_PANEL))  # 208


def r2_rmse(est, truth):
    d = est - truth
    ss_tot = np.sum((truth - truth.mean()) ** 2)
    degenerate = ss_tot < (REL_EPS * truth.mean()) ** 2 * len(truth)
    r2 = np.nan if degenerate else 1 - np.sum(d ** 2) / ss_tot
    rmse = float(np.sqrt(np.mean(d ** 2)))
    return float(r2), rmse


rows = []
for n in POOL_SIZES:
    for cov in DEPTHS:
        for seed in SEEDS:
            sim = f"{P231}/sims/cov{cov}_n{n}_g0_s{seed}_hotspots_p231_chr1"
            kmate_dir = f"{P231}/results/kmate_chrom_poolsize_depth/n{n}_cov{cov}_s{seed}"
            kmate_tsv = f"{kmate_dir}/p231_chrom_psd_n{n}_cov{cov}_s{seed}.tsv"
            truth_path = f"{sim}/recomb_truth_raw.tsv.gz"
            if not (os.path.exists(kmate_tsv) and os.path.exists(truth_path)):
                print(f"  [skip] n={n} cov={cov} s={seed}: missing output")
                continue

            KEYS = ["chrom", "pos", "ref_len", "alt_len"]
            est = pd.read_csv(kmate_tsv, sep="\t")
            truth = pd.read_csv(truth_path, sep="\t").dropna(subset=["truth_af"])
            # multi-allelic sites repeat (chrom,pos,ref_len,alt_len); cumcount as an
            # occurrence tie-breaker avoids a many-to-many cartesian blowup on merge.
            # NOTE: occ pairs by record order, not allele sequence (neither file carries
            # the allele), so same-length multiallelic sites are matched positionally.
            # Audited immaterial: max aggregate RMSE impact ~3e-6 (order matches at all
            # N=2 sites, differs at only ~2.6% of multiallelic sites by N=231).
            truth = truth.copy(); truth["occ"] = truth.groupby(KEYS).cumcount()
            est = est.copy(); est["occ"] = est.groupby(KEYS).cumcount()
            m = truth.merge(est[KEYS + ["occ", "alt_freq", "n_called"]], on=KEYS + ["occ"], how="inner")

            is_snp = (m.ref_len == 1) & (m.alt_len == 1)
            keep = m["n_called"] >= MISS_THRESH  # missingness filter (fair vs hapFIRE)
            snp_k = m[is_snp & keep]; nonsnp_k = m[(~is_snp) & keep]
            snp_a = m[is_snp]; nonsnp_a = m[~is_snp]

            # PRIMARY: >=90%-called (fair comparison basis)
            snp_r2, snp_rmse = r2_rmse(snp_k["alt_freq"].values, snp_k["truth_af"].values)
            nonsnp_r2, nonsnp_rmse = r2_rmse(nonsnp_k["alt_freq"].values, nonsnp_k["truth_af"].values)
            # PROVENANCE: all records, unfiltered
            snp_r2_a, snp_rmse_a = r2_rmse(snp_a["alt_freq"].values, snp_a["truth_af"].values)
            nonsnp_r2_a, nonsnp_rmse_a = r2_rmse(nonsnp_a["alt_freq"].values, nonsnp_a["truth_af"].values)

            rows.append(dict(N=n, depth=cov, seed=seed,
                             snp_R2=snp_r2, snp_RMSE=snp_rmse, snp_n=len(snp_k),
                             nonsnp_R2=nonsnp_r2, nonsnp_RMSE=nonsnp_rmse, nonsnp_n=len(nonsnp_k),
                             snp_R2_allrecords=snp_r2_a, snp_RMSE_allrecords=snp_rmse_a, snp_n_allrecords=len(snp_a),
                             nonsnp_R2_allrecords=nonsnp_r2_a, nonsnp_RMSE_allrecords=nonsnp_rmse_a, nonsnp_n_allrecords=len(nonsnp_a)))
            print(f"n={n:>3} cov={cov:>2}x s={seed}: [>=90% called] "
                  f"SNP R2={snp_r2:.4f} RMSE={snp_rmse:.4f} (n={len(snp_k):,}) | "
                  f"non-SNP R2={nonsnp_r2:.4f} RMSE={nonsnp_rmse:.4f} (n={len(nonsnp_k):,})")

df = pd.DataFrame(rows)
table_path = f"{OUT}/snp_vs_nonsnp_table.tsv"
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
TYPE_COLOR = {"snp": "#4c78a8", "nonsnp": "#f58518"}
TYPE_LABEL = {"snp": "SNP", "nonsnp": "indel/SV"}
JITTER_RNG = np.random.RandomState(0)


def box_panel(ax, sub, ycol_snp, ycol_nonsnp):
    x = np.arange(len(POOL_SIZES))
    w = 0.32
    for ti, (typ, col) in enumerate([("snp", ycol_snp), ("nonsnp", ycol_nonsnp)]):
        vals = [sub[sub.N == n][col].dropna().values for n in POOL_SIZES]
        pos = x + (ti - 0.5) * w
        ax.boxplot(vals, positions=pos, widths=w * 0.9, patch_artist=True,
                   showfliers=False, medianprops=dict(color="0.3", lw=1.3),
                   boxprops=dict(facecolor="none", edgecolor="0.3"),
                   whiskerprops=dict(color="0.3"), capprops=dict(color="0.3"))
        for p, v in zip(pos, vals):
            if len(v) == 0: continue
            jitter = JITTER_RNG.uniform(-w * 0.25, w * 0.25, size=len(v))
            ax.scatter(np.full(len(v), p) + jitter, v, color=TYPE_COLOR[typ],
                       s=20, alpha=0.85, edgecolor="white", linewidth=0.4, zorder=3)
    ax.set_xticks(x)
    ax.set_xticklabels([str(n) for n in POOL_SIZES])


fig, axes = plt.subplots(2, len(DEPTHS), figsize=(4.6 * len(DEPTHS), 6.6), squeeze=False, sharex="col", sharey="row")
for ci, depth in enumerate(DEPTHS):
    sub = df[df.depth == depth]
    box_panel(axes[0][ci], sub, "snp_R2", "nonsnp_R2")
    axes[0][ci].text(0.5, 1.04, f"{depth}×", transform=axes[0][ci].transAxes,
                     fontsize=10, color="#999999", ha="center", va="bottom")
    if ci == 0: axes[0][ci].set_ylabel("R² (AF, ≥90% called)", fontsize=10)
    box_panel(axes[1][ci], sub, "snp_RMSE", "nonsnp_RMSE")
    axes[1][ci].set_xlabel("number of pooled founders")
    if ci == 0: axes[1][ci].set_ylabel("RMSE (AF, ≥90% called)", fontsize=10)

from matplotlib.lines import Line2D
handles = [Line2D([0], [0], marker="o", color="none", markerfacecolor=TYPE_COLOR[t],
                  markeredgecolor="white", markersize=7, label=TYPE_LABEL[t])
          for t in ("snp", "nonsnp")]
fig.tight_layout()
fig.legend(handles=handles, loc="lower center", ncol=2, fontsize=9,
           bbox_to_anchor=(0.5, -0.03), frameon=False)
out_path = f"{OUT}/kmate_snp_vs_nonsnp_af_accuracy.png"
fig.savefig(out_path, dpi=140, bbox_inches="tight")
plt.close(fig)
print("saved", out_path)
