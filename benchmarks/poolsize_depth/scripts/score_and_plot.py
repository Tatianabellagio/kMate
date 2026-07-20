#!/usr/bin/env python3
"""Pool-size x depth accuracy sweep (Fig S13-style, hapFIRE paper): kMate accuracy
across pool size N and depth in {1x,10x}, 3 seeds each, on both panels (p80, p231).
N sweeps up to each panel's own founder count -- p80: {2,5,20,50,80}, p231:
{2,5,20,50,150} -- since N=150 on the 80-founder panel crosses into a different
sampling regime (with-replacement multinomial instead of gen0-no-replace) that
isn't a fair comparison point; see the PANELS comment below. Scores two metrics
per condition:
  - accession (founder-mixture h) accuracy: R2 and RMSE  [paper panel A]
  - allele-frequency accuracy: R2 and RMSE               [paper panels B/C,
    collapsed into one figure with depth as hue instead of a separate panel per depth]

Data: benchmarks/{p80,p231}/results/kmate_chrom_poolsize_depth/n{N}_cov{COV}_s{SEED}/
Produced by: benchmarks/p80/scripts/07j_run_kmate_poolsize_depth_p80.sh
             benchmarks/p231/scripts/07h_run_kmate_poolsize_depth_p231.sh

Run with the `basic` env:
    /global/home/users/tbellg/miniforge3/envs/basic/bin/python \
        benchmarks/poolsize_depth/scripts/score_and_plot.py
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = "/global/scratch/users/tbellg/kmate"
OUT = f"{ROOT}/benchmarks/poolsize_depth/results"
os.makedirs(OUT, exist_ok=True)

DEPTHS = [1, 10]
SEEDS = [42, 43, 44, 45, 46]
PANELS = {
    # p80's founder draw is `--gen0-no-replace`: at N=150 (> the 80-founder panel)
    # the draw switches to with-replacement multinomial and truth collapses to a
    # near-uniform 1-or-2-of-150 band; at N=80 (=panel size) it is *exactly*
    # uniform (1/80 for every founder). Both make truth variance ~0, so R2 (and
    # any metric normalized by sd(truth)) is either undefined or numerically
    # unstable there — a metric-normalization artifact, not estimator failure
    # (see poolsize_depth_table.tsv: af_R2 stays ~0.998 at N=150 throughout).
    # N=150 is excluded here as out-of-regime for p80 (it isn't a fair point of
    # comparison against p231, where N=150 < panel size 231 stays in-regime);
    # N=80 is kept as the (degenerate-truth) upper end of the in-regime sweep.
    "p80":  dict(res=f"{ROOT}/benchmarks/p80/results/kmate_chrom_poolsize_depth",
                 sims=f"{ROOT}/benchmarks/p80/sims",
                 sample_tpl="p80_chrom_psd_n{n}_cov{cov}_s{seed}",
                 truth_file="recomb_truth.tsv.gz",
                 pool_sizes=[2, 5, 20, 50, 80]),
    "p231": dict(res=f"{ROOT}/benchmarks/p231/results/kmate_chrom_poolsize_depth",
                 sims=f"{ROOT}/benchmarks/p231/sims",
                 sample_tpl="p231_chrom_psd_n{n}_cov{cov}_s{seed}",
                 truth_file="recomb_truth_raw.tsv.gz",
                 pool_sizes=[2, 5, 20, 50, 150]),
}
KEYS = ["chrom", "pos", "ref_len", "alt_len"]
# Below this fraction of truth's own mean, sd(truth) is floating-point noise
# around an exactly-flat vector, not real spread -- treat R2 as undefined
# rather than dividing by it (see p80 N=80: sd=1.7e-18 on a truth of 0.0125).
REL_EPS = 1e-6


def occ_key(df):
    return df.groupby(KEYS).cumcount()


def h_metrics(fo, h, truth_map):
    truth = np.array([truth_map.get(f, 0.0) for f in fo])
    d = h - truth
    ss_tot = np.sum((truth - truth.mean()) ** 2)
    degenerate = ss_tot < (REL_EPS * truth.mean()) ** 2 * len(truth)
    r2 = np.nan if degenerate else 1 - np.sum(d ** 2) / ss_tot
    rmse = float(np.sqrt(np.mean(d ** 2)))
    sd = truth.std()
    # RMSE/sd(truth): moves with R2 (both relative to truth's own spread),
    # unlike plain RMSE which shrinks in absolute terms as the mean founder
    # share ~1/N shrinks toward 0 -- see poolsize_depth_table.tsv N=150 p231:
    # RMSE keeps falling while R2 collapses. Matches the hapFIRE paper's own
    # "Normalized RMSE" (Fig S13) convention. Same degenerate guard as R2.
    rmse_norm = np.nan if degenerate else rmse / sd
    return dict(R2=float(r2), RMSE=rmse, RMSE_norm=float(rmse_norm))


def af_metrics(est_df, truth_df):
    truth_df = truth_df.dropna(subset=["truth_af"]).copy()
    truth_df["occ"] = occ_key(truth_df)
    est_df = est_df.copy()
    est_df["occ"] = occ_key(est_df)
    m = truth_df.merge(est_df[KEYS + ["occ", "alt_freq"]], on=KEYS + ["occ"], how="inner")
    d = (m["alt_freq"] - m["truth_af"]).values
    truth = m["truth_af"].values
    ss_tot = np.sum((truth - truth.mean()) ** 2)
    r2 = 1 - np.sum(d ** 2) / ss_tot if ss_tot > 0 else np.nan
    rmse = float(np.sqrt(np.mean(d ** 2)))
    return dict(R2=float(r2), RMSE=rmse, n=len(m))


rows = []
for panel, cfg in PANELS.items():
    for n in cfg["pool_sizes"]:
        for cov in DEPTHS:
            for seed in SEEDS:
                sample = cfg["sample_tpl"].format(n=n, cov=cov, seed=seed)
                d = f"{cfg['res']}/n{n}_cov{cov}_s{seed}"
                tsv, npz = f"{d}/{sample}.tsv", f"{d}/{sample}.h_per_chrom.npz"
                if not (os.path.exists(tsv) and os.path.exists(npz)):
                    print(f"  [skip] {panel} n={n} cov={cov} s={seed}: missing output")
                    continue
                sim = f"{cfg['sims']}/cov{cov}_n{n}_g0_s{seed}_hotspots_{panel}_chr1"
                # h accuracy
                hd = np.load(npz, allow_pickle=True)
                fo, h = hd["founders"].astype(str), hd["Chr1"].astype(float)
                tw = {}
                for r in pd.read_csv(f"{sim}/pool_weights.tsv", sep="\t").itertuples():
                    tw[str(r.founder)] = float(r.weight)
                hm = h_metrics(fo, h, tw)
                # AF accuracy
                est_df = pd.read_csv(tsv, sep="\t")
                truth_df = pd.read_csv(f"{sim}/{cfg['truth_file']}", sep="\t")
                afm = af_metrics(est_df, truth_df)
                rows.append(dict(panel=panel, N=n, depth=cov, seed=seed,
                                  h_R2=hm["R2"], h_RMSE=hm["RMSE"], h_RMSE_norm=hm["RMSE_norm"],
                                  af_R2=afm["R2"], af_RMSE=afm["RMSE"], af_n=afm["n"]))
                print(f"{panel} n={n:>3} cov={cov:>2}x s={seed}: "
                      f"h R2={hm['R2']:.4f} RMSE={hm['RMSE']:.4f} | "
                      f"AF R2={afm['R2']:.4f} RMSE={afm['RMSE']:.4f} (n={afm['n']:,})")

df = pd.DataFrame(rows)
table_path = f"{OUT}/poolsize_depth_table.tsv"
df.to_csv(table_path, sep="\t", index=False)
print(f"\nwrote {table_path} ({len(df)} rows)")

# ---------------------------------------------------------------------------
# Plotting: 2 figures (h-accuracy, AF-accuracy), each 2 rows (R2 / RMSE) x
# len(PANELS) columns, x = pool size (log scale, categorical positions),
# boxes colored by depth.
# ---------------------------------------------------------------------------
GREY = "#4d4d4d"
plt.rcParams.update({
    "text.color": GREY, "axes.labelcolor": GREY, "axes.titlecolor": GREY,
    "xtick.color": GREY, "ytick.color": GREY, "font.family": "sans-serif",
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.spines.left": False, "axes.spines.bottom": False,
    "axes.grid": True, "grid.color": "#dddddd", "grid.linewidth": 0.8,
    "axes.axisbelow": True, "xtick.bottom": False, "ytick.left": False,
})
DEPTH_COLOR = {1: "#dd8452", 10: "#4c78a8"}  # orange=1x, blue=10x


JITTER_RNG = np.random.RandomState(0)


def box_panel(ax, sub, ycol, panel_sizes):
    x = np.arange(len(panel_sizes))
    w = 0.32
    for di, depth in enumerate(DEPTHS):
        vals = [sub[(sub.N == n) & (sub.depth == depth)][ycol].dropna().values
                for n in panel_sizes]
        pos = x + (di - 0.5) * w
        ax.boxplot(vals, positions=pos, widths=w * 0.9, patch_artist=True,
                   showfliers=False, medianprops=dict(color="0.3", lw=1.3),
                   boxprops=dict(facecolor="none", edgecolor="0.3"),
                   whiskerprops=dict(color="0.3"), capprops=dict(color="0.3"))
        for p, v in zip(pos, vals):
            if len(v) == 0: continue
            jitter = JITTER_RNG.uniform(-w * 0.25, w * 0.25, size=len(v))
            ax.scatter(np.full(len(v), p) + jitter, v, color=DEPTH_COLOR[depth],
                       s=20, alpha=0.85, edgecolor="white", linewidth=0.4, zorder=3)
    ax.set_xticks(x)
    ax.set_xticklabels([str(n) for n in panel_sizes])


def make_figure(metric_r2, metric_rmse, rmse_label, out_path):
    fig, axes = plt.subplots(2, len(PANELS), figsize=(4.6 * len(PANELS), 6.6),
                             squeeze=False, sharex="col", sharey="row")
    for ci, panel in enumerate(PANELS):
        sub = df[df.panel == panel]
        panel_sizes = PANELS[panel]["pool_sizes"]
        box_panel(axes[0][ci], sub, metric_r2, panel_sizes)
        axes[0][ci].set_title(panel, fontsize=12, fontweight="bold", color=GREY)
        if ci == 0: axes[0][ci].set_ylabel("R²", fontsize=10)
        box_panel(axes[1][ci], sub, metric_rmse, panel_sizes)
        axes[1][ci].set_xlabel("number of pooled founders")
        if ci == 0: axes[1][ci].set_ylabel(rmse_label, fontsize=10)
    from matplotlib.lines import Line2D
    handles = [Line2D([0], [0], marker="o", color="none", markerfacecolor=DEPTH_COLOR[d],
                      markeredgecolor="white", markersize=7, label=f"{d}×") for d in DEPTHS]
    fig.tight_layout()
    fig.legend(handles=handles, loc="lower center", ncol=2, fontsize=9,
               bbox_to_anchor=(0.5, -0.03), frameon=False)
    fig.savefig(out_path, dpi=140, bbox_inches="tight")
    plt.close(fig)
    print("saved", out_path)


make_figure("h_R2", "h_RMSE_norm", "RMSE / sd(truth)",
            f"{OUT}/poolsize_depth_h_accuracy.png")
make_figure("af_R2", "af_RMSE", "RMSE",
            f"{OUT}/poolsize_depth_af_accuracy.png")
