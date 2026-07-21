#!/usr/bin/env python
"""Render MSR structure-null GEA figure (run in `basic` env).
Panel A: QQ of naive Spearman vs MSR-null p (both classes).
Panel B/C: Manhattan of MSR-null p (SNP / non-SNP) with GIF + Bonferroni/FDR lines.
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = "/global/scratch/users/tbellg/kmate/analysis/grenenet_gea/gea_newpanel/results/msr_kendall"
CLIM = "bio1"
CHR_ORDER = [f"Chr{i}" for i in range(1, 6)]
COL = {"snp": "#2c6fbb", "nonsnp": "#d1495b"}


def qq_xy(p):
    p = np.sort(p[np.isfinite(p)])
    n = p.size
    exp = -np.log10((np.arange(1, n + 1) - 0.5) / n)
    obs = -np.log10(np.clip(p, 1e-300, 1))
    # thin for plotting
    if n > 60000:
        idx = np.unique(np.linspace(0, n - 1, 60000).astype(int))
        return exp[idx], obs[idx]
    return exp, obs


def main():
    arr = {c: np.load(f"{OUT}/msr_arrays_{c}_{CLIM}.npz", allow_pickle=True)
           for c in ["snp", "nonsnp"]}
    man = {c: np.load(f"{OUT}/manhattan_{c}_{CLIM}.npz", allow_pickle=True)
           for c in ["snp", "nonsnp"]}

    fig = plt.figure(figsize=(13, 8))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.0], hspace=0.35, wspace=0.25)

    # ---- Panel A: QQ ----
    axq = fig.add_subplot(gs[0, 0])
    mx = 0
    for c in ["snp", "nonsnp"]:
        for pk, ls, lab in [("p_naive", ":", "naive"), ("p_msr", "-", "MSR")]:
            ex, ob = qq_xy(arr[c][pk])
            axq.plot(ex, ob, ls, color=COL[c], lw=1.4,
                     label=f"{c} {lab} (GIF {float(arr[c]['gif_'+('naive' if pk=='p_naive' else 'msr')]):.2f})")
            mx = max(mx, ob.max(), ex.max())
    axq.plot([0, mx], [0, mx], "k--", lw=0.8, alpha=0.6)
    axq.set_xlabel("expected -log10 p"); axq.set_ylabel("observed -log10 p")
    axq.set_title("QQ: naive Spearman (dotted) vs structure-preserving MSR null (solid)")
    axq.legend(fontsize=7, loc="upper left")

    # ---- Panel A2 (top-right): GIF bar summary ----
    axb = fig.add_subplot(gs[0, 1])
    labels = ["SNP naive", "SNP MSR", "nonSNP naive", "nonSNP MSR"]
    vals = [float(arr["snp"]["gif_naive"]), float(arr["snp"]["gif_msr"]),
            float(arr["nonsnp"]["gif_naive"]), float(arr["nonsnp"]["gif_msr"])]
    cols = [COL["snp"], COL["snp"], COL["nonsnp"], COL["nonsnp"]]
    hatch = ["", "//", "", "//"]
    bars = axb.bar(labels, vals, color=cols)
    for b, h in zip(bars, hatch):
        b.set_hatch(h)
    axb.axhline(1.0, color="k", ls="--", lw=0.8)
    axb.set_ylabel("genomic inflation factor (GIF)")
    axb.set_title("MSR null deflates inflation toward 1 (signal tail kept)")
    for b, v in zip(bars, vals):
        axb.text(b.get_x() + b.get_width() / 2, v + 0.03, f"{v:.2f}",
                 ha="center", fontsize=8)
    plt.setp(axb.get_xticklabels(), rotation=20, ha="right", fontsize=8)

    # ---- Panel B/C: Manhattan MSR ----
    for row_c, cls in [(gs[1, 0], "snp"), (gs[1, 1], "nonsnp")]:
        ax = fig.add_subplot(row_c)
        m = man[cls]
        ch = m["chrom"].astype(str); pos = m["pos"].astype(float); p = m["p_msr"].astype(float)
        offs = {}; cum = 0; ticks = []; ticklab = []
        for c in CHR_ORDER:
            sel = ch == c
            if not sel.any():
                continue
            mxp = pos[sel].max()
            offs[c] = cum
            ticks.append(cum + mxp / 2); ticklab.append(c.replace("Chr", ""))
            cum += mxp + 5e6
        for i, c in enumerate(CHR_ORDER):
            sel = ch == c
            if not sel.any():
                continue
            x = pos[sel] + offs[c]
            y = -np.log10(np.clip(p[sel], 1e-300, 1))
            ax.scatter(x, y, s=3, color=(COL[cls] if i % 2 == 0 else "#9aa7b5"),
                       alpha=0.5, linewidths=0)
        M = int(arr[cls]["p_msr"].size)
        bonf = -np.log10(0.05 / M)
        ax.axhline(bonf, color="k", ls="--", lw=0.8, label=f"Bonferroni 0.05 (M={M:,})")
        ax.axhline(-np.log10(1e-4), color="grey", ls=":", lw=0.9, label="p=1e-4 (MSR floor)")
        gif = float(arr[cls]["gif_msr"])
        np1e4 = int((arr[cls]["p_msr"] < 1e-4).sum())
        ax.set_title(f"{cls} MSR-null Manhattan (GIF {gif:.2f}; {np1e4} at p<1e-4; 0 FDR q<.05)")
        ax.set_xticks(ticks); ax.set_xticklabels(ticklab)
        ax.set_xlabel("chromosome"); ax.set_ylabel("-log10 p_MSR")
        ax.legend(fontsize=6, loc="upper right")

    fig.suptitle("Structure-preserving MSR-null rank climate-GEA (bio1, 31 sites): "
                 "SNP vs non-SNP", fontsize=12, y=0.98)
    png = f"{OUT}/msr_kendall_{CLIM}_figure.png"
    fig.savefig(png, dpi=140, bbox_inches="tight")
    print(f"-> {png}")


if __name__ == "__main__":
    main()
