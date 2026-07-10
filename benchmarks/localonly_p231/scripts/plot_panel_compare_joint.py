#!/usr/bin/env python3
"""Joint-density view of the panel & private-k-mer comparison.

One jointplot per series (p231 filt2inv / p80 filt2 / p80 unfiltered). The main
panel is a hexbin of nnz (log-x) vs per-block AF RMSE coloured by BLOCK COUNT
(log scale) — so you can see where the mass of blocks actually sits, which the
strip/median plots hide. Marginals: top = nnz supply histogram (this is the
thing that's otherwise invisible — how many blocks fall at each supply level,
arm vs centromere); right = RMSE histogram. Overlaid: running median, the
nnz=200 operating point, and the ~0.048 floor. Pools the same 3 scenarios."""
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.lines import Line2D

RES = Path("benchmarks/localonly_p231/results")
KEYS = ["chrom", "pos", "ref_len", "alt_len"]
SCEN = ["cov10_n50_g0_s42_hotspots", "cov10_n231_g1_s42_self97_hotspots",
        "cov10_n50_g3_s42_hotspots_dom500nr"]
EDGES = [1, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 100000]
CEN = (12_500_000, 17_500_000)
# reversed cmaps: LESS-dense hexes get the STRONGER colour (so the sparse tail
# doesn't get lost against the dense floor)
# (label, hue cmap, point/hist colour, floor_diag dir, pool suffix, variant, sims dir, truth fn)
CEN_COL = "0.45"  # grey for the centromere marginal
SERIES = [
    ("p231 filt2inv\n(231 mixed-tech, private dropped)", "Blues_r", "C0",
     "floor_diag", "_p231_chr1", "", "benchmarks/p231/sims", "recomb_truth_raw.tsv.gz"),
    ("p80 filt2\n(80 long-read, private dropped)", "Oranges_r", "C1",
     "floor_diag_p80", "_p80_chr1", "filt2", "benchmarks/p80/sims", "recomb_truth.tsv.gz"),
    ("p80 unfiltered\n(80 long-read, private KEPT)", "Greens_r", "C2",
     "floor_diag_p80", "_p80_chr1", "unfiltered", "benchmarks/p80/sims", "recomb_truth.tsv.gz"),
]


def per_block(fd, pool, variant, sims, truthfn, min_rec=20):
    base = Path("benchmarks/localonly_p231") / fd
    tag = f"{pool}_w10kb" + (f"_{variant}" if variant else "")
    bd, ef = base / f"{tag}.blockdiag.tsv", base / f"{tag}.recest.tsv"
    if not bd.exists() or not ef.exists():
        return None
    diag = pd.read_csv(bd, sep="\t"); est = pd.read_csv(ef, sep="\t")
    tr = pd.read_csv(Path(sims) / pool / truthfn, sep="\t").dropna(subset=["truth_af"])
    tr["occ"] = tr.groupby(KEYS).cumcount(); est["occ"] = est.groupby(KEYS).cumcount()
    m = tr.merge(est[KEYS + ["occ", "alt_freq", "block"]], on=KEYS + ["occ"], how="inner")
    m = m.merge(diag[["block", "nnz", "start"]], on="block", how="left")
    m = m[np.isfinite(m.alt_freq.values)].copy(); m["e2"] = (m.alt_freq - m.truth_af) ** 2
    g = pd.DataFrame({"rmse": np.sqrt(m.groupby("block").e2.mean()),
                      "n": m.groupby("block").alt_freq.size(),
                      "nnz": m.groupby("block").nnz.first(),
                      "start": m.groupby("block").start.first()})
    g = g[(g.n >= min_rec) & (g.nnz >= 1)].copy()
    g["cen"] = (g.start >= CEN[0]) & (g.start <= CEN[1])
    g["lnnz"] = np.log10(g.nnz.values)
    return g


XEDGES = np.log10([1, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 100000])
XTICKS = np.log10([1, 10, 100, 1000, 10000])
XTLAB = ["1", "10", "100", "1k", "10k"]
YMAX = 0.5
YBINS = np.linspace(0, YMAX, 41)


fig = plt.figure(figsize=(15, 5.6))
outer = fig.add_gridspec(1, 3, wspace=0.32, left=0.05, right=0.99, top=0.88, bottom=0.12)

for col, (lab, cmap, ccol, fd, suf, var, sims, tfn) in enumerate(SERIES):
    parts = [per_block(fd, s + suf, var, sims, tfn) for s in SCEN]
    parts = [p for p in parts if p is not None]
    if not parts:
        print(f"MISSING {lab}"); continue
    df = pd.concat(parts, ignore_index=True)
    print(f"{lab.splitlines()[0]:<16} n={len(df):>6} all-μ={df.rmse.mean():.3f} | "
          + " ".join(f"{t}:drop{100*(df.nnz<t).mean():4.1f}%/μ{df.loc[df.nnz>=t,'rmse'].mean():.3f}"
                     for t in (100, 200, 500)))

    inner = outer[0, col].subgridspec(2, 2, width_ratios=[4, 1], height_ratios=[1, 4],
                                      hspace=0.04, wspace=0.04)
    axm = fig.add_subplot(inner[1, 0])
    axt = fig.add_subplot(inner[0, 0], sharex=axm)
    axr = fig.add_subplot(inner[1, 1], sharey=axm)

    # main: hexbin density by block count (log colour)
    hb = axm.hexbin(df.lnnz, df.rmse.clip(upper=YMAX), gridsize=38, cmap=cmap,
                    bins="log", mincnt=1, extent=(0, np.log10(100000), 0, YMAX),
                    linewidths=0.15, edgecolors="face")
    # candidate min-kmers floors. VERTICAL dotted = % of scored blocks DROPPED (nnz<thr);
    # HORIZONTAL dotted (drawn rightward, over the KEPT region) = mean AF RMSE of the blocks
    # that survive the floor — the floor's payoff (median is tail-robust, so use the mean).
    GRY, XR = "0.35", np.log10(100000)
    for thr in (100, 200, 500):
        pct = 100.0 * (df.nnz < thr).mean()
        kept_mean = df.loc[df.nnz >= thr, "rmse"].mean()
        xt = np.log10(thr)
        axm.axvline(xt, color=GRY, ls=":", lw=0.9)
        axm.text(xt + 0.05, YMAX * 0.975, f"{thr}: {pct:.0f}% drop · mean RMSE {kept_mean:.3f}",
                 fontsize=6.3, rotation=90, va="top", ha="left", color=GRY,
                 bbox=dict(boxstyle="square,pad=0.1", fc="white", ec="none", alpha=0.7))
        axm.plot([xt, XR], [kept_mean, kept_mean], color=GRY, ls=":", lw=0.9)
    axm.set_xticks(XTICKS); axm.set_xticklabels(XTLAB)
    axm.set_xlim(0, np.log10(100000)); axm.set_ylim(0, YMAX)
    axm.set_xlabel("nonzero k-mers in block (nnz)")
    if col == 0:
        axm.set_ylabel("per-block AF RMSE")

    # top marginal: nnz supply, all blocks (series colour) + centromere (grey filled)
    axt.hist(df.lnnz, bins=XEDGES, color=ccol, alpha=0.55, edgecolor="white", lw=0.4)
    axt.hist(df.loc[df.cen, "lnnz"], bins=XEDGES, color=CEN_COL, alpha=0.75,
             edgecolor="white", lw=0.4)
    axt.set_yticks([]); axt.tick_params(labelbottom=False)
    axt.set_title(lab, fontsize=9.5, loc="left")
    axt.text(0.97, 0.78, f"n={len(df):,}\ncen={int(df.cen.sum()):,}",
             transform=axt.transAxes, ha="right", va="top", fontsize=7, color="dimgray")

    # right marginal: rmse, all blocks (series colour) + centromere (grey filled)
    axr.hist(df.rmse.clip(upper=YMAX), bins=YBINS, orientation="horizontal",
             color=ccol, alpha=0.55, edgecolor="white", lw=0.4)
    axr.hist(df.loc[df.cen, "rmse"].clip(upper=YMAX), bins=YBINS, orientation="horizontal",
             color=CEN_COL, alpha=0.75, edgecolor="white", lw=0.4)
    axr.set_xticks([]); axr.tick_params(labelleft=False)

from matplotlib.patches import Patch
fig.legend(handles=[Patch(facecolor=CEN_COL, alpha=0.75, label="centromere 12.5–17.5 Mb"),
                    Line2D([], [], color="0.35", ls=":", lw=0.9,
                           label="min-kmers floor: vert=% dropped, horiz=mean RMSE of kept blocks")],
           loc="upper right", fontsize=8, ncol=2, frameon=False, bbox_to_anchor=(0.99, 0.99))
fig.savefig(RES / "panel_compare_joint.png", dpi=140)
print(f"-> {RES/'panel_compare_joint.png'}  ({len(SERIES)} series)")
