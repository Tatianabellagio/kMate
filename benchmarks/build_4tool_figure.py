#!/usr/bin/env python3
"""Publication panel: RMSE (or R2/MAE) of kMate-global, kMate-block, hapFIRE, vg
across the p80 benchmark grid, from benchmark_table_4tool.tsv.

Layout: rows = var_class {SNP, SV}; one FIGURE per truth basis {fullcalled, allrec}
(fullcalled is the convention-free headline; allrec is the MAR-basis sensitivity).
x-axis = scenario (gen/mating/selection). Error bars = SD across seeds. 4 colored
bars/scenario. hapFIRE-SV is absent by design (no native SV) -> shown as a gap.

Usage: build_4tool_figure.py [--table ...] [--metric RMSE] [--out-prefix ...]
"""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

# tool/mode -> (label, color)
SERIES = [
    ("kMate", "global",  "kMate-global", "#2c6fbb"),
    ("kMate", "block",   "kMate-block",  "#7fb0e0"),
    ("hapFIRE", "default", "hapFIRE",     "#d1495b"),
    ("vg", "default",    "vg",           "#edae49"),
]
CLASSES = ["SNP", "SV"]


def scen_label(r):
    sel = {"balanced": "bal", "dom500": "dom", "dom500nr": "domNR"}.get(r.selection, r.selection)
    return f"g{int(r.generation)} {r.mating[:4]}\n{sel}"


def scen_order(r):
    selord = {"balanced": 0, "dom500": 1, "dom500nr": 2}
    return r.generation * 100 + (r.mating != "outcross") * 10 + selord.get(r.selection, 9)


def make_fig(t, metric, basis, out):
    sub0 = t[t.basis == basis].copy()
    sub0["scn"] = sub0.apply(scen_label, axis=1)
    sub0["_o"] = sub0.apply(scen_order, axis=1)
    g = (sub0.groupby(["var_class", "scn", "_o", "tool", "mode"])[metric]
              .agg(["mean", "std", "count"]).reset_index())
    scns = g.drop_duplicates("scn").sort_values("_o")["scn"].tolist()
    x = np.arange(len(scns)); w = 0.20

    fig, axes = plt.subplots(len(CLASSES), 1, figsize=(max(8, 1.1 * len(scns) + 3), 7.0),
                             squeeze=False, sharex=True)
    for ci, vc in enumerate(CLASSES):
        ax = axes[ci][0]
        gv = g[g.var_class == vc]
        for si, (tool, mode, lab, col) in enumerate(SERIES):
            s = (gv[(gv.tool == tool) & (gv["mode"] == mode)]
                 .set_index("scn").reindex(scns))
            ax.bar(x + (si - 1.5) * w, s["mean"].values, w, yerr=s["std"].values,
                   color=col, capsize=2, error_kw=dict(lw=0.8, alpha=0.7), label=lab)
        ax.set_ylabel(metric)
        ax.text(0.01, 0.93, vc, transform=ax.transAxes, fontsize=12, fontweight="bold",
                va="top", bbox=dict(fc="white", ec="none", alpha=0.7, pad=1.5))
        ax.grid(axis="y", ls=":", alpha=0.4); ax.margins(x=0.01)
    axes[-1][0].set_xticks(x); axes[-1][0].set_xticklabels(scns, fontsize=8)
    handles = [Patch(fc=c, label=l) for _, _, l, c in SERIES]
    fig.legend(handles=handles, loc="upper center", ncol=4, fontsize=10, frameon=False,
               bbox_to_anchor=(0.5, 1.0))
    better = "lower = better" if metric in ("RMSE", "MAE") else "higher = better"
    fig.suptitle(f"p80 accuracy: {metric} ({better}) — basis={basis}  "
                 f"(hapFIRE has no native SV; error bars = SD across seeds)",
                 y=1.03, fontsize=10, color="0.3")
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print(f"saved -> {out}  ({len(g)} cells, {len(scns)} scenarios)")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--table", default="benchmarks/benchmark_table_4tool.tsv")
    ap.add_argument("--metric", default="RMSE", choices=["RMSE", "MAE", "R2"])
    ap.add_argument("--out-prefix", default="benchmarks/benchmark_4tool")
    a = ap.parse_args()
    t = pd.read_csv(a.table, sep="\t")
    for basis in ["fullcalled", "allrec"]:
        make_fig(t, a.metric, basis, f"{a.out_prefix}_{a.metric}_{basis}.png")


if __name__ == "__main__":
    main()
