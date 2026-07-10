#!/usr/bin/env python3
"""Publication panel figure from benchmark_table.tsv: RMSE barplots across whichever
kMate modes are present for a given founder count (current: chrom/ld under --unit;
retained for provenance: global/block, the pre-2026-07-07 --block-mode labels),
faceted by variation class (rows) x panel (cols), x-axis = scenario.
Error bars = SD across seeds (g3 scenarios are multi-seed; single-seed -> no bar).
Lower RMSE = better.

Usage: build_benchmark_figure.py [--table ...] [--out ...] [--gens 3] [--metric RMSE]
"""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

MODE_COLOR = {"global": "#4C72B0", "block": "#DD8452",  # blue / orange (pre-2026-07-07 --block-mode)
              "chrom": "#55A868", "ld": "#8172B2"}       # green / purple (current --unit)
MODE_ORDER = ["global", "block", "chrom", "ld"]
CLASSES = ["SNP", "indel", "SV"]


def scen_label(r):
    sel = {"balanced": "bal", "dom500": "dom", "dom500nr": "domNR"}.get(r.selection, r.selection)
    return f"g{int(r.generation)} {r.mating[:4]}\n{sel}"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--table", default=str(Path(__file__).resolve().parent / "benchmark_table.tsv"))
    ap.add_argument("--out", default=str(Path(__file__).resolve().parent / "benchmark_panel_rmse.png"))
    ap.add_argument("--metric", default="RMSE", choices=["RMSE", "MAE", "R2"])
    ap.add_argument("--gens", default="0,1,3", help="generations to include (comma list)")
    ap.add_argument("--n-founders", type=int, default=50,
                    help="pin founder count so error bars reflect SEEDS only (not n)")
    a = ap.parse_args()

    t = pd.read_csv(a.table, sep="\t")
    gens = [int(x) for x in a.gens.split(",")]
    t = t[t.generation.isin(gens) & t.var_class.isin(CLASSES)
          & (t.n_founders == a.n_founders)].copy()
    t["scn"] = t.apply(scen_label, axis=1)
    # stable scenario order: by generation, mating (outcross first), selection
    selord = {"balanced": 0, "dom500": 1, "dom500nr": 2}
    t["_o"] = (t.generation * 100 + (t.mating != "outcross") * 10
               + t.selection.map(selord).fillna(9))

    # aggregate mean/sd across seeds
    g = (t.groupby(["panel", "var_class", "scn", "mode", "_o"])[a.metric]
           .agg(["mean", "std", "count"]).reset_index())

    modes = [m for m in MODE_ORDER if m in g["mode"].unique()]
    nm = max(len(modes), 1)
    panels = ["p80", "p231"]
    fig, axes = plt.subplots(len(CLASSES), len(panels), figsize=(7.0 * len(panels), 3.2 * len(CLASSES)),
                             squeeze=False, sharex="col")
    for ci, vc in enumerate(CLASSES):
        for pi, pan in enumerate(panels):
            ax = axes[ci][pi]
            sub = g[(g.var_class == vc) & (g.panel == pan)].sort_values("_o")
            scns = sub.drop_duplicates("scn").sort_values("_o")["scn"].tolist()
            x = np.arange(len(scns)); w = 0.8 / nm
            for mi, mode in enumerate(modes):
                s = sub[sub["mode"] == mode].set_index("scn").reindex(scns)
                ax.bar(x + (mi - (nm - 1) / 2) * w, s["mean"].values, w, yerr=s["std"].values,
                       color=MODE_COLOR[mode], capsize=2.5, error_kw=dict(lw=1, alpha=0.7),
                       label=mode)
            ax.set_xticks(x); ax.set_xticklabels(scns, fontsize=7)
            ax.set_ylabel(a.metric if pi == 0 else "")
            if ci == 0:
                ax.set_title(pan, fontsize=12, fontweight="bold")
            ax.text(0.015, 0.92, vc, transform=ax.transAxes, fontsize=10, fontweight="bold",
                    va="top", bbox=dict(fc="white", ec="none", alpha=0.7, pad=1.5))
            ax.grid(axis="y", ls=":", alpha=0.4)
            ax.margins(x=0.01)
    handles = [Patch(fc=MODE_COLOR[m], label=f"kMate {m}") for m in modes]
    fig.legend(handles=handles, loc="upper center", ncol=nm, fontsize=11,
               frameon=False, bbox_to_anchor=(0.5, 1.0))
    better = "lower = better" if a.metric in ("RMSE", "MAE") else "higher = better"
    fig.suptitle(f"kMate accuracy ({a.metric}, {better}); error bars = SD across seeds",
                 y=1.03, fontsize=10, color="0.3")
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(a.out, dpi=150, bbox_inches="tight")
    print(f"saved -> {a.out}  ({len(g)} scenario×mode×class cells)")


if __name__ == "__main__":
    main()
