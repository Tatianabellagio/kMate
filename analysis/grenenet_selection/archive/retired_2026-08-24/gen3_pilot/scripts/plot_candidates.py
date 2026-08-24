#!/usr/bin/env python
"""Poster figures for candidate climate-adaptive SVs (mirrors the phase-1 panels).

For each candidate non-SNP record, two panels:
  B  Frequency change (p3 - p0) vs site temperature (bio1) across the gen-3 pools,
     points colored by temperature, with the regression line — the climate
     gradient (the Kendall-tau signal made visual).
  C  Frequency trajectory over years 0->3 (year 0 = founding p0), one faint line
     per plot, split into Cold vs Warm gardens — the up-in-warm / flat-in-cold
     selection signature.

Run in the `plotting` conda env (matplotlib). Reads the pool matrices, p0, and
the candidate list.

Usage:
  PLOTPY=/global/home/users/tbellg/miniforge3/envs/plotting/bin/python
  $PLOTPY analysis/grenenet_selection/r2_gea_nonsnp/plot_candidates.py --candidates analysis/grenenet_selection/r2_gea_nonsnp/results/gea/poster_candidates.csv --n 6
"""
from __future__ import annotations
import argparse, os, sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import cm, colors
from scipy.stats import linregress, kendalltau
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import lib

POOLDIR = f"{lib.GEA}/common/results/pool_matrices"
NORM = colors.Normalize(vmin=5, vmax=22)       # temperature color scale (figure)
CMAP = cm.get_cmap("RdBu_r")


def _load(kind="nonsnp"):
    p0 = np.load(f"{lib.AF_STORE}/p0_{kind}.npy")
    gens = {}
    for g in (1, 2, 3):
        af = np.load(f"{POOLDIR}/pool_gen{g}_{kind}_af.npy", mmap_mode="r")
        mt = pd.read_csv(f"{POOLDIR}/pool_gen{g}_{kind}.meta.csv")
        gens[g] = (af, mt)
    return p0, gens


def plot_candidate(ax_b, ax_cold, ax_warm, rec, p0, gens, title):
    af3, mt3 = gens[3]
    temp = mt3.bio1.to_numpy()
    dp = af3[:, rec].astype(float) - p0[rec]
    # ---- Panel B: Δp vs temperature ----
    ax_b.axhline(0, color="0.8", lw=.8, zorder=0)
    ax_b.scatter(temp, dp * 100, c=temp, cmap=CMAP, norm=NORM, s=26,
                 edgecolor="0.3", linewidth=.3, zorder=3)
    m = np.isfinite(dp)
    lr = linregress(temp[m], dp[m])
    xs = np.linspace(temp.min(), temp.max(), 50)
    ax_b.plot(xs, (lr.intercept + lr.slope * xs) * 100, color="0.25", lw=1.6)
    tau, pv = kendalltau(temp[m], dp[m])
    ax_b.set_title(f"{title}\nτ={tau:.2f}  p={pv:.1e}", fontsize=8)
    ax_b.set_xlabel("Temperature (°C)"); ax_b.set_ylabel("Frequency change (p₃−p₀)")
    # ---- Panel C: trajectories, year 0..3, cold vs warm ----
    tmed = np.median(np.unique(temp))
    # per-plot freq over years: year0 = p0 (all plots), years 1-3 from pool matrices
    traj = {}   # (site,plot) -> {year: freq, temp}
    for g in (1, 2, 3):
        af, mt = gens[g]
        col = af[:, rec].astype(float)
        for i, r in mt.iterrows():
            key = (int(r.site), int(r["plot"]))
            traj.setdefault(key, {"temp": r.bio1, "yrs": {0: p0[rec] * 100}})
            traj[key]["yrs"][g] = col[i] * 100
    for ax, warm in [(ax_cold, False), (ax_warm, True)]:
        ax.axhline(p0[rec] * 100, color="0.85", lw=.8, zorder=0)
        for key, d in traj.items():
            if (d["temp"] >= tmed) != warm:
                continue
            ys = sorted(d["yrs"].items())
            xx = [y for y, _ in ys]; yy = [v for _, v in ys]
            ax.plot(xx, yy, "-", color=CMAP(NORM(d["temp"])), alpha=.45, lw=1,
                    marker="o", ms=2.5)
        ax.set_xticks([0, 1, 2, 3]); ax.set_xlabel("Years")
        ax.set_title("Warm gardens" if warm else "Cold gardens",
                     color=CMAP(NORM(20 if warm else 8)), fontsize=9)
    ax_cold.set_ylabel("Frequency (pₜ)  %")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--candidates", default=f"{lib.GEA}/r2_gea_nonsnp/results/gea/poster_candidates.csv")
    ap.add_argument("--n", type=int, default=6)
    ap.add_argument("--kind", default="nonsnp")
    ap.add_argument("--out", default=f"{lib.GEA}/r2_gea_nonsnp/results/gea/figures")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)
    cands = pd.read_csv(args.candidates).head(args.n)
    p0, gens = _load(args.kind)
    for _, c in cands.iterrows():
        rec = int(c.rec_index)
        fig = plt.figure(figsize=(9, 3))
        gs = fig.add_gridspec(1, 4, width_ratios=[1.5, 0.05, 1, 1], wspace=.45)
        ax_b = fig.add_subplot(gs[0]); ax_cold = fig.add_subplot(gs[2])
        ax_warm = fig.add_subplot(gs[3], sharey=ax_cold)
        title = (f"{c.chrom}:{int(c.pos)}  {int(c.sv_size)}bp SV  "
                 f"(p₀={c.p0:.2f})")
        plot_candidate(ax_b, ax_cold, ax_warm, rec, p0, gens, title)
        sm = cm.ScalarMappable(norm=NORM, cmap=CMAP); sm.set_array([])
        fig.colorbar(sm, cax=fig.add_subplot(gs[1]), label="Temp (°C)")
        tag = f"{c.chrom}_{int(c.pos)}_{int(c.sv_size)}bp"
        fig.savefig(f"{args.out}/cand_{tag}.png", dpi=160, bbox_inches="tight")
        plt.close(fig)
        print(f"  wrote cand_{tag}.png  (τ={c.tau:.2f})")
    print(f"-> {args.out}/")


if __name__ == "__main__":
    main()
