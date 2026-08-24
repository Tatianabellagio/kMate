#!/usr/bin/env python
"""Across-garden allele-frequency dynamics of the CARK-block lead SV
(Chr2:13,127,635, 5.1 kb insertion, alt_len 5129) over the 3 experimental
generations, split into COLD vs WARM gardens and coloured by each garden's mean
annual temperature (bio1). Gen 0 = founding (SEEDMIX p0), shared by all gardens.

This is the per-SITE view behind the bio1 GEA hit: the SV is consistently higher
in warm gardens (positive climate slope) and lower/purged in cold gardens
(negative slope) -- the spatial climate signal LFMM detects (see README / the
combined locus figure). Data: per-generation AF matrices (gen_matrices/), which
stack the compact af_store per-sample vectors; site bio1 from ERA5 (lib.load_climate).

env: kmate.  Run on a compute node.
"""
from __future__ import annotations
import os, sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from scipy import stats

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
import lib

HERE = os.path.dirname(os.path.abspath(__file__))
GM = f"{lib.GEA}/gen_matrices"
CHROM, POS, RL, AL = "Chr2", 13_127_635, 1, 5129
KEY = f"{CHROM}:{POS}:{RL}:{AL}"


def sv_site_traj():
    """Per (site, generation) mean AF for the SV, + founding p0 as gen 0."""
    idx = np.load(f"{lib.AF_STORE}/index_nonsnp.npz", allow_pickle=True)
    ri = np.where((idx["chrom"] == CHROM) & (idx["pos"] == POS) &
                  (idx["ref_len"] == RL) & (idx["alt_len"] == AL))[0]
    if not len(ri):
        raise SystemExit(f"{KEY} not in nonsnp store")
    ri = int(ri[0])
    clim = lib.load_climate()["bio1"]

    rows = []
    for g in (1, 2, 3):
        M = np.load(f"{GM}/gen{g}_nonsnp_af.npy", mmap_mode="r")
        af = np.asarray(M[:, ri], dtype=np.float32)              # one column = all samples
        print(f"  gen{g}: read {len(af)} samples", flush=True)
        rm = pd.read_csv(f"{GM}/gen{g}.rowmeta.csv")
        rm["af"] = af
        rm = rm.dropna(subset=["af"])
        rm["site"] = rm["site"].astype(int)
        agg = (rm.groupby("site")
                 .agg(af=("af", "mean"), n=("af", "size"),
                      cov=("coverage", "mean")).reset_index())
        agg["generation"] = g
        rows.append(agg)
    traj = pd.concat(rows, ignore_index=True)
    traj["bio1"] = traj["site"].map(clim)
    traj = traj.dropna(subset=["bio1"])

    # gen 0 = founding p0 (same for every garden); pull from the fast group_means
    # cache (SEEDMIX mean) rather than re-reading the 8 SEEDMIX TSVs.
    z = np.load(f"{lib.GEA}/group_means.npz", allow_pickle=False)
    pm = np.where((z["chrom"] == CHROM) & (z["pos"] == POS) &
                  (z["ref_len"] == RL) & (z["alt_len"] == AL))[0]
    p0 = float(z["p0"][pm[0]]) if len(pm) else np.nan
    g0 = pd.DataFrame({"site": traj.site.unique()})
    g0["generation"] = 0
    g0["af"] = p0
    g0["bio1"] = g0["site"].map(clim)
    g0["n"] = np.nan
    g0["cov"] = np.nan
    return pd.concat([g0, traj], ignore_index=True), p0


def main():
    traj, p0 = sv_site_traj()
    traj.to_csv(f"{HERE}/cark_sv_site_traj.csv", index=False)
    nsite = traj.site.nunique()
    print(f"{KEY}: {nsite} gardens, p0={p0:.4f}, bio1 "
          f"{traj.bio1.min():.1f}..{traj.bio1.max():.1f} C")

    # cold vs warm split at the median garden bio1
    site_bio1 = traj.groupby("site")["bio1"].first()
    split = site_bio1.median()
    cold_sites = set(site_bio1[site_bio1 < split].index)
    warm_sites = set(site_bio1[site_bio1 >= split].index)

    vmin, vmax = traj.bio1.min(), traj.bio1.max()
    cmap = mcolors.LinearSegmentedColormap.from_list(
        "gardens", ["#2166ac", "#4a98c9", "#c9dCe8", "#f4c9b8", "#d6604d", "#8c1515"])
    norm = mcolors.Normalize(vmin, vmax)

    fig, axes = plt.subplots(1, 2, figsize=(11, 5.2), sharey=True)
    for ax, sites, title, c in [(axes[0], cold_sites, "Cold Gardens", "#2166ac"),
                                (axes[1], warm_sites, "Warm Gardens", "#b2182b")]:
        sub = traj[traj.site.isin(sites)]
        for site, s in sub.groupby("site"):
            s = s.sort_values("generation")
            col = cmap(norm(s.bio1.iloc[0]))
            ax.plot(s.generation, s.af, "-", color=col, lw=1.1, alpha=0.55, zorder=2)
            szs = np.where(np.isnan(s["n"]), 18, 12 + s["cov"].fillna(0).clip(0, 60))
            ax.scatter(s.generation, s.af, s=szs, color=col, alpha=0.85,
                       edgecolors="white", linewidths=0.3, zorder=3)
        # climate slope across all garden-gen points in the panel
        sl, ic, r, pv, se = stats.linregress(sub.generation, sub.af)
        # panel identity as an in-panel corner annotation (project convention: no titles)
        ax.text(0.5, 0.98, title, transform=ax.transAxes, va="top", ha="center",
                fontsize=15, fontweight="bold", color=c)
        ax.text(0.03, 0.88, f"Slope = {sl:.3f}\nP-value = {pv:.2e}",
                transform=ax.transAxes, va="top", ha="left", fontsize=9, color="#333")
        ax.set_xlabel("Generations", fontsize=13)
        ax.set_xticks([0, 1, 2, 3])
        ax.set_ylim(-0.02, 1.02)
        ax.grid(True, lw=0.3, c="0.9")
        ax.set_axisbelow(True)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
    axes[0].set_ylabel("Allele frequency", fontsize=13)

    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    cb = fig.colorbar(sm, ax=axes, fraction=0.025, pad=0.02)
    cb.set_label("Gardens Mean Annual Temperature (°C)", fontsize=10)
    cb.ax.tick_params(labelsize=8)

    out = f"{HERE}/cark_sv_garden_dynamics"
    fig.savefig(out + ".png", dpi=170, bbox_inches="tight")
    fig.savefig(out + ".pdf", bbox_inches="tight")
    print(f"wrote {out}.png/.pdf ; cold={len(cold_sites)} warm={len(warm_sites)} "
          f"gardens (split at bio1={split:.1f})")


if __name__ == "__main__":
    main()
