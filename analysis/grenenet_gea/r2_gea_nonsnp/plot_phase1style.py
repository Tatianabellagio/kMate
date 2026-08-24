#!/usr/bin/env python
"""Render a candidate SV in the GrENE-net phase-1 figure style (3 panels):

  A  site-level frequency change (gen3 site-mean − p0) vs climate (bio1), with the
     regression line + Pearson r — the headline GEA readout at honest N≈20 sites.
  B  per-plot allele-frequency trajectories gen0(=shared p0)→1→2→3 in COLD gardens.
  C  same in WARM gardens.

Reusable: pass any SV by --chrom/--pos (+--p0 to disambiguate co-located records).
Outputs <out>.pdf (vector, for the poster) and <out>.png.

  PY=/global/home/users/tbellg/miniforge3/envs/basic/bin/python
  $PY analysis/grenenet_gea/r2_gea_nonsnp/plot_phase1style.py --chrom Chr4 --pos 15467237 \
      --p0 0.1215 --gene AT4G31980 --sym AT4G31980 \
      --note "80 bp deletion; rises in WARM gardens" --out-prefix candidate_AT4G31980_phase1style
"""
import argparse, glob, os, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import lib
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import cm
# keep text as editable text (not outlined glyphs) when the PDF is opened in Illustrator
matplotlib.rcParams["pdf.fonttype"] = 42      # TrueType (editable), not Type-3 outlines
matplotlib.rcParams["ps.fonttype"] = 42
matplotlib.rcParams["svg.fonttype"] = "none"

PM = f"{lib.GEA}/pool_matrices"; STORE = lib.AF_STORE
FIGDIR = f"{lib.GEA}/gea/figures"
SITES_NAMES = "/global/scratch/users/tbellg/gea_grene-net/key_files/sites_simple_names.csv"


def sv_column(chrom, pos, p0=None):
    """Map an SV (chrom,pos[,p0]) to its column in the *_nonsnp_af.npy matrices."""
    idx = np.load(f"{STORE}/index_nonsnp.npz")
    rl = idx["ref_len"].astype(np.int64); al = idx["alt_len"].astype(np.int64)
    size = np.abs(al - rl)
    nc = np.asarray(np.load(sorted(glob.glob(f"{STORE}/nc_nonsnp/*.npy"))[0]))
    mask = (size > 50) & (nc >= 150)
    nonsnp_idx = np.where(mask)[0]                       # col in matrices == nonsnp record index
    c = idx["chrom"].astype("U5")[mask]; p = idx["pos"][mask]
    p0all = np.load(f"{STORE}/p0_nonsnp.npy")[mask]
    cand = np.where((c == chrom) & (p == int(pos)))[0]
    if not len(cand):
        raise SystemExit(f"no SV at {chrom}:{pos}")
    j = cand[np.argmin(np.abs(p0all[cand] - p0))] if p0 is not None else cand[0]
    return nonsnp_idx[j], float(p0all[j]), int(rl[mask][j]), int(al[mask][j])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--chrom", required=True); ap.add_argument("--pos", type=int, required=True)
    ap.add_argument("--p0", type=float, default=None)
    ap.add_argument("--gene", default=""); ap.add_argument("--sym", default="")
    ap.add_argument("--note", default="")
    ap.add_argument("--split", type=float, default=None, help="bio1 cold/warm split (default=median site)")
    ap.add_argument("--avg-only", action="store_true", help="trajectory panels: site-mean lines only, no per-plot dots")
    ap.add_argument("--combined", action="store_true", help="single trajectory panel with ALL sites (cold+warm) colored by bio1")
    ap.add_argument("--label-sites", choices=["none", "city", "country", "city_country"], default="none",
                    help="panel A: site location names repelled off points with leader lines")
    ap.add_argument("--out-prefix", required=True)
    a = ap.parse_args()

    col, p0, rl, al = sv_column(a.chrom, a.pos, a.p0)
    kind = "insertion" if al > rl else "deletion"; svlen = abs(al - rl)

    # per-pool AF for this SV at gen1,2,3 + meta (site, plot, bio1, flowers)
    af = {0: None}; meta = {}
    for g in (1, 2, 3):
        P = np.load(f"{PM}/pool_gen{g}_nonsnp_af.npy", mmap_mode="r")
        af[g] = np.asarray(P[:, col], float)
        meta[g] = pd.read_csv(f"{PM}/pool_gen{g}_nonsnp.meta.csv")
        meta[g]["key"] = meta[g].site.astype(str) + "_" + meta[g]["plot"].astype(str)
    sites = np.array(sorted(meta[3].site.unique()))
    bio1 = meta[3].groupby("site").bio1.first()
    split = a.split if a.split is not None else float(np.median(bio1.loc[sites]))

    # Panel A: site-mean gen3 AF (flower-weighted) − p0, vs bio1
    rows = []
    for s in sites:
        m = meta[3][meta[3].site == s]
        f = m.total_flowers.to_numpy(float); v = af[3][m.index.to_numpy()]
        ok = np.isfinite(v) & (f > 0)
        if ok.sum() == 0:
            continue
        mean3 = np.sum(f[ok] * v[ok]) / np.sum(f[ok])
        rows.append((s, float(bio1.loc[s]), mean3 - p0))
    A = pd.DataFrame(rows, columns=["site", "bio1", "dp"])
    r = np.corrcoef(A.bio1, A.dp)[0, 1]

    # per-site mean trajectory (thick lines) + per-replicate-plot trajectories (thin faint
    # lines, in the background), gen0 = shared p0
    def panel_data(sites_sub):
        plot_lines = []        # (gens, afs, bio1) per replicate plot -> thin faint lines
        site_lines = []        # (gens, means, bio1)                  -> thick lines
        for s in sites_sub:
            sb = float(bio1.loc[s]); xs = [0]; ms = [p0]
            traj = {}                                       # key -> {gen: af}
            for g in (1, 2, 3):
                m = meta[g][meta[g].site == s]
                v = af[g][m.index.to_numpy()]; f = m.total_flowers.to_numpy(float)
                ok = np.isfinite(v) & (f > 0)
                if ok.sum():
                    xs.append(g); ms.append(np.sum(f[ok] * v[ok]) / np.sum(f[ok]))
                for kk, vi, fi in zip(m.key.to_numpy(), v, f):
                    if np.isfinite(vi) and fi > 0:
                        traj.setdefault(kk, {0: p0})[g] = float(vi)
            for d in traj.values():
                gg = sorted(d)
                if len(gg) > 1:
                    plot_lines.append(([g for g in gg], [d[g] for g in gg], sb))
            if len(xs) > 1:
                site_lines.append((xs, ms, sb))
        return plot_lines, site_lines

    # per-replicate-plot endpoint Δp (gen3 − p0) for the climate-gradient panel background
    def panelA_pts():
        pts = []; m = meta[3]
        for s in sites:
            mm = m[m.site == s]
            v = af[3][mm.index.to_numpy()]; f = mm.total_flowers.to_numpy(float)
            ok = np.isfinite(v) & (f > 0)
            for vi in v[ok]:
                pts.append((int(s), float(vi) - p0))
        return pts

    cold_sites = bio1.index[bio1 <= split]; warm_sites = bio1.index[bio1 > split]
    norm = plt.Normalize(float(bio1.loc[sites].min()), float(bio1.loc[sites].max()))
    cmap = cm.coolwarm

    if a.combined:
        fig, ax = plt.subplots(1, 2, figsize=(9.6, 6.3), gridspec_kw={"width_ratios": [1.1, 1]})
        traj_panels = [(sites, "Allele-frequency trajectories (all sites)")]
    else:
        fig, ax = plt.subplots(1, 3, figsize=(15, 3.6), gridspec_kw={"width_ratios": [1.15, 1, 1]})
        traj_panels = [(cold_sites, f"Colder sites (bio1≤{split:.1f})"),
                       (warm_sites, f"Warmer sites (bio1>{split:.1f})")]
    # A — climate-gradient panel: points at their TRUE temperature (x = bio1)
    if not a.avg_only:
        pa = panelA_pts()
        if pa:
            bx = np.array([float(bio1.loc[s]) for s, _ in pa]); by = np.array([d for _, d in pa])
            ax[0].scatter(bx, by, c=bx, cmap=cmap, norm=norm, s=24, alpha=.30, lw=0, zorder=1)
    sca = ax[0].scatter(A.bio1, A.dp, c=A.bio1, cmap=cmap, norm=norm, s=80, edgecolor="k", lw=.5, zorder=3)
    b, a0 = np.polyfit(A.bio1, A.dp, 1); xs = np.array([A.bio1.min(), A.bio1.max()])
    ax[0].plot(xs, a0 + b * xs, "k-", lw=1.3)
    ax[0].axhline(0, color="grey", lw=.6, ls=":")
    ax[0].set_ylabel("Frequency change (p₃ − p₀)")
    ax[0].set_title(f"Frequency change vs climate (r={r:.2f})")
    fig.colorbar(sca, ax=ax[0], label="bio1 (°C)")
    if a.label_sites == "none":
        ax[0].set_xlabel("Temperature (°C, bio1)")
    else:
        # location names EVENLY SPACED below the axis, each joined by a thin leader line to
        # its TRUE temperature tick (positions stay correct; spacing tweakable in Illustrator)
        sn = pd.read_csv(SITES_NAMES, encoding="utf-8-sig").set_index("site")
        loc = (sn.city.astype(str) + ", " + sn.country.astype(str)) if a.label_sites == "city_country" \
            else sn[a.label_sites].astype(str)
        o = A.sort_values("bio1").reset_index(drop=True); n = len(o)
        ylo, yhi = ax[0].get_ylim(); yr = yhi - ylo
        lab_x = np.linspace(float(A.bio1.min()), float(A.bio1.max()), n)   # evenly spaced names
        ax[0].tick_params(axis="x", labelbottom=False)   # keep temp ticks (-> vertical grid), hide numbers
        for i in range(n):
            s = int(o.site[i])
            ax[0].annotate(str(loc.get(s, s)), xy=(float(o.bio1[i]), ylo),
                           xytext=(lab_x[i], ylo - 0.05 * yr), textcoords="data",
                           rotation=90, ha="center", va="top", fontsize=7,
                           arrowprops=dict(arrowstyle="-", color="0.6", lw=.4),
                           annotation_clip=False)
        ax[0].set_xlabel("")   # names below serve as the axis (cold → warm, left → right)
    # trajectory panel(s) — thin faint per-replicate lines (unless --avg-only) under thick site-mean lines
    for axi, (sub, lab) in zip(ax[1:], traj_panels):
        plines, slines = panel_data(sub)
        if not a.avg_only:
            for xs_, ys_, sb in plines:
                axi.plot(xs_, ys_, "-", lw=.5, alpha=.18, color=cmap(norm(sb)), zorder=1)
        for xs_, ms_, sb in slines:
            axi.plot(xs_, ms_, "-o", lw=2.3, ms=3.5, alpha=1.0, color=cmap(norm(sb)), zorder=3)
        axi.axhline(p0, color="k", ls=":", lw=1, zorder=2)
        axi.set_xlabel("Generation"); axi.set_ylabel("Allele frequency")
        axi.set_xticks([0, 1, 2, 3]); axi.set_ylim(-0.02, 1.02); axi.set_title(lab)
        if a.combined:
            fig.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap), ax=axi, label="site bio1 (°C)")

    # clean "whitegrid" look: light grid behind data, no box (spines), gray ticks
    for axx in ax:
        axx.set_axisbelow(True)
        axx.grid(True, color="0.88", lw=.6, zorder=0)
        for sp in axx.spines.values():
            sp.set_visible(False)
        axx.tick_params(colors="0.5", labelcolor="0.3", length=0)

    gtag = f"{a.sym or a.gene}" + (f" ({a.gene})" if a.gene and a.sym and a.sym != a.gene else "")
    note = f"; {a.note}" if a.note else ""
    fig.suptitle(f"{a.chrom}:{a.pos:,}  ({svlen} bp {kind} at {gtag}; p₀={p0:.3f}{note})",
                 fontsize=11, y=1.02)
    fig.tight_layout()
    os.makedirs(FIGDIR, exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(f"{FIGDIR}/{a.out_prefix}.{ext}", bbox_inches="tight", dpi=150)
    print(f"wrote {FIGDIR}/{a.out_prefix}.pdf (+.png) | r={r:.2f}, p0={p0:.3f}, "
          f"{svlen}bp {kind}, {A.shape[0]} sites")


if __name__ == "__main__":
    main()
