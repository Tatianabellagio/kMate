#!/usr/bin/env python
"""Per-garden allele-frequency trajectory of ONE variant, as small multiples.

One small panel per garden, ordered cold -> warm by site bio1 in reading order, each
showing the variant's frequency across the 3 experimental generations from the shared
founding p0 (gen 0 = SEEDMIX mean, identical in every garden).

Generalises `plot_sv_garden_dynamics.py`, which is hardwired to the CARK-block SV and
collapses everything into two cold/warm panels. Two panels hide which gardens carry
the signal; that is exactly the question for a per-garden GWAS hit, so here each
garden keeps its own axes.

What each panel shows:
  * faint dots      one per plot replicate (the actual sequenced pools)
  * bold line       the flower-weighted-free site mean per generation
  * dashed grey     the founding p0, so rise/fall is read against the start
  * coloured        by that garden's mean annual temperature (bio1), blue -> red
  * red panel frame the garden is Bonferroni-significant for this variant in the
                    per-garden GEMMA scan named by `sigscan` (default: the widest)

Per project convention there are NO titles: garden identity is an in-panel corner
annotation.

Usage (JSON config as argv[1]):
  PY plot_variant_garden_grid.py '{"chrom":"Chr5","pos":1146373,"ref_len":95,
     "alt_len":1,"sym":"CML50","label":"94 bp intronic deletion"}'

Optional keys: `ncol` (default 6), `sigscan` ("sv"|"nonsnp"|"snp"|"auto", default
auto), `out` (basename override).

env: kmate.  Run on a compute node.
"""
from __future__ import annotations
import os
import sys
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(HERE, "..", "..")))
import lib                                                    # noqa: E402

GM = f"{lib.GEA}/common/results/gen_matrices"
GWASD = f"{lib.GEA}/r3_persite_gwas/results/gemma_gwas"
# per-allele significant gardens for the non-SNP scans (genes/convergence/resolve_gwas_alleles.py)
RESOLVED = f"{lib.GEA}/genes/convergence/results/gwas_hits_allele_resolved.csv"
OUT = f"{HERE}/results/plots"
GENS = (1, 2, 3)
# per-class genome-wide Bonferroni from each GEMMA scan's own marker count
BONF_GWAS = {"snp": 7.545, "nonsnp": 7.021, "sv": 5.408}
# Cool->warm WITHOUT a near-white midpoint. The previous ramp passed through #c9dce8 /
# #f4c9b8, which are both ~90% luminance, so the 13 of 30 gardens sitting in the middle
# of the bio1 range (6.0-21.6 C, most gardens 11-15 C) rendered nearly invisible on a
# white background. Routing cool->warm through saturated teal and amber instead of pale
# blue and pale peach keeps the temperature reading and holds luminance roughly constant.
CMAP = mcolors.LinearSegmentedColormap.from_list(
    "gardens", ["#123f6d", "#2166ac", "#3d90b8", "#5aa9a0",
                "#c08a2e", "#d1603c", "#b02418", "#6e0f0f"])


def store_row(cfg) -> int:
    """Row index of the variant in the compact non-SNP AF store.

    ⚠ `(chrom, pos, ref_len, alt_len)` is NOT unique -- two ALT sequences of the same
    length at one position collide, affecting 14.1% of non-SNP records. Taking `[0]` here
    while a caller keys the same variant with a dict (which keeps the LAST collider) makes
    the figure and the caller's table describe different variants. Measured on SEC10a:
    p0 0.000716 here vs 0.1138 in the caller.

    So an explicit `store_row` in the config wins, and a silent guess is never made: an
    ambiguous key without `store_row` is a hard error, not a coin flip.
    """
    if "store_row" in cfg:
        return int(cfg["store_row"])
    idx = np.load(f"{lib.AF_STORE}/index_nonsnp.npz", allow_pickle=True)
    m = ((idx["chrom"] == cfg["chrom"]) & (idx["pos"] == cfg["pos"])
         & (idx["ref_len"] == cfg["ref_len"]) & (idx["alt_len"] == cfg["alt_len"]))
    r = np.where(m)[0]
    if not len(r):
        raise SystemExit(f"{cfg['chrom']}:{cfg['pos']} "
                         f"({cfg['ref_len']}/{cfg['alt_len']}) not in the nonsnp store")
    if len(r) > 1:
        raise SystemExit(
            f"{cfg['chrom']}:{cfg['pos']} ({cfg['ref_len']}/{cfg['alt_len']}) matches "
            f"{len(r)} af_store records {list(r)} -- ambiguous key. Pass an explicit "
            f'"store_row" in the config to say which one.')
    return int(r[0])


def founding_p0(cfg, row: int | None = None) -> float:
    """Founding frequency for the SAME record the AF column came from.

    ⚠ Must be resolved positionally when `store_row` is given. This function used to key on
    (chrom,pos,ref_len,alt_len) and take `[0]`, which is a different record from the one
    `store_row()` returns whenever the key is ambiguous (14.1% of non-SNP records). The
    result was a figure plotting the RIGHT allele-frequency column against the WRONG
    founding line -- two different variants in one panel. Measured on SEC10a: p0 0.000716
    drawn under the AF trajectory of the record whose p0 is 0.1138.

    Verified once: the non-SNP subset of group_means.npz is row-for-row identical to
    index_nonsnp.npz on all four fields, so af_store row i is group_means row
    where(non_snp)[0][i].
    """
    z = np.load(f"{lib.GEA}/common/results/group_means.npz", allow_pickle=False)
    if row is not None:
        ns = np.where((z["ref_len"] != 1) | (z["alt_len"] != 1))[0]
        j = ns[row]
        assert (str(z["chrom"][j]) == cfg["chrom"] and int(z["pos"][j]) == int(cfg["pos"])
                and int(z["ref_len"][j]) == int(cfg["ref_len"])
                and int(z["alt_len"][j]) == int(cfg["alt_len"])), \
            "group_means/af_store row alignment broken -- re-verify before use"
        return float(z["p0"][j])
    m = np.where((z["chrom"] == cfg["chrom"]) & (z["pos"] == cfg["pos"])
                 & (z["ref_len"] == cfg["ref_len"])
                 & (z["alt_len"] == cfg["alt_len"]))[0]
    if len(m) > 1:
        raise SystemExit(f"{cfg['chrom']}:{cfg['pos']} matches {len(m)} group_means records "
                         f"-- ambiguous key; pass an explicit \"store_row\".")
    return float(z["p0"][m[0]]) if len(m) else np.nan


def per_plot_af(cfg) -> pd.DataFrame:
    """One row per POOL (site x plot x generation), flower-weighted over timepoints.

    `gen_matrices` is the RAW per-SAMPLE layer: one row per sequencing timepoint, so a
    plot sampled four times contributes four rows. Plotting those directly drew 70 dots
    for garden 13 gen2, which has 11 plots, and made the bold line a sample-mean that
    over-weights the most-sampled plots. The analysis unit everywhere else in this tree
    is the pool, with timepoints merged flower-weighted:

        p_pool = sum_t (flowers_t * p_t) / sum_t flowers_t

    which is what `build_selection_trait.py` and `build_pool_matrix.py` both do. Weight
    fallback (non-finite or <=0 -> 1.0) matches build_selection_trait.py exactly.
    """
    ri = store_row(cfg)
    frames = []
    for g in GENS:
        M = np.load(f"{GM}/gen{g}_nonsnp_af.npy", mmap_mode="r")
        rm = pd.read_csv(f"{GM}/gen{g}.rowmeta.csv")
        rm["af"] = np.asarray(M[:, ri], dtype=np.float32)
        rm = rm.dropna(subset=["af"]).copy()
        w = rm["flowerscollected"].to_numpy(float)
        rm["w"] = np.where(np.isfinite(w) & (w > 0), w, 1.0)
        rm["wa"] = rm["w"] * rm["af"]
        pool = (rm.groupby(["site", "plot"], as_index=False)
                  .agg(wa=("wa", "sum"), w=("w", "sum"),
                       coverage=("coverage", "mean"), n_timepoints=("af", "size")))
        pool["af"] = pool["wa"] / pool["w"]
        pool["generation"] = g
        frames.append(pool.drop(columns=["wa", "w"]))
        print(f"  gen{g}: {len(rm)} samples -> {len(pool)} pools", flush=True)
    D = pd.concat(frames, ignore_index=True)
    D["site"] = D["site"].astype(int)
    return D


def sig_gardens(cfg) -> tuple[set, str, int]:
    """Significant gardens, the scan used, and how many gardens that scan ran in.

    The scan covers 30 gardens; the AF store has 31 (one garden has AF but no GWAS),
    so the denominator here must come from the scan, not from the plotted panels.
    """
    want = cfg.get("sigscan", "auto")
    best, bestn, bestcls, ntot = set(), -1, "", 0
    is_snp = int(cfg["ref_len"]) == 1 and int(cfg["alt_len"]) == 1
    # a SNP and an indel can share a position: never label an indel with SNP-scan hits
    auto = ("snp",) if is_snp else ("nonsnp", "sv")
    for cls in (auto if want == "auto" else (want,)):
        f = f"{GWASD}/persite_gwas_{cls}.npz"
        if not os.path.exists(f):
            continue
        z = np.load(f, allow_pickle=True)
        i = np.where((z["chrom"] == cfg["chrom"].lower()) & (z["pos"] == cfg["pos"]))[0]
        if not len(i):
            continue
        if cls != "snp" and os.path.exists(RESOLVED):
            # the npz carries no allele: at a multiallelic position its first row can be a
            # different allele than the one plotted. Take the gardens from the allele-resolved
            # hit table instead (absent there = significant nowhere).
            r = pd.read_csv(RESOLVED)
            r = r[(r.cls == cls) & (r.chrom == cfg["chrom"].lower()) & (r.pos == int(cfg["pos"]))
                  & (r.ref_len == int(cfg["ref_len"])) & (r.alt_len == int(cfg["alt_len"]))]
            if "store_row" in cfg and "key_rank" in r.columns:
                # same-length ALT alleles share the 4-field key; the store lists them in
                # panel order, so the plotted row's rank among them names the allele
                idx = np.load(f"{lib.AF_STORE}/index_nonsnp.npz", allow_pickle=True)
                same = np.where((idx["chrom"] == cfg["chrom"]) & (idx["pos"] == cfg["pos"])
                                & (idx["ref_len"] == cfg["ref_len"])
                                & (idx["alt_len"] == cfg["alt_len"]))[0]
                r = r[r.key_rank == int(np.searchsorted(same, int(cfg["store_row"])))]
            s = {int(g) for v in r.gardens for g in str(v).split(",")}
        else:
            nlp = -np.log10(np.clip(z["P"][int(i[0])], 1e-300, None))
            s = {int(g) for g, v in zip(z["sites"], nlp) if v >= BONF_GWAS[cls]}
        if len(s) > bestn:
            best, bestn, bestcls, ntot = s, len(s), cls, len(z["sites"])
    return best, bestcls, ntot


# Per-garden corner annotation and colourbar label for whichever climate axis orders the
# grid. Units: bio1-11 temperature-derived (degC; bio4 is sd*100, bio3 a ratio*100),
# bio12-19 precipitation (mm; bio15 is a CV in %), pc1-3 unitless z-scores.
AXIS_LABEL = {
    "bio1": "garden mean annual temperature, bio1 (°C)",
    "bio3": "garden isothermality, bio3 (bio2/bio7 x100)",
    "bio4": "garden temperature seasonality, bio4 (sd x100)",
    "bio6": "garden min temperature of coldest month, bio6 (°C)",
    "bio8": "garden mean temperature of wettest quarter, bio8 (°C)",
    "bio10": "garden mean temperature of warmest quarter, bio10 (°C)",
    "bio11": "garden mean temperature of coldest quarter, bio11 (°C)",
    "bio12": "garden annual precipitation, bio12 (mm)",
    "bio15": "garden precipitation seasonality, bio15 (CV %)",
    "bio17": "garden precipitation of driest quarter, bio17 (mm)",
    "pc1": "garden climate PC1 (z)", "pc2": "garden climate PC2 (z)",
    "pc3": "garden climate PC3 (z)",
}


def _fmt_axis(cax: str, v: float) -> str:
    if cax.startswith("pc"):
        return f"{cax} {v:+.2f}"
    n = int(cax[3:])
    if n in (3, 4, 15):
        return f"{v:.0f}"
    return f"{v:.0f} mm" if n >= 12 else f"{v:.1f} °C"


def main():
    cfg = json.loads(sys.argv[1])
    sym = cfg.get("sym", f"{cfg['chrom']}_{cfg['pos']}")
    ncol = int(cfg.get("ncol", 6))
    D = per_plot_af(cfg)
    p0 = founding_p0(cfg, cfg.get("store_row"))
    # Order/colour gardens by the axis this candidate's signal is actually on. Default
    # bio1, but 53 of the 89 screen candidates associate with bio12 (precipitation): for
    # those, sorting the grid by temperature scatters the gradient across the page and the
    # figure cannot be used to judge the climate criterion it was made to judge.
    cax = cfg.get("climate", "bio1")
    clim_tab = lib.load_climate()
    if cax in clim_tab.columns:
        clim = clim_tab[cax]
    else:
        # pc1-3 are not in the bioclim file: recompute them exactly as the GEA did
        # (axis_clusters.site_climate), which drops the site index -> restore it.
        sys.path.insert(0, os.path.join(HERE, "..", "convergence"))
        import axis_clusters as ac                                      # noqa: E402
        A = ac.site_climate()
        A.index = sorted(pd.read_csv(ac.POOLMETA).site.astype(int).unique())
        clim = A[cax]
    D["bio1"] = D.site.map(clim)          # `bio1` kept as the internal column name
    D = D.dropna(subset=["bio1"])
    sig, sigcls, n_scan = sig_gardens(cfg)

    site_b1 = D.groupby("site")["bio1"].first().sort_values()
    gardens = list(site_b1.index)
    nrow = int(np.ceil(len(gardens) / ncol))
    norm = mcolors.Normalize(site_b1.min(), site_b1.max())
    print(f"{sym} {cfg['chrom']}:{cfg['pos']}  p0={p0:.4f}  "
          f"{len(gardens)} gardens, {cax} {site_b1.min():.1f}..{site_b1.max():.1f}  "
          f"| significant in {len(sig)}/{n_scan} ({sigcls} scan)")

    # PLAIN LINEAR frequency, full range, NOTHING CLIPPED (user decision 2026-09-15).
    # Every pool is drawn at its real value on a shared linear axis. Trade-off, stated
    # rather than worked around: pool AF spans 0.0008-0.56, so the bulk of the data sits
    # low in the panel. Do NOT re-introduce clipping to compensate -- hiding pools that
    # still contribute to the plotted mean is what made garden 45 gen2 (pools 0.5605,
    # 0.4386, 0.0194, 0.0117, 0.0060) show a mean dot at 0.207 with no dot near it.
    site_mean_max = float(D.groupby(["site", "generation"]).af.mean().max())
    yhi = float(max(D.af.max(), p0)) * 1.06
    ylo = -yhi * 0.04
    print(f"  y-range {ylo:.3f}-{yhi:.3f} (LINEAR, full data, nothing clipped); "
          f"{len(D)} pools, AF {D.af.min():.4f}-{D.af.max():.4f}, "
          f"largest site mean {site_mean_max:.3f}")
    fig, axes = plt.subplots(nrow, ncol, figsize=(2.05 * ncol, 1.72 * nrow),
                             sharex=True, sharey=True)
    axes = np.atleast_1d(axes).ravel()

    for ax, site in zip(axes, gardens):
        s = D[D.site == site]
        b1 = float(s.bio1.iloc[0])
        col = CMAP(norm(b1))
        ax.axhline(p0, ls="--", lw=0.7, color="#999", zorder=1)
        # every sequenced pool, so the site mean is read against its own scatter
        jx = s.generation + np.random.default_rng(site).uniform(-.09, .09, len(s))
        ax.scatter(jx, s.af, s=7, color=col, alpha=0.38, edgecolors="none", zorder=2)
        m = s.groupby("generation")["af"].mean()
        m = pd.concat([pd.Series({0: p0}), m]).sort_index()
        ax.plot(m.index, m.values, "-o", color=col, lw=1.8, ms=4.2,
                mec="white", mew=0.7, zorder=4)
        # panel identity as a corner annotation (project convention: no subplot titles)
        ax.annotate(f"garden {site}", xy=(0.055, 0.93), xycoords="axes fraction",
                    fontsize=7.5, fontweight="bold", va="top", ha="left", color="#222")
        ax.annotate(_fmt_axis(cax, b1),
                    xy=(0.055, 0.78), xycoords="axes fraction",
                    fontsize=7, va="top", ha="left", color=col)
        if site in sig:
            ax.annotate("★", xy=(0.955, 0.93), xycoords="axes fraction", fontsize=10,
                        va="top", ha="right", color="#b2182b")
            for sp in ax.spines.values():
                sp.set_color("#b2182b"); sp.set_linewidth(1.5)
        else:
            for sp in ("top", "right"):
                ax.spines[sp].set_visible(False)
        ax.set_xticks([0, 1, 2, 3])
        ax.set_ylim(ylo, yhi)
        ax.tick_params(labelsize=7)
        ax.grid(True, lw=0.3, c="0.92")
        ax.set_axisbelow(True)

    for ax in axes[len(gardens):]:
        ax.axis("off")

    lab = cfg.get("label", f"{cfg['chrom']}:{cfg['pos']}")
    lab2 = "deletion" if cfg["ref_len"] > cfg["alt_len"] else "insertion"
    # gardens in the GWAS scan (30) vs in the AF store (31) differ
    # denominator = gardens the GWAS scan ran in (30), not panels drawn (31)
    fig.supxlabel("generation  (0 = founding SEEDMIX)", fontsize=11, y=0.005)
    fig.supylabel(f"alt allele frequency — {sym} {lab}", fontsize=11, x=0.055)
    sm = plt.cm.ScalarMappable(cmap=CMAP, norm=norm)
    cb = fig.colorbar(sm, ax=axes.tolist(), fraction=0.014, pad=0.012)
    cb.set_label(AXIS_LABEL.get(cax, f"garden {cax}"), fontsize=9)
    cb.ax.tick_params(labelsize=7.5)
    fig.text(0.055, 1.005,
             f"ONE variant: {cfg['chrom']}:{cfg['pos']:,}  ref {cfg['ref_len']} / alt "
             f"{cfg['alt_len']}  ({abs(cfg['alt_len']-cfg['ref_len'])} bp {sym} {lab2}) "
             f"— the marker significant in {len(sig)}/{n_scan} gardens of the {sigcls} scan",
             fontsize=9.5, color="#222", va="bottom", ha="left", fontweight="bold")
    fig.text(0.055, 0.982,
             "every dot is that SAME variant in one sequenced pool (one garden plot, one "
             "generation);  bold line = mean over the garden's pools;  "
             f"dashed = founding frequency p0 = {p0:.3f};  "
             f"★ / red frame = Bonferroni-significant in that garden;  every pool drawn at its real value (no clipping)",
             fontsize=8, color="#555", va="bottom", ha="left")

    os.makedirs(OUT, exist_ok=True)
    out = cfg.get("out", f"{OUT}/{sym}_garden_trajectories")
    fig.savefig(out + ".png", dpi=200, bbox_inches="tight")
    fig.savefig(out + ".pdf", bbox_inches="tight")
    print(f"wrote {out}.png/.pdf")


if __name__ == "__main__":
    main()
