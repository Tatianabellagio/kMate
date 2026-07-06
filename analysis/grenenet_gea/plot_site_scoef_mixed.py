#!/usr/bin/env python
"""Site-4 Manhattan of the RANDOM-SLOPE selection coefficient — plots as independent.

Model (per haplotype, per site):  logit(p_plot,t) ~ gen + (gen | plot)
i.e. each plot is its own trajectory; the site slope beta1 is the population slope and
its uncertainty comes from how much the per-plot slopes DISAGREE (the random-slope
variance = the among-replicate drift yardstick).

For this balanced design (every persistent plot has the same generations 0..3) the
mixed-model fixed-effect slope = MEAN of the per-plot OLS slopes, and its SE = the
among-plot SE of those slopes. So we compute, per haplotype:
    slope_p = Σ_t (t-1.5)·logit(p_{p,t}) / 5        (per-plot OLS slope, t=0..3)
    beta1   = mean_p slope_p
    SE      = sd_p(slope_p)/sqrt(n_plots)            (random-slope / among-plot SE)
    z = beta1/SE,  p = 2·t.sf(|z|, df=n_plots-1)
Gen 0 = shared seedmix anchor (same p0 in every plot), mirroring the pooled plot.

Contrast with site4_scoef_manhattan.png (pooled mean + model-based Wald): there the SE
ignored replicate disagreement and the panel hit -log10p~40; here the SE is set by the
plots, so drift-driven movers lose significance.  Env: kmate. SITE via env.
"""
import os, sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lib

H = "results/grenenet_gea/hapfreq"
B = f"{H}/pipelineB"
SITE = int(os.environ.get("SITE", 4))
EPS = 1e-3
Tg = np.array([0.0, 1.0, 2.0, 3.0])
CHROMS = ["Chr1", "Chr2", "Chr3", "Chr4", "Chr5"]


def logit(p):
    p = np.clip(p, EPS, 1 - EPS)
    return np.log(p / (1 - p))


def main():
    mat = np.load(f"{H}/hapfreq_matrix.npy")
    samples = [l.strip() for l in open(f"{H}/hapfreq_samples.txt") if l.strip()]
    sidx = {s: i for i, s in enumerate(samples)}
    p0 = np.load(f"{H}/hapfreq_p0_seedmix.npy").astype(np.float64)
    reg = pd.read_csv(f"{H}/hapfreq_registry.csv")

    pt = lib.pool_table()
    pt = pt[(pt.site == SITE) & pt.sampleid.astype(str).isin(sidx)].copy()
    pt["row"] = pt.sampleid.astype(str).map(sidx)

    # per (gen,plot): flower-weighted hapfreq across timepoints
    cell = {}
    for (gen, plot), g in pt.groupby(["generation", "plot"]):
        w = g.flowerscollected.to_numpy(float)
        w = np.where(np.isfinite(w) & (w > 0), w, 1.0)
        cell[(int(gen), int(plot))] = (mat[g.row.to_numpy()] * w[:, None]).sum(0) / w.sum()

    present = {}
    for (gen, plot) in cell:
        present.setdefault(plot, set()).add(gen)
    plots = sorted(pl for pl, gs in present.items() if {1, 2, 3} <= gs)
    n = len(plots)

    # per-plot 4-point trajectory on log-odds; gen0 = shared p0 anchor
    nhap = mat.shape[1]
    y = np.empty((n, 4, nhap))
    for pi, pl in enumerate(plots):
        y[pi, 0, :] = logit(p0)
        for gen in (1, 2, 3):
            y[pi, gen, :] = logit(cell[(gen, pl)])

    # per-plot OLS slope (t=0..3, Σ(t-1.5)²=5), then mean ± among-plot SE
    coef = (Tg - 1.5)[None, :, None]                     # (1,4,1)
    slope_p = (coef * y).sum(1) / 5.0                    # (n_plots, nhap)
    beta1 = slope_p.mean(0)
    se = slope_p.std(0, ddof=1) / np.sqrt(n)
    z = beta1 / np.clip(se, 1e-12, None)
    p = 2 * stats.t.sf(np.abs(z), df=n - 1)

    d = pd.DataFrame({"chrom": reg.chrom, "start": reg.unit_start, "end": reg.unit_end,
                      "panel_freq": reg.panel_freq, "beta1": beta1, "se": se,
                      "z": z, "p": p})
    d = d[reg.covered & reg.panel_freq.between(0.05, 0.95)].reset_index(drop=True)
    d["mid"] = (d.start + d.end) / 2
    d["nlp"] = -np.log10(np.clip(d.p, 1e-300, 1))

    off, centers, x = 0.0, [], np.zeros(len(d))
    for ch in CHROMS:
        m = (d.chrom == ch).to_numpy()
        if not m.any():
            continue
        cmax = d.loc[m, "mid"].max()
        x[m] = d.loc[m, "mid"] + off
        centers.append(off + cmax / 2); off += cmax * 1.02
    d["x"] = x
    up = d.beta1 > 0
    bonf = 0.05 / len(d)

    fig, (a0, a1) = plt.subplots(2, 1, figsize=(11, 7), sharex=True)
    a0.scatter(d.x[up], d.beta1[up], s=5, c="#c0392b", alpha=.5, edgecolors="none",
               rasterized=True, label="rose")
    a0.scatter(d.x[~up], d.beta1[~up], s=5, c="#2471a3", alpha=.5, edgecolors="none",
               rasterized=True, label="fell")
    a0.axhline(0, color="k", lw=.6)
    a0.set_ylabel(r"random-slope $\beta_1$  (log-odds / gen)")
    a0.set_title(f"Site {SITE} — random-slope $\\beta_1$, {n} plots as independent replicates "
                 f"({len(d):,} testable haplotypes)", loc="left", fontsize=11)
    a0.legend(fontsize=8, frameon=False, loc="upper right")

    a1.scatter(d.x[up], d.nlp[up], s=5, c="#c0392b", alpha=.5, edgecolors="none", rasterized=True)
    a1.scatter(d.x[~up], d.nlp[~up], s=5, c="#2471a3", alpha=.5, edgecolors="none", rasterized=True)
    a1.axhline(-np.log10(bonf), color="firebrick", lw=.9, ls="--",
               label=f"Bonferroni 0.05/{len(d)}")
    a1.set_ylabel(r"$-\log_{10}\,p$  (among-plot $t$, df=%d)" % (n - 1))
    a1.set_title("Significance from REPLICATE disagreement (drift-controlled), not model SE",
                 loc="left", fontsize=10)
    a1.set_xticks(centers); a1.set_xticklabels(CHROMS); a1.set_xlabel("genome position")
    a1.legend(fontsize=8, frameon=False, loc="upper right")
    for ax in (a0, a1):
        ax.spines[["top", "right"]].set_visible(False)

    fig.tight_layout()
    out = f"{B}/site{SITE}_scoef_mixed_manhattan.png"
    fig.savefig(out, dpi=150)
    print(f"[done] {out}")
    print(f"  site {SITE}: {n} plots, {len(d):,} testable haplotypes")
    print(f"  max -log10p = {d.nlp.max():.1f}  (pooled-mean Wald version maxed ~40)")
    print(f"  Bonferroni hits (p<{bonf:.1e}): {int((d.p<bonf).sum()):,}")
    print(f"  raw p<0.05: {int((d.p<0.05).sum()):,}  ({100*(d.p<0.05).mean():.0f}%)")
    top = d.reindex(d.nlp.sort_values(ascending=False).index).head(8)
    print(top[["chrom", "start", "end", "panel_freq", "beta1", "se", "z", "p"]].to_string(index=False))


if __name__ == "__main__":
    main()
