#!/usr/bin/env python
"""Ecotype-level consistency-vs-drift selection at ONE site (default site 4).

The recombination check showed block-local h == genome-wide h (r~0.99): haploblocks
just project ~25 ecotype trajectories. So run the selection test where it belongs --
on the founder/ecotype frequencies (global h), low-dim and properly replicated.

Per founder f, per replicate plot j at the site:
  trajectory  h_{f,j,t}  over t = 0 (seedmix p0), 1, 2, 3  (genome-wide global h, mean
  over Chr1..5; per-(gen,plot) flower-weighted across timepoints)
Site statistics: dfreq_f = mean_j (h_{f,j,3} - p0_f)  (net change, interpretable),
  s_bar_f = mean_j logit-slope (selection coefficient), concordance = #plots same sign.

DRIFT NULL (the point): replicates start identical at p0 and, if NEUTRAL, drift
INDEPENDENTLY -> dfreq ~ 0 with spread set by drift+sampling. We:
  * estimate one pooled effective size Ne by matching the observed among-plot variance
    of gen-3 frequencies to Wright-Fisher expectation p(1-p)[1-(1-1/Ne)^3]
    (using only well-resolved founders, freq in 0.02..0.4);
  * simulate neutral binomial-WF trajectories from p0 (Ne, 3 gens, n_plots reps), NSIM
    times, with shared founding error injected from the 8 seedmix reps (so concordance
    from a mis-estimated founding freq is NOT counted as selection) -> null dfreq;
  * two-sided p = tail of |dfreq| under the null; BH-FDR across founders.

The founding seed mix is near-uniform/diverse (eff_n~160, max ~1.7%); the WINNERS start
rare (~0.004) and rise -- so keep founders present at founding OR reaching gen-3 signal.

Output: ranked winners/losers with effect, concordance, drift p; figure + csv. Env: kmate.
"""
import os, sys, glob
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lib

# Repointed 2026-07-07 to the --unit chrom cohort (was the stale window store
# results/grenenet_kmate_window[_seedmix]). genome_h below reads the chrom format.
WIN = lib.OUT
SEED = lib.SEEDMIX
OUT = "analysis/grenenet_gea/archive/window_hapfreq_retired/hapfreq"
SITE = int(os.environ.get("SITE", 4))
CHROMS = ["Chr1", "Chr2", "Chr3", "Chr4", "Chr5"]
P0MIN = 0.002         # present at founding
G3MIN = 0.005         # or reaches >0.5% by gen 3 (catches the rare-start winners)
EPS = 1e-3
NSIM = 4000
Tg = np.array([0.0, 1.0, 2.0, 3.0])
rng = np.random.RandomState(0)


def logit(p):
    p = np.clip(p, EPS, 1 - EPS)
    return np.log(p / (1 - p))


def genome_h(samp, base):
    """genome-wide founder h for one sample = mean over chroms of the per-chrom h.

    Reads the --unit chrom cohort format `{base}/{samp}_{ch}.h_per_chrom.npz`, key
    `{ch}` (a single 231-vector per chrom). (Was the window store's
    `_ch.h_blocks_per_chrom.npz` key `{ch}_global_h` — same quantity, different file.)
    """
    gs = []
    for ch in CHROMS:
        f = f"{base}/{samp}_{ch}.h_per_chrom.npz"
        if not os.path.exists(f):
            return None
        gs.append(np.load(f, allow_pickle=True)[ch].astype(np.float64))
    return np.mean(gs, 0)


def slope(y):
    """OLS slope of y (… x 4 over t=0..3); Σ(t-1.5)²=5."""
    return ((Tg - 1.5) * y).sum(-1) / 5.0


def main():
    founders = np.load(glob.glob(f"{SEED}/*_Chr1.h_per_chrom.npz")[0],
                       allow_pickle=True)["founders"].astype(str)
    nF = len(founders)

    # ---- founding p0 + uncertainty over the 8 seedmix reps ----
    seeds = sorted({p.split("/")[-1].split("_Chr")[0]
                    for p in glob.glob(f"{SEED}/*_Chr1.h_per_chrom.npz")})
    P0 = np.vstack([genome_h(s, SEED) for s in seeds])
    p0 = P0.mean(0); v0 = P0.var(0, ddof=1)
    print(f"founding seed mix: eff_n {1/(p0**2).sum():.0f}, max {p0.max():.3f} "
          f"({len(seeds)} reps)")

    # ---- site replicate trajectories: per plot, gen 1/2/3 genome-wide h ----
    pt = lib.pool_table()
    s = pt[pt.site == SITE].copy()
    s = s[s.sampleid.astype(str).apply(
        lambda x: os.path.exists(f"{WIN}/{x}_Chr1.h_per_chrom.npz"))]
    cell = {}
    for (gen, plot), g in s.groupby(["generation", "plot"]):
        hs, ws = [], []
        for _, r in g.iterrows():
            h = genome_h(str(r.sampleid), WIN)
            if h is not None:
                w = r.flowerscollected if np.isfinite(r.flowerscollected) and r.flowerscollected > 0 else 1.0
                hs.append(h); ws.append(w)
        if hs:
            ws = np.array(ws)
            cell[(int(gen), int(plot))] = (np.vstack(hs) * ws[:, None]).sum(0) / ws.sum()
    present = {}
    for (gen, plot) in cell:
        present.setdefault(plot, set()).add(gen)
    plots = sorted(pl for pl, gs in present.items() if {1, 2, 3} <= gs)
    n = len(plots)

    H = np.zeros((n, 4, nF)); H[:, 0, :] = p0[None, :]
    for j, pl in enumerate(plots):
        for gen in (1, 2, 3):
            H[j, gen, :] = cell[(gen, pl)]
    print(f"site {SITE}: {n} persistent plots, {nF} founders")

    # ---- observed site statistics ----
    freq3 = H[:, 3, :].mean(0)
    dplot = H[:, 3, :] - p0[None, :]                            # per-plot net change
    dfreq = dplot.mean(0)
    nup = (dplot > 0).sum(0)
    conc = np.maximum(nup, n - nup) / n
    s_bar = slope(logit(H).transpose(0, 2, 1)).mean(0)
    keep = (p0 > P0MIN) | (freq3 > G3MIN)

    # ---- pooled Ne from among-plot variance of WELL-RESOLVED founders ----
    res = keep & (freq3 > 0.02) & (freq3 < 0.4)
    pbar = freq3[res]
    A = float(np.clip(np.median(H[:, 3, res].var(0, ddof=1) / (pbar * (1 - pbar))), 1e-4, 0.99))
    Ne = float(np.clip(1.0 / (1.0 - (1.0 - A) ** (1.0 / 3.0)), 5, 1e5))
    print(f"  Ne ~ {Ne:.0f} (3-gen drift accrual {A:.3f}; from {int(res.sum())} resolved founders)")

    # ---- WF drift null: neutral sims from p0 (+shared founding error) ----
    def sim_dfreq():
        p0e = np.clip(p0 + rng.normal(0, np.sqrt(v0)), 0, 1)
        cur = np.repeat(p0e[None, :], n, 0)
        for _ in (1, 2, 3):
            cur = rng.binomial(int(Ne), np.clip(cur, 0, 1)) / Ne
        return cur.mean(0) - p0
    null = np.vstack([sim_dfreq() for _ in range(NSIM)])
    p_drift = ((np.abs(null) >= np.abs(dfreq)[None, :]).sum(0) + 1) / (NSIM + 1)
    null_sd = null.std(0)

    d = pd.DataFrame({"founder": founders, "p0": p0, "freq_gen3": freq3, "d_freq": dfreq,
                      "s_bar": s_bar, "n_up": nup, "n_plots": n, "concordance": conc,
                      "null_sd": null_sd, "p_drift": p_drift})[keep].copy()
    d = d.sort_values("d_freq", ascending=False).reset_index(drop=True)
    m = len(d); o = d.p_drift.values.argsort()
    q = np.empty(m); q[o] = np.minimum.accumulate((d.p_drift.values[o] * m / (np.arange(m) + 1))[::-1])[::-1]
    d["q"] = np.clip(q, 0, 1)
    d.to_csv(f"{OUT}/site{SITE}_ecotype_selection.csv", index=False)

    sig = d[d.q < 0.05]
    cols = ["founder", "p0", "freq_gen3", "d_freq", "s_bar", "n_up", "concordance", "p_drift", "q"]
    print(f"\n  {m} founders analysed; beyond drift (q<0.05): {len(sig)} "
          f"({int((sig.d_freq>0).sum())} winners, {int((sig.d_freq<0).sum())} losers)")
    print("\n  TOP WINNERS (rose most):"); print(d.head(10)[cols].to_string(index=False))
    print("\n  TOP LOSERS (fell most):");  print(d.tail(10)[cols].to_string(index=False))

    # ---- figure ----
    fig, (a0, a1) = plt.subplots(1, 2, figsize=(13, 6.5), gridspec_kw={"width_ratios": [1.6, 1]})
    dd = d.sort_values("d_freq").reset_index(drop=True)
    y = np.arange(len(dd))
    col = np.where(dd.q < 0.05, np.where(dd.d_freq > 0, "#c0392b", "#2471a3"), "#bdbdbd")
    a0.errorbar(dd.d_freq, y, xerr=1.96 * dd.null_sd, fmt="none", ecolor="#e5e5e5", zorder=1)
    a0.scatter(dd.d_freq, y, c=col, s=22, zorder=2)
    a0.axvline(0, color="k", lw=.6)
    for _, r in pd.concat([dd.head(8), dd.tail(8)]).iterrows():
        yi = int(dd.index[dd.founder == r.founder][0])
        a0.annotate(r.founder, (r.d_freq, yi), fontsize=6, va="center",
                    ha="right" if r.d_freq < 0 else "left")
    a0.set_xlabel(r"site mean frequency change  $\Delta f$  (gen 0$\to$3)")
    a0.set_ylabel(f"ecotype rank ({m} founders)")
    a0.set_title(f"Site {SITE}: which ecotypes won/lost beyond drift\n"
                 f"(red/blue = q<0.05; grey band = drift-null 95% CI; Ne~{Ne:.0f})",
                 loc="left", fontsize=10)
    a0.spines[["top", "right"]].set_visible(False)

    a1.scatter(dd.concordance, dd.d_freq, c=col, s=22)
    a1.axhline(0, color="k", lw=.6)
    a1.set_xlabel(f"replicate concordance (frac of {n} plots same direction)")
    a1.set_ylabel(r"$\Delta f$")
    a1.set_title("effect vs replicate concordance", loc="left", fontsize=10)
    a1.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    out = f"{OUT}/site{SITE}_ecotype_selection.png"
    fig.savefig(out, dpi=150)
    print(f"\n[done] {out}\n[done] {OUT}/site{SITE}_ecotype_selection.csv")


if __name__ == "__main__":
    main()
