#!/usr/bin/env python
"""Step A — per-plot table joining kMate founder composition, realized census fitness
(cumulative flowering-individual count from survival.csv, cap-aware / rank-normalized),
and survival; plus a first-pass test of whether inferred selection s_f predicts realized
demographic fitness.

Fitness = cumulative census Sum_t flowerstotal over gens with survival in {0,1} (death=real 0,
drop nan/-1), then rank/quantile-normalized (immune to the ~100 collection-cap censoring;
uses only ordering). See memory grenenet-survival-flowers-table.

Inputs (all cached, light):
  - results/grenenet_gea/fitness/sample_genome_h.npz   (build_sample_h_cache.py)
  - results/grenenet_gea/hapfreq/cross_site_S_matrix_gen1.npz  (per-site s_f, 30 sites)
  - /global/home/users/tbellg/scratch/grene/data/survival.csv
  - lib.pool_table()  (sampleid -> site,plot,generation,flowerscollected)

Outputs: results/grenenet_gea/fitness/{plot_fitness_composition.csv, plot_h.npz, stepA_summary.json}
Env: kmate.
"""
import os, sys, json
import numpy as np
import pandas as pd
from scipy import stats
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lib

FIT = "results/grenenet_gea/fitness"
H = "results/grenenet_gea/hapfreq"
SURV = "/global/home/users/tbellg/scratch/grene/data/survival.csv"
GENS_SEQ = (1, 2, 3)
GENS_CENSUS = (1, 2, 3, 4, 5)


def qn(x):
    """rank-based inverse-normal (quantile) transform of a 1-D array (NaNs preserved)."""
    x = np.asarray(x, float); out = np.full_like(x, np.nan); m = np.isfinite(x)
    r = stats.rankdata(x[m]); out[m] = stats.norm.ppf((r - 0.5) / m.sum())
    return out


def load_fitness():
    """Per-plot lifetime census fitness + survival summary from survival.csv."""
    sv = pd.read_csv(SURV)
    sv = sv[pd.to_numeric(sv["site"], errors="coerce").notna()
            & pd.to_numeric(sv["plot"], errors="coerce").notna()].copy()
    sv["site"] = sv["site"].astype(float).astype(int); sv["plot"] = sv["plot"].astype(float).astype(int)
    cum = np.zeros(len(sv)); ngen = np.zeros(len(sv)); capped = np.zeros(len(sv), bool)
    ever_died = np.zeros(len(sv), bool); last_state = np.full(len(sv), np.nan)
    for g in GENS_CENSUS:
        s = pd.to_numeric(sv[f"{g}_survival"], errors="coerce").values
        f = pd.to_numeric(sv[f"{g}_flowerstotal"], errors="coerce").values
        use = np.isin(s, [0, 1]) & np.isfinite(f)
        cum[use] += f[use]; ngen[use] += 1
        capped |= use & (f >= 100)
        ever_died |= (s == 0)
        seen = np.isin(s, [0, 1]); last_state[seen] = s[seen]      # last observed 0/1
    out = pd.DataFrame({"site": sv["site"].values, "plot": sv["plot"].values, "cum_flowers": cum,
                        "n_obs_gen": ngen.astype(int), "ever_capped": capped,
                        "ever_died": ever_died, "terminal_survival": last_state})
    return out[out.n_obs_gen > 0].reset_index(drop=True)


def main():
    cache = f"{FIT}/sample_genome_h.npz"
    if not os.path.exists(cache):
        sys.exit(f"[wait] per-sample h cache not ready: {cache} (run build_sample_h_cache.py)")
    d = np.load(cache, allow_pickle=True)
    samp_ix = {s: i for i, s in enumerate(d["samples"].astype(str))}
    Hmat = d["H"]; founders = d["founders"].astype(str)
    Sd = np.load(f"{H}/cross_site_S_matrix_gen1.npz", allow_pickle=True)
    assert np.array_equal(Sd["founders"].astype(str), founders), "founder order mismatch"
    Sf = {int(s): Sd["S"][i] for i, s in enumerate(Sd["sites"])}       # site -> s_f (231,)

    pt = lib.pool_table()
    pt = pt[pt.sampleid.astype(str).isin(samp_ix)].copy()
    fit = load_fitness()

    # ---- per (site,plot): flower-weighted mean composition per sequenced gen ----
    rows = []; plot_h = {}
    for (site, plot), g in pt.groupby(["site", "plot"]):
        if int(site) not in Sf:
            continue
        frow = fit[(fit["site"] == int(site)) & (fit["plot"] == int(plot))]
        if frow.empty:
            continue
        comp = {}
        for gen, gg in g.groupby("generation"):
            idx = [samp_ix[str(s)] for s in gg.sampleid.astype(str)]
            w = gg.flowerscollected.to_numpy(float); w = np.where(np.isfinite(w) & (w > 0), w, 1.0)
            comp[int(gen)] = (Hmat[idx] * w[:, None]).sum(0) / w.sum()
        seq_gens = sorted(comp)
        h_early = comp[seq_gens[0]]; h_final = comp[seq_gens[-1]]
        sf = Sf[int(site)]; top = sf >= np.quantile(sf, 0.9)          # site's top-decile winners
        z = (sf - sf.mean()) / (sf.std() + 1e-12)
        r = frow.iloc[0]
        rows.append(dict(site=int(site), plot=int(plot),
                         seq_gens="".join(map(str, seq_gens)),
                         cum_flowers=r.cum_flowers, n_obs_gen=int(r.n_obs_gen),
                         ever_capped=bool(r.ever_capped), ever_died=bool(r.ever_died),
                         terminal_survival=r.terminal_survival,
                         winner_share_early=float(h_early[top].sum()),
                         winner_share_final=float(h_final[top].sum()),
                         bv_early=float(h_early @ z), bv_final=float(h_final @ z),
                         eff_n_final=float(1.0 / np.sum(h_final ** 2))))
        plot_h[f"{site}_{plot}_final"] = h_final; plot_h[f"{site}_{plot}_early"] = h_early
    df = pd.DataFrame(rows)

    # ---- fitness transforms: global QN + within-site rank ----
    df["fit_qn"] = qn(df.cum_flowers.values)
    df["fit_rank_site"] = df.groupby("site").cum_flowers.rank(pct=True)
    os.makedirs(FIT, exist_ok=True)
    df.to_csv(f"{FIT}/plot_fitness_composition.csv", index=False)
    np.savez(f"{FIT}/plot_h.npz", **plot_h, founders=founders)

    # ---- diagnostics ----------------------------------------------------------
    print(f"\n{len(df)} sequenced plots with fitness across {df.site.nunique()} sites")
    print("plots/site:", df.groupby("site").size().describe()[["min", "50%", "max"]].to_dict())

    # (a) among-plot composition variance within sites (is there anything to test?)
    cvs = []
    for site, g in df.groupby("site"):
        hs = np.vstack([plot_h[f"{site}_{p}_final"] for p in g["plot"]])
        if len(hs) < 3:
            continue
        big = hs.mean(0) > 0.02
        if big.sum():
            cvs.append(np.nanmean(hs[:, big].std(0) / (hs[:, big].mean(0) + 1e-9)))
    print(f"(a) within-site among-plot composition CV (major founders): "
          f"median {np.nanmedian(cvs):.3f} over {len(cvs)} sites "
          f"(low CV => plots converged => low within-site power)")

    # (b) does winner-enrichment / selection-alignment predict census fitness? (within-site)
    def pooled_within_site(xcol, ycol, kind="spearman"):
        zs, ns = [], []
        for site, g in df.groupby("site"):
            x = g[xcol].values; y = g[ycol].values
            m = np.isfinite(x) & np.isfinite(y)
            if m.sum() < 5 or np.std(x[m]) == 0:
                continue
            r = (stats.spearmanr(x[m], y[m])[0] if kind == "spearman"
                 else stats.pointbiserialr(y[m].astype(int), x[m])[0])
            if np.isfinite(r):
                zs.append(np.arctanh(np.clip(r, -.999, .999))); ns.append(m.sum() - 3)
        zs = np.array(zs); ns = np.array(ns)
        zbar = np.sum(zs * ns) / ns.sum(); se = 1 / np.sqrt(ns.sum())
        return dict(n_sites=len(zs), mean_r=float(np.tanh(zbar)),
                    z=float(zbar / se), p=float(2 * stats.norm.sf(abs(zbar / se))))

    res = {}
    for xcol in ["winner_share_final", "bv_final", "winner_share_early", "bv_early"]:
        res[f"{xcol}~fit"] = pooled_within_site(xcol, "cum_flowers")
        print(f"(b) within-site Spearman({xcol}, census fitness): "
              f"mean r={res[f'{xcol}~fit']['mean_r']:+.3f}  p={res[f'{xcol}~fit']['p']:.3g}  "
              f"({res[f'{xcol}~fit']['n_sites']} sites)")

    # (c) survival ~ selection alignment (does composition predict the plot persisting?)
    surv = df[df.terminal_survival.isin([0, 1])].copy()
    for xcol in ["winner_share_final", "bv_final"]:
        zs, ns = [], []
        for site, g in surv.groupby("site"):
            if g.terminal_survival.nunique() < 2 or len(g) < 6:
                continue
            r = stats.pointbiserialr(g.terminal_survival.astype(int), g[xcol])[0]
            if np.isfinite(r):
                zs.append(np.arctanh(np.clip(r, -.999, .999))); ns.append(len(g) - 3)
        if ns:
            zbar = np.sum(np.array(zs) * ns) / np.sum(ns); se = 1 / np.sqrt(np.sum(ns))
            res[f"{xcol}~surv"] = dict(n_sites=len(zs), mean_r=float(np.tanh(zbar)),
                                       z=float(zbar / se), p=float(2 * stats.norm.sf(abs(zbar / se))))
            print(f"(c) within-site pt-biserial({xcol}, terminal survival): "
                  f"mean r={np.tanh(zbar):+.3f}  p={2*stats.norm.sf(abs(zbar/se)):.3g}  ({len(zs)} sites)")

    json.dump(res, open(f"{FIT}/stepA_summary.json", "w"), indent=2)
    print(f"\n[done] {FIT}/plot_fitness_composition.csv ({len(df)} plots) + plot_h.npz + stepA_summary.json")


if __name__ == "__main__":
    main()
