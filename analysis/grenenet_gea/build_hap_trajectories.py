#!/usr/bin/env python
"""Pipeline B v2 — STEP 1 + STEP 2 on haplotype-frequency trajectories.

Builds, for every haplotype unit (70,384) and every site g, the per-generation
pooled site frequency p̄_{t,g} (Step 1) and its log-odds point variance V_{t,g}
(Step 2), following analysis/grenenet_gea/PIPELINE_B_POOLED_MODEL.md exactly.

Currency = haplotype frequency (each hapfreq column = one one-vs-rest allele).
Gen 0 = the SEEDMIX founding anchor (p0/v0 over the 8 reps), shared by all sites.

Locked choices honoured: persistent plots (present at gen 1,2,3); flower-weighted
pooling; ε=1e-3 logit clip; variance-components max(drift, sampling floor) on the
log-odds scale; drift shrinkage toward the per-haplotype pool, prior df d0=4.

Inputs (results/grenenet_gea/hapfreq/):
  hapfreq_matrix.npy (2168 × 70384), hapfreq_samples.txt, hapfreq_p0/v0_seedmix.npy
  + lib.pool_table() (sample → site/plot/generation/flowerscollected/coverage)
Outputs (results/grenenet_gea/hapfreq/pipelineB/):
  traj_p.npy  (n_sites × 4 × 70384)  pooled site freq p̄_{t,g}  (t=0..3)
  traj_v.npy  (n_sites × 4 × 70384)  log-odds point variance V_{t,g}
  traj_sites.npy (n_sites,)  site ids   traj_nplots.npy (n_sites,)  persistent n_g

Env: kmate.
"""
import os, sys
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lib

H = os.environ.get("HF_DIR", "results/grenenet_gea/hapfreq")   # hapfreq store (clq90: hapfreq_clq90)
# MODE: "persistent" = plots present at ALL of gen 1,2,3 (19 sites, 4-pt trajectories);
#       "varlen" = each site uses the LONGEST available consistent trajectory (gen-0 anchor
#       + whatever evolved gens it has a fixed plot set for) → more sites, shorter series.
MODE = os.environ.get("TRAJ_MODE", "persistent")
OUT = f"{H}/pipelineB" if MODE == "persistent" else f"{H}/pipelineB_{MODE}"
EPS = 1e-3
D0 = 4.0                         # drift-variance shrinkage prior df (locked)
# longest/earliest-first preference for the per-site evolved-generation set (varlen)
GEN_PRIORITY = [{1, 2, 3}, {1, 2}, {2, 3}, {1, 3}, {1}, {2}, {3}]


def logit(p):
    p = np.clip(p, EPS, 1 - EPS)
    return np.log(p / (1 - p))


def main():
    os.makedirs(OUT, exist_ok=True)
    mat = np.load(f"{H}/hapfreq_matrix.npy")                     # (2168, nhap)
    samples = [l.strip() for l in open(f"{H}/hapfreq_samples.txt") if l.strip()]
    sidx = {s: i for i, s in enumerate(samples)}
    nhap = mat.shape[1]
    p0 = np.load(f"{H}/hapfreq_p0_seedmix.npy").astype(np.float64)
    v0 = np.load(f"{H}/hapfreq_v0_seedmix.npy").astype(np.float64)

    pt = lib.pool_table()
    pt = pt[pt.sampleid.astype(str).isin(sidx)].copy()
    pt["row"] = pt.sampleid.astype(str).map(sidx)

    # --- collapse timepoints within each (site,gen,plot): flower-weighted hapfreq,
    #     total flowers (sum), mean coverage ---------------------------------------
    plotcell = {}           # (site,gen,plot) -> dict(p=vec, flowers=, cov=)
    for (site, gen, plot), grp in pt.groupby(["site", "generation", "plot"]):
        w = grp.flowerscollected.to_numpy(float)
        w = np.where(np.isfinite(w) & (w > 0), w, 1.0)          # guard 0/NaN flowers
        rows = grp.row.to_numpy()
        p = (mat[rows] * w[:, None]).sum(0) / w.sum()           # flower-wt timepoint merge
        cov = np.nanmean(grp.coverage.to_numpy(float))
        plotcell[(int(site), int(gen), int(plot))] = dict(
            p=p.astype(np.float64), flowers=float(w.sum()),
            cov=float(cov if np.isfinite(cov) and cov > 0 else 1.0))

    # --- per-site evolved-generation set + fixed plot set ---------------------------
    # persistent: require a plot present at ALL of gen 1,2,3.
    # varlen: take the longest (then earliest) generation subset for which the site has a
    #         fixed (persistent) plot set; gen-0 anchor always included.
    present = {}                                                # site -> {plot: set(gens)}
    for (site, gen, plot) in plotcell:
        present.setdefault(site, {}).setdefault(plot, set()).add(gen)
    site_gens, site_plots = {}, {}
    for site, d in present.items():
        if MODE == "persistent":
            pls = [pl for pl, gs in d.items() if {1, 2, 3} <= gs]
            if pls:
                site_gens[site], site_plots[site] = [1, 2, 3], pls
        else:
            for Sset in GEN_PRIORITY:                           # longest/earliest first
                pls = [pl for pl, gs in d.items() if Sset <= gs]
                if pls:
                    site_gens[site], site_plots[site] = sorted(Sset), pls
                    break
    sites = sorted(site_gens)
    nplots = np.array([len(site_plots[s]) for s in sites])
    ngens = np.array([len(site_gens[s]) for s in sites])        # evolved gens used (1..3)

    P = np.zeros((len(sites), 4, nhap), np.float64)
    P[:, 0, :] = p0                                             # gen-0 SEEDMIX anchor
    V = np.full((len(sites), 4, nhap), np.inf, np.float64)      # unused gens → ω=0 downstream
    V[:, 0, :] = v0
    drift_raw = np.full((len(sites), 4, nhap), np.nan)          # σ²_{t,g} (pre-shrink)
    vsamp_bar = np.zeros((len(sites), 4, nhap))                 # mean plot sampling var

    # --- STEP 1: pooled site freq p̄_{t,g} over the site's used gens ----------------
    for gi, site in enumerate(sites):
        pls = site_plots[site]
        for gen in site_gens[site]:
            ti = gen                                            # time index == generation
            cells = [plotcell[(site, gen, pl)] for pl in pls]
            fw = np.array([c["flowers"] for c in cells])
            Pp = np.vstack([c["p"] for c in cells])             # (n_g, nhap) per-plot freq
            pbar = (Pp * fw[:, None]).sum(0) / fw.sum()         # Step 1 pooled site freq
            P[gi, ti, :] = pbar
            if len(cells) >= 2:                                 # drift among plots
                drift_raw[gi, ti, :] = logit(Pp).var(0, ddof=1)
            Neff_inv = np.array([1.0 / (2 * c["flowers"]) + 1.0 / c["cov"] for c in cells])
            pc = np.clip(pbar, EPS, 1 - EPS)
            vsamp_bar[gi, ti, :] = (Neff_inv[:, None] / (pc * (1 - pc))[None, :]).mean(0)

    # --- STEP 2: per-haplotype drift pool, shrink, combine, /n_g (used gens only) ----
    sigma_pool = np.nanmean(drift_raw.reshape(-1, nhap), axis=0)
    sigma_pool = np.where(np.isfinite(sigma_pool), sigma_pool, 0.0)
    for gi, site in enumerate(sites):
        n_g = nplots[gi]
        for gen in site_gens[site]:
            ti = gen
            raw = np.where(np.isfinite(drift_raw[gi, ti, :]), drift_raw[gi, ti, :], sigma_pool)
            shrunk = ((n_g - 1) * raw + D0 * sigma_pool) / ((n_g - 1) + D0)
            V[gi, ti, :] = np.maximum(shrunk, vsamp_bar[gi, ti, :]) / n_g

    os.makedirs(OUT, exist_ok=True)
    np.save(f"{OUT}/traj_p.npy", P.astype(np.float32))
    np.save(f"{OUT}/traj_v.npy", V.astype(np.float32))
    np.save(f"{OUT}/traj_sites.npy", np.array(sites))
    np.save(f"{OUT}/traj_nplots.npy", nplots)
    np.save(f"{OUT}/traj_ngens.npy", ngens)
    from collections import Counter
    lendist = Counter(ngens + 1)                                # trajectory length incl. gen 0
    print(f"[done MODE={MODE}] {len(sites)} sites (n_g {nplots.min()}–{nplots.max()}, "
          f"median {int(np.median(nplots))}); P,V = {P.shape} -> {OUT}/")
    print(f"  trajectory length (incl gen0): " +
          ", ".join(f"{L}pt×{n}" for L, n in sorted(lendist.items())))
    print(f"  drift σ²_pool median {np.median(sigma_pool):.4f}")


if __name__ == "__main__":
    main()
