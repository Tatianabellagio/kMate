#!/usr/bin/env python
"""Pipeline B v2 — option B: permutation-calibrated block-WZA (Weighted-Z Analysis).

Replaces the per-block min-p Šidák aggregation (build_hap_gea.py) with Booker &
Whitlock's weighted-Z, calibrated by the SAME site-permutation null. For each climate
assignment (real + N_perm permutations of bio1 across the 19 sites) we:
  1. recompute the per-haplotype climate slope β1 (Step 4, fixed random-effects weights),
  2. robust-standardize β1 across the kept haplotypes → per-hap z (≈ GEA Z-score),
  3. aggregate to each block:  WZA_b = Σ_i w_i z_i / sqrt(Σ_i w_i²),  weights
     w_i = panel het = panel_freq(1-panel_freq), over the block's kept (k-1
     independent one-vs-rest) haplotypes.
The permutation distribution of WZA_b self-calibrates for block size + within-block LD
(a big block has a big WZA under BOTH real and permuted climate), so there is no
window-size correction to choose and no normal-theory assumption. block p = empirical
two-sided tail; BH-FDR across blocks.

NOTE this does NOT remove the N=19 climate–structure confounding (λ≈4); it gives the
correctly size-calibrated block test and the honest polygenic-enrichment readout.

Inputs: pipelineB/traj_{p,v,sites}.npy + hapfreq_registry.csv + lib.load_climate().
Outputs: pipelineB/block_wza.csv, manhattan_wza_bio1.png.  Env: kmate; run via sbatch.
"""
import os, sys
import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lib

H = os.environ.get("HF_DIR", "results/grenenet_gea/hapfreq")   # hapfreq store (clq90: hapfreq_clq90)
B = os.environ.get("PB_DIR", f"{H}/pipelineB")
EPS = 1e-3
N_PERM = int(os.environ.get("N_PERM", 5000))
CLIM = "bio1"
Tg = np.array([0.0, 1.0, 2.0, 3.0])


def logit(p):
    p = np.clip(p, EPS, 1 - EPS)
    return np.log(p / (1 - p))


def main():
    P = np.load(f"{B}/traj_p.npy").astype(np.float64)
    V = np.load(f"{B}/traj_v.npy").astype(np.float64)
    sites = np.load(f"{B}/traj_sites.npy")
    reg = pd.read_csv(f"{H}/hapfreq_registry.csv")
    S, _, nhap = P.shape

    c = lib.load_climate()[CLIM].reindex([int(s) for s in sites]).to_numpy(float)
    c = (c - c.mean()) / c.std(ddof=1)

    # Step 3 (s_g, SE) + Step 4 weights (same as build_hap_gea.py)
    y = logit(P); w = 1.0 / np.clip(V, 1e-12, None)
    sw = w.sum(1); tbar = (w * Tg[None, :, None]).sum(1) / sw
    dt = Tg[None, :, None] - tbar[:, None, :]
    Sxx = (w * dt**2).sum(1)
    s_g = (w * dt * y).sum(1) / Sxx
    se_g = np.sqrt(1.0 / Sxx)
    w0 = 1.0 / np.clip(se_g**2, 1e-12, None); sw0 = w0.sum(0)
    sbar0 = (w0 * s_g).sum(0) / sw0
    Q = (w0 * (s_g - sbar0[None, :])**2).sum(0)
    denom = sw0 - (w0**2).sum(0) / sw0
    tau2 = np.clip((Q - (S - 1)) / np.clip(denom, 1e-12, None), 0, None)
    wg = 1.0 / (se_g**2 + tau2[None, :])

    def beta1(cvec):
        swg = wg.sum(0)
        cbar = (wg * cvec[:, None]).sum(0) / swg
        sbar = (wg * s_g).sum(0) / swg
        dc = cvec[:, None] - cbar[None, :]
        return (wg * dc * (s_g - sbar[None, :])).sum(0) / np.clip((wg * dc**2).sum(0), 1e-12, None)

    # ---- block membership: kept = (k-1) independent one-vs-rest haps (drop max panel_freq) ----
    reg["unit"] = reg.chrom + ":" + reg.unit_start.astype(str) + "-" + reg.unit_end.astype(str)
    testable = (reg.covered & reg.panel_freq.between(0.05, 0.95)).to_numpy()
    keep = np.zeros(nhap, bool)
    for u, grp in reg[testable].groupby("unit"):
        kept = grp.index.drop(grp.panel_freq.idxmax()) if len(grp) > 1 else grp.index
        keep[kept.to_numpy()] = True
    ki = np.where(keep)[0]                                       # kept haplotype indices
    units = reg.loc[ki, "unit"].to_numpy()
    ucodes, binv = np.unique(units, return_inverse=True)        # block code per kept hap
    nblk = len(ucodes)
    wi = (reg.loc[ki, "panel_freq"].to_numpy() * (1 - reg.loc[ki, "panel_freq"].to_numpy()))
    den = np.sqrt(np.bincount(binv, weights=wi**2, minlength=nblk))   # constant Σw² per block

    def wza(cvec):
        b = beta1(cvec)[ki]
        med = np.median(b); mad = np.median(np.abs(b - med)) * 1.4826 + 1e-12
        z = (b - med) / mad                                     # robust per-hap Z
        num = np.bincount(binv, weights=wi * z, minlength=nblk)
        return num / np.clip(den, 1e-12, None)

    wza_obs = wza(c)
    rng = np.random.RandomState(0)
    cnt = np.zeros(nblk)
    aobs = np.abs(wza_obs)
    for _ in range(N_PERM):
        cnt += (np.abs(wza(c[rng.permutation(S)])) >= aobs)
    p_two = (cnt + 1) / (N_PERM + 1)

    # per-block summary: best (most extreme |z|) kept haplotype as the lead
    bdf = reg.loc[ki, ["unit", "chrom", "unit_start", "unit_end", "cluster", "panel_freq"]].copy()
    bdf["beta1"] = beta1(c)[ki]
    lead = bdf.loc[bdf.groupby("unit").beta1.apply(lambda s: s.abs().idxmax())]
    blk = pd.DataFrame({"unit": ucodes, "wza": wza_obs, "block_p": p_two,
                        "n_hap": np.bincount(binv, minlength=nblk)})
    blk = blk.merge(lead.rename(columns={"cluster": "lead_cluster", "beta1": "lead_beta1"})
                        [["unit", "chrom", "unit_start", "unit_end", "lead_cluster", "lead_beta1"]],
                    on="unit")
    m = len(blk); order = blk.block_p.values.argsort()
    q = np.empty(m); ranked = blk.block_p.values[order] * m / (np.arange(m) + 1)
    q[order] = np.minimum.accumulate(ranked[::-1])[::-1]
    blk["q"] = np.clip(q, 0, 1)
    blk = blk.sort_values("block_p").reset_index(drop=True)
    blk.to_csv(f"{B}/block_wza.csv", index=False)

    p = blk.block_p.clip(1e-6, 1).values
    lam = np.median(stats.chi2.isf(p, 1)) / stats.chi2.ppf(0.5, 1)
    print(f"[done] block-WZA {CLIM}: {nblk:,} blocks, {len(ki):,} kept haps, {N_PERM} perms")
    print(f"  λ(block-WZA) = {lam:.3f}")
    print(f"  q<0.05: {int((blk.q<0.05).sum())}; q<0.10: {int((blk.q<0.10).sum())}; "
          f"raw p<0.01: {int((blk.block_p<0.01).sum())} (null≈{int(0.01*m)})")
    print(blk.head(12).to_string(index=False))


if __name__ == "__main__":
    main()
