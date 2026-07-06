#!/usr/bin/env python
"""Pipeline B v2 — STEP 3 + STEP 4 on haplotype trajectories (the climate-GEA).

STEP 3 (per site g, per haplotype): weighted log-odds slope of p̄_{t,g} vs gen,
weights ω_{t,g}=1/V_{t,g} over t∈{0,1,2,3} → selection coefficient s_g + SE(s_g),
both from the one weighted fit (PIPELINE_B_POOLED_MODEL.md Step 3).

STEP 4 (per haplotype, across the 19 sites): random-effects meta-regression
s_g = β0 + β1·bio1_g, weights w_g = 1/(SE(s_g)²+τ̂²), τ̂² by DerSimonian–Laird →
β1 = climate effect. Significance by SITE-LEVEL PERMUTATION: shuffle the bio1
labels across the 19 sites, refit β1 (fixed weights), build the null → two-sided
p (any association) + one-sided p (up-in-warm).

BLOCK-LEVEL TEST UNIT: each hapfreq column = a one-vs-rest allele. k=2 blocks keep
only the minor haplotype (the two columns are complements); k>2 keep all one-vs-rest
but correct multiple testing at BLOCK resolution: within-block Šidák over (k-1)
independent haplotypes → block p; BH-FDR across blocks.

Inputs: results/grenenet_gea/hapfreq/pipelineB/traj_{p,v,sites,nplots}.npy + registry.
Outputs (same dir): hap_gea.csv (per testable haplotype), block_gea.csv (per block, q).

Env: kmate. Permutation is the heavy step → run via sbatch.
"""
import os, sys
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lib

H = os.environ.get("HF_DIR", "results/grenenet_gea/hapfreq")   # hapfreq store (clq90: hapfreq_clq90)
B = os.environ.get("PB_DIR", f"{H}/pipelineB")
EPS = 1e-3
N_PERM = int(os.environ.get("N_PERM", 5000))
CLIM = "bio1"                                   # headline climate axis (mean annual temp)
T = np.array([0.0, 1.0, 2.0, 3.0])              # generations


def logit(p):
    p = np.clip(p, EPS, 1 - EPS)
    return np.log(p / (1 - p))


def main():
    P = np.load(f"{B}/traj_p.npy").astype(np.float64)           # (S,4,nhap)
    V = np.load(f"{B}/traj_v.npy").astype(np.float64)
    sites = np.load(f"{B}/traj_sites.npy")
    reg = pd.read_csv(f"{H}/hapfreq_registry.csv")
    S, _, nhap = P.shape

    # climate vector for the 19 sites, standardized
    clim_tab = lib.load_climate()[CLIM]
    c = clim_tab.reindex([int(s) for s in sites]).to_numpy(float)
    assert np.isfinite(c).all(), f"missing {CLIM} for some site: {sites[~np.isfinite(c)]}"
    c = (c - c.mean()) / c.std(ddof=1)

    # ---- STEP 3: per-site weighted log-odds slope s_g + SE ------------------------
    y = logit(P)                                               # (S,4,nhap)
    w = 1.0 / np.clip(V, 1e-12, None)                          # ω_{t,g}
    sw = w.sum(1)                                              # (S,nhap)  Σ_t ω
    tbar = (w * T[None, :, None]).sum(1) / sw                  # (S,nhap)  weighted mean gen
    dt = T[None, :, None] - tbar[:, None, :]                   # (S,4,nhap)
    Sxx = (w * dt**2).sum(1)                                   # (S,nhap)
    s_g = (w * dt * y).sum(1) / Sxx                            # (S,nhap)  selection coef
    se_g = np.sqrt(1.0 / Sxx)                                  # (S,nhap)  SE(s_g)

    # ---- STEP 4: DerSimonian–Laird τ², random-effects β1 -------------------------
    w0 = 1.0 / np.clip(se_g**2, 1e-12, None)                   # fixed-effect weights (S,nhap)
    sw0 = w0.sum(0)
    sbar0 = (w0 * s_g).sum(0) / sw0
    Q = (w0 * (s_g - sbar0[None, :])**2).sum(0)
    denom = sw0 - (w0**2).sum(0) / sw0
    tau2 = np.clip((Q - (S - 1)) / np.clip(denom, 1e-12, None), 0, None)   # (nhap,)
    wg = 1.0 / (se_g**2 + tau2[None, :])                       # random-effects weights (S,nhap)

    def beta1(cvec, wg, s_g):
        swg = wg.sum(0)
        cbar = (wg * cvec[:, None]).sum(0) / swg
        sbar = (wg * s_g).sum(0) / swg
        dc = cvec[:, None] - cbar[None, :]
        return (wg * dc * (s_g - sbar[None, :])).sum(0) / np.clip((wg * dc**2).sum(0), 1e-12, None)

    b_obs = beta1(c, wg, s_g)                                  # (nhap,)

    # site-level permutation null (fixed weights; permute bio1 across the 19 sites)
    rng = np.random.RandomState(0)
    cnt_two = np.zeros(nhap); cnt_one = np.zeros(nhap)
    ab = np.abs(b_obs)
    for _ in range(N_PERM):
        cp = c[rng.permutation(S)]
        bp = beta1(cp, wg, s_g)
        cnt_two += (np.abs(bp) >= ab)
        cnt_one += (bp >= b_obs)
    p_two = (cnt_two + 1) / (N_PERM + 1)
    p_one = (cnt_one + 1) / (N_PERM + 1)

    out = reg.copy()
    out["beta1"] = b_obs
    out["tau2"] = tau2
    out["s_mean"] = s_g.mean(0)                                # cross-site mean selection coef
    out["p_two"] = p_two
    out["p_one"] = p_one
    out.to_csv(f"{B}/hap_gea.csv", index=False)

    # ---- BLOCK-LEVEL test unit + FDR ---------------------------------------------
    # testable: covered units, panel_freq in [0.05,0.95]; k=2 keep the MINOR haplotype.
    out["unit"] = out.chrom + ":" + out.unit_start.astype(str) + "-" + out.unit_end.astype(str)
    testable = out.covered & out.panel_freq.between(0.05, 0.95)
    # Keep the k-1 INDEPENDENT one-vs-rest haplotypes per block: drop the single
    # highest-panel_freq member (the 'rest' reference). k=2 -> keep the minor only.
    # (Audit fix: the old code kept all k haps but used a (k-1) Šidák exponent in the
    # min -> anti-conservative for k>2; now #kept == exponent, exactly Šidák.)
    keep = pd.Series(False, index=out.index)
    for u, grp in out[testable].groupby("unit"):
        kept = grp.index.drop(grp.panel_freq.idxmax()) if len(grp) > 1 else grp.index
        keep.loc[kept] = True
    cand = out[keep].copy()
    nkept = cand.groupby("unit").cluster.transform("size")     # = independent haps tested
    cand["p_sidak"] = 1 - (1 - cand.p_two) ** nkept            # Šidák over the kept set
    blk = (cand.sort_values("p_sidak")
              .groupby("unit")
              .agg(chrom=("chrom", "first"), unit_start=("unit_start", "first"),
                   unit_end=("unit_end", "first"), n_hap_tested=("cluster", "size"),
                   best_cluster=("cluster", "first"), best_beta1=("beta1", "first"),
                   best_s_mean=("s_mean", "first"), block_p=("p_sidak", "min"))
              .reset_index())
    # BH-FDR across blocks
    m = len(blk); order = blk.block_p.values.argsort()
    q = np.empty(m); ranked = blk.block_p.values[order] * m / (np.arange(m) + 1)
    q[order] = np.minimum.accumulate(ranked[::-1])[::-1]
    blk["q"] = np.clip(q, 0, 1)
    blk = blk.sort_values("block_p").reset_index(drop=True)
    blk.to_csv(f"{B}/block_gea.csv", index=False)

    print(f"[done] {CLIM} GEA over {nhap:,} haplotypes, {S} sites, {N_PERM} perms")
    print(f"  testable haplotypes: {int(keep.sum()):,}; blocks tested: {m:,}")
    print(f"  blocks q<0.05: {int((blk.q<0.05).sum())}; q<0.10: {int((blk.q<0.10).sum())}; "
          f"raw block_p<0.01: {int((blk.block_p<0.01).sum())}")
    print(f"  top blocks:\n{blk.head(10).to_string(index=False)}")


if __name__ == "__main__":
    main()
