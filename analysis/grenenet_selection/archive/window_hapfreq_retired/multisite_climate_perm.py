#!/usr/bin/env python
"""STAGE 2 — site-climate permutation null for the multisite founder-GWAS CLIMATE contrast.

The CLIMATE contrast (climate-differential / local-adaptation selection) is the only one a
site->climate label permutation tests: GLOBAL (1-vector) and JOINT (z^T C^-1 z) are invariant to
which site is which, so permuting bio1 across sites is the exact exchangeable null for "does a
block's cross-site selection effect track climate". STAGE 1 reported parametric N(0,1) p-values
(lambda_clim ~ 1); this calibrates them genome-wide without re-running EMMAX.

Loads multisite_founder_gwas.npz (Z = M x S genomic-controlled effects, C cross-site covariance,
bio1, marker coords). For each of N_PERM permutations of bio1 across sites, recompute the
C-orthogonalized climate contrast and its per-block z. Reports:
  (1) FWER (max-statistic): null dist of max_b |clim_z|; genome-wide-significant blocks = those whose
      observed |clim_z| exceeds the (1-alpha) quantile of the null max. FWER p of the top block.
  (2) GLOBAL EXCESS: is there ANY diffuse climate signal? statistic = #blocks with |clim_z|>2.5
      (and var of clim_z); observed vs permutation null.
  (3) cross_site_winners PC1-vs-bio1 r: permutation p (permute bio1 vs the derived winning-clade PC1).
Writes multisite_climate_perm.{json,png}. Env: kmate. N_PERM (10000), PERM_SEED (0).
"""
import os, sys, json
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

H = "analysis/grenenet_gea/archive/window_hapfreq_retired/hapfreq"
N_PERM = int(os.environ.get("N_PERM", 10000))
SEED = int(os.environ.get("PERM_SEED", 0))
THR = 2.5                                                    # |z| cutoff for the diffuse-excess count
SUFFIX = os.environ.get("OUT_SUFFIX", "")    # "_clq90" pairs with the clq0.9 GWAS run


def clim_contrast(ZC, Cinv, one, dg, b):
    """C-orthogonalized (vs GLOBAL) standardized-bio1 contrast -> per-block z. Matches STAGE 1."""
    c0 = (b - b.mean()) / b.std()
    c = c0 - (float(one @ Cinv @ c0) / dg) * one             # orthogonalize against the 1-vector
    dc = float(c @ Cinv @ c)
    return (ZC @ c) / np.sqrt(dc)                            # ZC = Z @ Cinv precomputed (M x S)


def main():
    d = np.load(f"{H}/multisite_founder_gwas{SUFFIX}.npz", allow_pickle=True)
    Z, C, bio1 = d["Z"], d["C"], d["bio1"]
    M, S = Z.shape
    Cinv = np.linalg.pinv(C); one = np.ones(S); dg = float(one @ Cinv @ one)
    ZC = Z @ Cinv                                            # (M x S), constant across permutations
    obs = clim_contrast(ZC, Cinv, one, dg, bio1)

    # sanity: reproduce STAGE 1 z_clim value-set (csv rows are reordered by p_joint and (chrom,start)
    # is non-unique across allele-markers, so compare sorted value sets, not row-paired).
    csv = pd.read_csv(f"{H}/multisite_founder_gwas{SUFFIX}.csv")
    rerr = float(np.max(np.abs(np.sort(obs) - np.sort(csv["z_clim"].to_numpy()))))
    print(f"loaded {M:,} blocks x {S} sites | reproduce STAGE1 z_clim value-set max|err|={rerr:.2e}")

    obs_max = float(np.abs(obs).max())
    obs_nexc = int((np.abs(obs) > THR).sum())
    obs_var = float(obs.var())

    rng = np.random.default_rng(SEED)
    null_max = np.empty(N_PERM); null_nexc = np.empty(N_PERM, int); null_var = np.empty(N_PERM)
    for k in range(N_PERM):
        zc = clim_contrast(ZC, Cinv, one, dg, rng.permutation(bio1))
        a = np.abs(zc)
        null_max[k] = a.max(); null_nexc[k] = int((a > THR).sum()); null_var[k] = zc.var()

    fwer_thr = float(np.quantile(null_max, 0.95))
    n_fwer = int((np.abs(obs) > fwer_thr).sum())
    p_top = float((null_max >= obs_max).mean())              # FWER p of the strongest block
    p_excess = float((null_nexc >= obs_nexc).mean())         # diffuse-signal p
    p_var = float((null_var >= obs_var).mean())

    # (3) cross_site_winners PC1 vs bio1 permutation p (independent analysis, 17 sites)
    cs_p1 = cs_r1 = None
    csw_path = f"{H}/cross_site_winners_30{SUFFIX}.csv"
    if os.path.exists(csw_path):
        csw = pd.read_csv(csw_path); m = np.isfinite(csw.bio1) & np.isfinite(csw.win_pc1)
        x, y = csw.bio1[m].to_numpy(), csw.win_pc1[m].to_numpy()
        cs_r1 = float(np.corrcoef(x, y)[0, 1])
        rng2 = np.random.default_rng(SEED)
        nd = np.array([abs(np.corrcoef(rng2.permutation(x), y)[0, 1]) for _ in range(N_PERM)])
        cs_p1 = float((nd >= abs(cs_r1)).mean())

    print(f"\nCLIMATE contrast permutation null (N_PERM={N_PERM}, {S} sites):")
    print(f"  (1) FWER: obs max|z|={obs_max:.2f}; null-max 95% thr={fwer_thr:.2f} -> "
          f"{n_fwer} genome-wide-sig blocks; FWER p(top)={p_top:.4f}")
    print(f"  (2) EXCESS: #|z|>{THR} obs={obs_nexc} vs null median {np.median(null_nexc):.0f} -> p={p_excess:.4f}; "
          f"var(z) obs={obs_var:.3f} vs null {np.median(null_var):.3f} -> p={p_var:.4f}")
    if cs_r1 is not None:
        print(f"  (3) cross_site winning-clade PC1 vs bio1: r={cs_r1:.2f}, perm p={cs_p1:.4f}")

    out = dict(n_perm=N_PERM, seed=SEED, n_blocks=M, n_sites=S, thr=THR,
               obs_max_z=obs_max, fwer_thr_95=fwer_thr, n_fwer_sig=n_fwer, fwer_p_top=p_top,
               obs_n_excess=obs_nexc, null_n_excess_med=float(np.median(null_nexc)), p_excess=p_excess,
               obs_var_z=obs_var, null_var_z_med=float(np.median(null_var)), p_var=p_var,
               crosssite_pc1_bio1_r=cs_r1, crosssite_pc1_bio1_perm_p=cs_p1,
               note="site->climate permutation calibrates ONLY the CLIMATE contrast (GLOBAL/JOINT "
                    "invariant to site labels). Does NOT remove the variable-trait-window x climate "
                    "structural confound (use TRAIT_GENS=1 run to check that).")
    json.dump(out, open(f"{H}/multisite_climate_perm{SUFFIX}.json", "w"), indent=2)

    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    ax[0].hist(null_max, bins=50, color="#bcc7e6", edgecolor="none")
    ax[0].axvline(obs_max, color="firebrick", lw=1.6, label=f"observed max|z|={obs_max:.2f}")
    ax[0].axvline(fwer_thr, color="k", lw=1, ls="--", label=f"null 95% = {fwer_thr:.2f}")
    ax[0].set_xlabel("genome-wide max |climate z| (per permutation)"); ax[0].set_ylabel("permutations")
    ax[0].set_title(f"(1) FWER max-statistic null\nFWER p(top block) = {p_top:.3f}", fontsize=10, loc="left")
    ax[0].legend(fontsize=8, frameon=False); ax[0].spines[["top", "right"]].set_visible(False)
    ax[1].hist(null_nexc, bins=40, color="#bcc7e6", edgecolor="none")
    ax[1].axvline(obs_nexc, color="firebrick", lw=1.6, label=f"observed = {obs_nexc}")
    ax[1].set_xlabel(f"# blocks with |climate z| > {THR}"); ax[1].set_ylabel("permutations")
    ax[1].set_title(f"(2) diffuse climate-signal excess\np = {p_excess:.3f}", fontsize=10, loc="left")
    ax[1].legend(fontsize=8, frameon=False); ax[1].spines[["top", "right"]].set_visible(False)
    fig.suptitle(f"Multisite CLIMATE contrast — site-climate permutation null ({S} sites, {N_PERM} perms)", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.94]); fig.savefig(f"{H}/multisite_climate_perm{SUFFIX}.png", dpi=150)
    print(f"\n[done] {H}/multisite_climate_perm{SUFFIX}.json + .png")


if __name__ == "__main__":
    main()
