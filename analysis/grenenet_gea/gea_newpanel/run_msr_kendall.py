#!/usr/bin/env python
"""Structure-preserving-null rank climate-GEA (SNP vs non-SNP), site level.

GOAL: keep the monotone RANK statistic (Spearman rho, the tractable analog of the
already-computed site Kendall tau-b; both are monotone-equivalent at n=31) but
control population structure via the NULL DISTRIBUTION rather than by shrinking the
statistic. p-values then stop being inflated yet real signal survives. We then ask
whether the surviving signal is the SAME or DIFFERENT between SNP and non-SNP.

Unit = 31 flower-weighted SITE means (one climate value = one independent obs),
exactly as run_site_binom_kendall.py / build_site_matrices.py.

Structure model: 31x31 among-site covariance Sigma (BayPass-Omega analog) estimated
from the SNP site-AF matrix (structure is structure; the SNP-derived Sigma is used as
the reference for BOTH classes). Sigma = (1/M) Xc @ Xc^T with Xc = SNP site_af
standardized per variant across the 31 sites.  Eigendecompose Sigma = V diag(lam) V^T.

Structure-preserving NULL = MORAN SPECTRAL RANDOMIZATION (sign-flip variant, Wagner &
Dray 2015): express the observed (centered) site-env in the Sigma eigenbasis,
c = V^T e; each surrogate flips the sign of each spectral coefficient independently,
e_surr = V @ (s .* c), s in {-1,+1}^31.  This preserves EXACTLY the projection
magnitude |c_k| onto every structure eigen-direction (hence identical among-site
autocorrelation to the real bio1) while randomizing alignment with the genome.
(MVN(0,Sigma) only matches structure in expectation and injects sampling noise; the
sign-flip MSR is exact per draw, so we use it.)

Statistic (vectorizable across M variants): Spearman rho between the 31 site-mean AFs
and env = Pearson correlation of their average-ranks over the 31 sites.
Pre-rank+center+unit-norm the site_af rows once -> U_af [M x 31]; rank+center+unit-norm
each surrogate env -> E [31 x P]; then rho for ALL variants at once = U_af @ e.

Empirical structure-null p (two-sided, +1 smoothed): (1+#{|rho_null|>=|rho_obs|})/(1+P).
P set to 10000 so the empirical-p floor 1/(P+1)~=1e-4 can actually reach the
reported p<1e-4 threshold (naive P=2000 could not).

Outputs to analysis/grenenet_gea/gea_newpanel/results/msr_kendall/:
  msr_kendall_{cls}_site_{clim}.csv : chrom,pos,ref_len,alt_len,site_maf,block,
                                       rho,p_naive,p_msr
  msr_arrays_{cls}_{clim}.npz       : rho,p_naive,p_msr,site_maf,gif_* for plotting
  sigma_{clim}.npz                  : Sigma, eigvals, eigvecs, env, env_rank
"""
from __future__ import annotations
import argparse, os, sys, time
import numpy as np
import pandas as pd
from scipy.stats import rankdata, chi2, t as tdist
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import lib

NSITES = 31


def gif_from_p(p):
    p = np.asarray(p, float)
    p = p[np.isfinite(p)]
    p = np.clip(p, 1e-300, 1.0)
    obs_chi2 = chi2.isf(p, 1)                       # observed 1-df chi2
    return float(np.median(obs_chi2) / chi2.ppf(0.5, 1))


def bh_fdr(p):
    p = np.asarray(p, float)
    n = p.size
    order = np.argsort(p)
    ranked = p[order] * n / (np.arange(n) + 1)
    q = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.clip(q, 0, 1)
    return out


def build_site_af(cmdir, cls, gen, sites, site_of, w):
    """Flower-weighted site-mean AF [n_site x M] + record table."""
    recs = pd.read_csv(f"{cmdir}/{cls}_gen{gen}.records.csv")
    af = np.load(f"{cmdir}/{cls}_gen{gen}_af.npy")          # [355 x M] float32
    assert af.shape[1] == len(recs), (af.shape, len(recs))
    site_af = np.empty((len(sites), af.shape[1]), dtype=np.float64)
    for i, s in enumerate(sites):
        m = site_of == s
        ws = w[m]; wn = ws / ws.sum()
        site_af[i] = wn @ af[m].astype(np.float64)
    del af
    return site_af, recs


def rank_unit_cols(site_af):
    """[n_site x M] AF -> [M x n_site] avg-rank, per variant centered + unit-norm."""
    R = rankdata(site_af, method="average", axis=0)        # rank each variant / sites
    R = R - R.mean(axis=0, keepdims=True)
    nrm = np.sqrt((R * R).sum(axis=0))
    nrm[nrm == 0] = np.inf
    U = (R / nrm[None, :]).T                                # [M x n_site]
    return np.ascontiguousarray(U.astype(np.float32))


def rank_unit_vec(x):
    r = rankdata(x, method="average").astype(np.float64)
    r = r - r.mean()
    r = r / np.sqrt((r * r).sum())
    return r


def naive_spearman_p(rho, n):
    """Analytic two-sided Spearman p via t-approx (n-2 df), vectorized."""
    rho = np.clip(rho, -0.999999, 0.999999)
    tstat = rho * np.sqrt((n - 2) / (1.0 - rho * rho))
    return 2.0 * tdist.sf(np.abs(tstat), n - 2)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--classes", default="snp,nonsnp")
    ap.add_argument("--gen", type=int, default=9)
    ap.add_argument("--climate", default="bio1")
    ap.add_argument("--P", type=int, default=10000, help="# MSR surrogate envs")
    ap.add_argument("--sigma-class", default="snp", help="class used to build Sigma")
    ap.add_argument("--seed", type=int, default=12345)
    ap.add_argument("--cmdir", default=f"{lib.GEA}/phase1_replication/results/class_matrices")
    ap.add_argument("--out", default=f"{lib.GEA}/gea_newpanel/results/msr_kendall")
    ap.add_argument("--batch", type=int, default=500, help="surrogate batch size")
    args = ap.parse_args()
    classes = args.classes.split(",")
    os.makedirs(args.out, exist_ok=True)
    rng = np.random.default_rng(args.seed)

    pools = pd.read_csv(f"{args.cmdir}/gen{args.gen}.pools.csv")
    sites = np.sort(pools.site.unique())
    assert len(sites) == NSITES, len(sites)
    site_of = pools.site.to_numpy()
    w = pools.total_flowers.to_numpy(float)
    clim = (pools.drop_duplicates("site").set_index("site")
            .loc[sites, args.climate].to_numpy(float))
    print(f"[setup] {NSITES} sites x classes {classes}; env {args.climate} "
          f"range {clim.min():.1f}-{clim.max():.1f}; P={args.P}", flush=True)

    # --- Sigma (among-site structure covariance) from the reference class ---
    t0 = time.time()
    ref_af, _ = build_site_af(args.cmdir, args.sigma_class, args.gen, sites, site_of, w)
    Xc = ref_af - ref_af.mean(axis=0, keepdims=True)         # center per variant
    sd = Xc.std(axis=0, ddof=0)
    keep = sd > 0
    Xc = Xc[:, keep] / sd[keep][None, :]                     # standardize per variant
    Sigma = (Xc @ Xc.T) / Xc.shape[1]                        # [31 x 31]
    lam, V = np.linalg.eigh(Sigma)                           # ascending
    lam = lam[::-1]; V = V[:, ::-1]
    print(f"[Sigma from {args.sigma_class}] M={Xc.shape[1]:,} "
          f"eig top5={np.round(lam[:5],3)} tail={np.round(lam[-3:],4)} "
          f"({time.time()-t0:.0f}s)", flush=True)
    del ref_af, Xc

    # --- observed + MSR surrogate envs, rank-normalized ---
    e = clim - clim.mean()
    c = V.T @ e                                              # spectral coeffs
    u_obs = rank_unit_vec(e)                                 # ranked obs env
    signs = rng.choice([-1.0, 1.0], size=(args.P, NSITES))
    E_surr = (V @ (signs * c[None, :]).T)                    # [31 x P] continuous
    # rank + center + unit-norm each surrogate column
    Er = rankdata(E_surr, method="average", axis=0)
    Er = Er - Er.mean(axis=0, keepdims=True)
    Er = Er / np.sqrt((Er * Er).sum(axis=0, keepdims=True))
    E = np.ascontiguousarray(Er.astype(np.float32))          # [31 x P]
    # sanity: surrogate autocorr wrt Sigma matches observed (MSR preserves it exactly)
    quad_obs = float(e @ Sigma @ e)
    qs = np.array([E_surr[:, k] @ Sigma @ E_surr[:, k] for k in range(min(args.P, 200))])
    print(f"[env] obs e'Sigma e={quad_obs:.4f}; surrogate mean={qs.mean():.4f} "
          f"(should match: MSR preserves structure)", flush=True)

    np.savez(f"{args.out}/sigma_{args.climate}.npz", Sigma=Sigma, eigvals=lam,
             eigvecs=V, env=e, env_rank=u_obs)

    summ = {}
    for cls in classes:
        t0 = time.time()
        site_af, recs = build_site_af(args.cmdir, cls, args.gen, sites, site_of, w)
        mean_af = site_af.mean(axis=0)
        site_maf = np.minimum(mean_af, 1.0 - mean_af)
        U_af = rank_unit_cols(site_af)                       # [M x 31] f32
        del site_af
        M = U_af.shape[0]
        rho_obs = U_af @ u_obs.astype(np.float32)            # [M]
        arho = np.abs(rho_obs)
        p_naive = naive_spearman_p(rho_obs.astype(np.float64), NSITES)
        # empirical two-sided structure-null p, batched over surrogates
        ge = np.zeros(M, dtype=np.int32)
        for b0 in range(0, args.P, args.batch):
            b1 = min(b0 + args.batch, args.P)
            null_block = U_af @ E[:, b0:b1]                  # [M x b]
            ge += (np.abs(null_block) >= arho[:, None]).sum(axis=1).astype(np.int32)
        p_msr = (1 + ge) / (1 + args.P)
        del U_af

        gif_naive = gif_from_p(p_naive)
        gif_msr = gif_from_p(p_msr)
        q_msr = bh_fdr(p_msr)
        n_p1e4 = int((p_msr < 1e-4).sum())
        n_fdr = int((q_msr < 0.05).sum())
        print(f"[{cls}] M={M:,} | GIF naive={gif_naive:.2f} MSR={gif_msr:.2f} | "
              f"p_msr<1e-4: {n_p1e4:,} | FDR q<.05: {n_fdr:,} | "
              f"min p_msr={p_msr.min():.2e} min q={q_msr.min():.3f} "
              f"({time.time()-t0:.0f}s)", flush=True)

        out = recs[["chrom", "pos", "ref_len", "alt_len", "block"]].copy()
        out["site_maf"] = site_maf
        out["rho"] = rho_obs
        out["p_naive"] = p_naive
        out["p_msr"] = p_msr
        csv = f"{args.out}/msr_kendall_{cls}_site_{args.climate}.csv"
        out.to_csv(csv, index=False)
        np.savez(f"{args.out}/msr_arrays_{cls}_{args.climate}.npz",
                 rho=rho_obs, p_naive=p_naive, p_msr=p_msr, q_msr=q_msr,
                 site_maf=site_maf, gif_naive=gif_naive, gif_msr=gif_msr,
                 chrom=out.chrom.to_numpy().astype(str), pos=out.pos.to_numpy(),
                 block=out.block.to_numpy().astype(str))
        summ[cls] = dict(M=M, gif_naive=gif_naive, gif_msr=gif_msr,
                         n_p1e4=n_p1e4, n_fdr=n_fdr,
                         min_p=float(p_msr.min()), min_q=float(q_msr.min()))
        print(f"  -> {csv}", flush=True)

    print("\n=== MSR STRUCTURE-NULL SUMMARY ===")
    for cls, d in summ.items():
        print(f"{cls:8s} M={d['M']:>9,} GIFnaive={d['gif_naive']:.2f} "
              f"GIFmsr={d['gif_msr']:.2f} p<1e-4={d['n_p1e4']:>6,} "
              f"FDRq<.05={d['n_fdr']:>6,} minP={d['min_p']:.1e} minQ={d['min_q']:.3f}")


if __name__ == "__main__":
    main()
