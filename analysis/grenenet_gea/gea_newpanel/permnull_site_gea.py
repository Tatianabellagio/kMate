#!/usr/bin/env python
"""Site-permutation empirical null for the site-level climate GEA.

HYPOTHESIS under test: the GIF/genomic-control calibration used in the site-level
LFMM is OVER-conservative, hiding real climate-associated loci. We replace GC with a
proper SITE-PERMUTATION empirical null (exchangeable unit = the 31 sites) and ask how
many loci/blocks become genome-wide significant.

Method (fast, structure-corrected, does NOT refit LFMM per perm):
  1. Y = [31 sites x n_var] Delta-p matrix (per class). Center Y per variant over
     sites. SVD over sites -> top-K left singular vectors U_K (site-level PCs =
     LFMM latent-factor analog). Residualize: Yr = Yc - U_K (U_K^T Yc). Unit-norm
     each variant column of Yr.
  2. Standardize climate X over 31 sites (mean 0, unit norm). Per-locus statistic =
     cosine(Yr_locus, X) = partial correlation controlling for the K site-PCs
     (numerator identical whether X is residualized on U_K or not, since Yr_locus is
     orthogonal to span(U_K)).
  3. Permute the 31-site X vector (>=NPERM). Per perm recompute all |stat|:
       (a) genome-wide MAX |stat| -> 5% of the max-null = permutation FWER threshold;
       (b) per-locus empirical p = (1 + #{|stat_perm| >= |stat_obs|}) / (1 + NPERM).
     BH-FDR on the per-locus empirical p.
  4. site-MAF >= 0.05 filter (site_af = Delta-p + p0[col], mean over sites, folded).
     Assign clq0.9 blocks. Annotate significant blocks with genes.

Runs all {class} x {axis} x {K} combos. Outputs under
results/grenenet_gea/gea_newpanel/permnull/.
"""
from __future__ import annotations
import os, sys, json, time
import numpy as np
import pandas as pd

sys.path.insert(0, "/global/scratch/users/tbellg/kmate/analysis/grenenet_gea")
sys.path.insert(0, "/global/scratch/users/tbellg/kmate/analysis/grenenet_gea/gea_newpanel")
from blocks_clq09 import assign_clq09_blocks
from lib import annotate_svs, load_genes

ROOT = "/global/scratch/users/tbellg/kmate/results/grenenet_gea"
LFMM = f"{ROOT}/gea_newpanel/lfmm_site"
ENVD = f"{ROOT}/gea_newpanel/env_site"
RECD = f"{ROOT}/phase1_replication/class_matrices"
OUTD = f"{ROOT}/gea_newpanel/permnull"
os.makedirs(OUTD, exist_ok=True)

NSITES = 31
NPERM = int(os.environ.get("PERMNULL_NPERM", "2000"))
KS = [int(k) for k in os.environ.get("PERMNULL_KS", "1,3").split(",")]
AXES = os.environ.get("PERMNULL_AXES", "bio5,pc1").split(",")
CLASSES = os.environ.get("PERMNULL_CLASSES", "snp,nonsnp").split(",")
MAF_MIN = 0.05
SEED = 12345


def bh_fdr(p):
    p = np.asarray(p, float)
    n = p.size
    order = np.argsort(p)
    ranked = p[order] * n / (np.arange(n) + 1)
    q = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.clip(q, 0, 1)
    return out


def load_class(cls):
    dims = open(f"{LFMM}/lfmm_{cls}_site_dims.txt").read().split()
    ns, nv = int(dims[0]), int(dims[1])
    assert ns == NSITES, (ns, NSITES)
    Y = np.fromfile(f"{LFMM}/lfmm_{cls}_site_Y.f64", dtype=np.float64).reshape(ns, nv)
    rec = pd.read_csv(f"{RECD}/{cls}_gen9.records.csv")
    assert len(rec) == nv, (len(rec), nv)
    p0 = np.load(f"{ROOT}/af_store/p0_{cls}.npy").astype(np.float64)
    col = rec["col"].to_numpy()
    site_af = Y + p0[col][None, :]                 # [31 x nv] per-site AF
    mean_af = site_af.mean(axis=0)
    maf = np.minimum(mean_af, 1.0 - mean_af)
    return Y, rec, maf


def residualize(Y, K):
    """Center per variant over sites, remove top-K site-PCs, unit-norm columns."""
    Yc = Y - Y.mean(axis=0, keepdims=True)
    # SVD over sites: U is [nsites x nsites]; left singular vectors = site-PC basis
    # economy SVD on Yc (31 x nv) -> U (31 x 31), no need for Vt
    U, S, _ = np.linalg.svd(Yc, full_matrices=False)
    Uk = U[:, :K]                                   # [31 x K]
    Yr = Yc - Uk @ (Uk.T @ Yc)                      # residual, orthogonal to Uk
    nrm = np.sqrt((Yr * Yr).sum(axis=0))
    nrm[nrm == 0] = np.inf
    Yr_unit = Yr / nrm[None, :]
    return Yr_unit


def std_env(axis):
    x = pd.read_csv(f"{ENVD}/env_site_{axis}.csv").iloc[:, 0].to_numpy(float)
    assert x.size == NSITES
    x = x - x.mean()
    x = x / np.sqrt((x * x).sum())
    return x


def run_combo(Yr_m, x, rng):
    """Yr_m = unit-norm residual Y restricted to the tested (MAF-passing) loci.
    Returns obs |stat|, per-perm genome-wide MAX over the tested family, per-locus
    empirical p. FWER max-null is over exactly the tested loci."""
    nv = Yr_m.shape[1]
    stat_obs = np.abs(Yr_m.T @ x)                   # partial corr, |.|
    ge = np.zeros(nv, dtype=np.int32)               # #{|perm| >= |obs|}
    maxnull = np.empty(NPERM)
    for b in range(NPERM):
        xp = x[rng.permutation(NSITES)]
        sp = np.abs(Yr_m.T @ xp)
        maxnull[b] = sp.max()
        ge += (sp >= stat_obs)
    emp_p = (1 + ge) / (1 + NPERM)
    return stat_obs, maxnull, emp_p


def main():
    genes = load_genes()
    summary = []
    for cls in CLASSES:
        t0 = time.time()
        Y, rec, maf = load_class(cls)
        maf_mask = maf >= MAF_MIN
        # block assignment (once per class)
        blocks = assign_clq09_blocks(rec["chrom"].to_numpy(), rec["pos"].to_numpy())
        rec = rec.assign(site_maf=maf, block_clq09=blocks)
        print(f"[{cls}] loaded Y{Y.shape} maf>=.05 keeps {maf_mask.sum():,}/{len(maf):,} "
              f"({time.time()-t0:.0f}s)", flush=True)
        idx = np.where(maf_mask)[0]
        sub_base = rec.iloc[idx].reset_index(drop=True)
        for K in KS:
            Yr = residualize(Y, K)
            Yr_m = np.ascontiguousarray(Yr[:, idx])   # tested family only
            for axis in AXES:
                rng = np.random.default_rng(SEED + K * 100 + hash(axis) % 97 + (cls == "snp"))
                x = std_env(axis)
                obs_m, maxnull, emp_p_m = run_combo(Yr_m, x, rng)
                thr_fwer = np.quantile(maxnull, 0.95)   # over the tested family
                q = bh_fdr(emp_p_m)
                pass_fwer = obs_m >= thr_fwer
                pass_fdr = q < 0.05
                # assemble hit table
                hit_mask = pass_fwer | pass_fdr
                sub = sub_base.copy()
                sub["stat"] = obs_m
                sub["emp_p"] = emp_p_m
                sub["fdr_q"] = q
                sub["pass_fwer"] = pass_fwer
                sub["pass_fdr"] = pass_fdr
                hits = sub[hit_mask].copy()
                if len(hits):
                    hits = annotate_svs(hits, flank=2000, genes=genes)
                tag = f"{cls}_{axis}_K{K}"
                hits.sort_values("emp_p").to_csv(f"{OUTD}/hits_{tag}.csv", index=False)
                # calibration: bulk emp_p uniformity (fraction < 0.05 should be ~0.05)
                frac_p05 = float((emp_p_m < 0.05).mean())
                nblk_fwer = sub.loc[pass_fwer, "block_clq09"].nunique()
                nblk_fdr = sub.loc[pass_fdr, "block_clq09"].nunique()
                row = dict(cls=cls, axis=axis, K=K, n_loci=int(len(idx)),
                           thr_fwer=float(thr_fwer), max_obs_stat=float(obs_m.max()),
                           n_fwer=int(pass_fwer.sum()), n_fdr=int(pass_fdr.sum()),
                           n_blocks_fwer=int(nblk_fwer), n_blocks_fdr=int(nblk_fdr),
                           min_emp_p=float(emp_p_m.min()), min_fdr_q=float(q.min()),
                           frac_emp_p_lt05=frac_p05)
                summary.append(row)
                print(f"  {tag}: n={len(idx):,} thrFWER={thr_fwer:.4f} "
                      f"maxObs={obs_m.max():.4f} nFWER={pass_fwer.sum()} "
                      f"nFDR={pass_fdr.sum()} minP={emp_p_m.min():.2e} "
                      f"minQ={q.min():.3f} frac_p<.05={frac_p05:.3f}", flush=True)
        del Y, Yr
    sdf = pd.DataFrame(summary)
    sdf.to_csv(f"{OUTD}/summary.csv", index=False)
    with open(f"{OUTD}/summary.json", "w") as f:
        json.dump(summary, f, indent=2)
    print("\n=== SUMMARY ===")
    print(sdf.to_string(index=False))


if __name__ == "__main__":
    main()
