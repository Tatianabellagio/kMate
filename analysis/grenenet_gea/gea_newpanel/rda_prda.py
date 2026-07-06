#!/usr/bin/env python
"""Partial Redundancy Analysis (pRDA) climate GEA on the honest 31-site unit.

Multivariate, polygenic-oriented complement to the univariate site-level LFMM
(which found no permutation-defensible per-locus climate hits). Follows
Capblancq & Forester (2021) / Forester (2018):

  * RESPONSE  Y = [31 sites x nvar] site-level Delta-p (per class).
  * CONSTRAIN X = a compact multivariate climate set (temperature + precip axes).
  * CONDITION Z = a few site-level structure PCs (SVD of the Delta-p matrix;
                  the LFMM latent-factor analog) -> partial RDA.

Candidate detection = the `rdadapt` Mahalanobis method: loci loadings on the
climate-constrained RDA axes -> robust-free Mahalanobis distance D^2 ~ chi^2_K ->
genomic-control (GIF lambda) -> p-values -> BH-FDR. Mahalanobis distance using the
sample covariance of the loadings is AFFINE-INVARIANT, so vegan's loading-scaling
convention is irrelevant (a Python pRDA is exactly equivalent here).

Reports the genomic inflation factor (lambda) of the candidate p-values (aim ~1-1.2,
NOT deflated) and the clq0.9 blocks that are candidates at FDR. Annotates SV hits
with genes. Runs all classes x (structure-K) x (retained-axes-K) combos.

NOTE: permutation-null validation is intentionally NOT run here (already known null).

Outputs under results/grenenet_gea/gea_newpanel/rda/ .
"""
from __future__ import annotations
import os, sys, json, time
import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, "/global/scratch/users/tbellg/kmate/analysis/grenenet_gea")
sys.path.insert(0, "/global/scratch/users/tbellg/kmate/analysis/grenenet_gea/gea_newpanel")
from blocks_clq09 import assign_clq09_blocks
from lib import annotate_svs, load_genes

ROOT = "/global/scratch/users/tbellg/kmate/results/grenenet_gea"
LFMM = f"{ROOT}/gea_newpanel/lfmm_site"
ENVD = f"{ROOT}/gea_newpanel/env_site"
RECD = f"{ROOT}/phase1_replication/class_matrices"
OUTD = f"{ROOT}/gea_newpanel/rda"
os.makedirs(OUTD, exist_ok=True)

NSITES = 31
MAF_MIN = 0.05
# multivariate climate: max-temp, min-temp, annual precip, precip seasonality
CLIMSET = os.environ.get("RDA_CLIM", "bio5,bio6,bio12,bio15").split(",")
KSTRUCT = [int(k) for k in os.environ.get("RDA_KSTRUCT", "0,2,3").split(",")]  # cond PCs
CLASSES = os.environ.get("RDA_CLASSES", "snp,nonsnp,sv").split(",")
# p0 store: SV records index into the non-SNP store
P0_FOR = {"snp": "snp", "nonsnp": "nonsnp", "sv": "nonsnp"}


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
    p0 = np.load(f"{ROOT}/af_store/p0_{P0_FOR[cls]}.npy").astype(np.float64)
    col = rec["col"].to_numpy()
    site_af = Y + p0[col][None, :]          # per-site AF = Delta-p + baseline
    mean_af = site_af.mean(axis=0)
    maf = np.minimum(mean_af, 1.0 - mean_af)
    return Y, rec, maf


def load_climate():
    cols = []
    for a in CLIMSET:
        x = pd.read_csv(f"{ENVD}/env_site_{a}.csv").iloc[:, 0].to_numpy(float)
        assert x.size == NSITES
        cols.append(x)
    X = np.column_stack(cols)                # [31 x q]
    X = X - X.mean(axis=0, keepdims=True)
    X = X / X.std(axis=0, ddof=0, keepdims=True)
    return X


def struct_pcs(Yc, k):
    """Top-k site-level PCs (left singular vectors of centered Y). Orthonormal [31 x k]."""
    if k == 0:
        return np.zeros((NSITES, 0))
    U, S, _ = np.linalg.svd(Yc, full_matrices=False)
    return U[:, :k]


def prda_loadings(Yc, X, Z, kax):
    """Partial RDA. Yc = [31 x p] centered response. X=[31 x q] climate. Z=[31 x c]
    orthonormal condition PCs. Returns loci loadings on first kax constrained axes
    [p x kax] and the number of constrained axes available."""
    # residualize climate on the condition (Z orthonormal -> projection subtract)
    if Z.shape[1] > 0:
        Xr = X - Z @ (Z.T @ X)
        Yr = Yc - Z @ (Z.T @ Yc)
    else:
        Xr, Yr = X, Yc
    # hat matrix of climate residuals (pinv handles any collinearity in Xr)
    Hx = Xr @ np.linalg.pinv(Xr)             # [31 x 31]
    Yhat = Hx @ Yr                           # [31 x p] fitted (constrained) response
    U, S, Vt = np.linalg.svd(Yhat, full_matrices=False)
    rank = int((S > 1e-9 * S.max()).sum()) if S.size else 0
    k = min(kax, max(rank, 1))
    L = Vt[:k].T                             # [p x k] loci loadings (affine-invariant)
    return L, rank, k


def mahalanobis_pvals(L):
    """Capblancq rdadapt: Mahalanobis D^2 of loadings, GIF-corrected chi^2 p-values."""
    k = L.shape[1]
    mu = L.mean(axis=0)
    Lc = L - mu
    cov = (Lc.T @ Lc) / (L.shape[0] - 1)
    covinv = np.linalg.pinv(cov)
    D2 = np.einsum("ij,jk,ik->i", Lc, covinv, Lc)   # [p]
    lam = np.median(D2) / stats.chi2.ppf(0.5, k)     # genomic inflation factor
    p_raw = stats.chi2.sf(D2, k)
    p_gc = stats.chi2.sf(D2 / lam, k) if lam > 0 else p_raw
    return D2, lam, p_raw, p_gc


def main():
    genes = load_genes()
    X = load_climate()
    print(f"climate set = {CLIMSET}  X{X.shape}", flush=True)
    summary = []
    for cls in CLASSES:
        t0 = time.time()
        Y, rec, maf = load_class(cls)
        keep = maf >= MAF_MIN
        blocks = assign_clq09_blocks(rec["chrom"].to_numpy(), rec["pos"].to_numpy())
        rec = rec.assign(site_maf=maf, block_clq09=blocks)
        idx = np.where(keep)[0]
        sub_base = rec.iloc[idx].reset_index(drop=True)
        Yc = Y - Y.mean(axis=0, keepdims=True)
        Yc_m = np.ascontiguousarray(Yc[:, idx])
        print(f"[{cls}] Y{Y.shape} maf>=.05 keeps {keep.sum():,}/{len(maf):,} "
              f"({time.time()-t0:.0f}s)", flush=True)
        for c in KSTRUCT:
            Z = struct_pcs(Yc_m, c)
            L, rank, kax = prda_loadings(Yc_m, X, Z, kax=X.shape[1])
            D2, lam, p_raw, p_gc = mahalanobis_pvals(L)
            q = bh_fdr(p_gc)
            hit = q < 0.05
            sub = sub_base.copy()
            sub["D2"] = D2
            sub["p_raw"] = p_raw
            sub["p_gc"] = p_gc
            sub["fdr_q"] = q
            sub["cand"] = hit
            hits = sub[hit].copy()
            if len(hits):
                hits = annotate_svs(hits, flank=2000, genes=genes)
            tag = f"{cls}_c{c}"
            hits.sort_values("fdr_q").to_csv(f"{OUTD}/cand_{tag}.csv", index=False)
            nblk = sub.loc[hit, "block_clq09"].nunique()
            row = dict(cls=cls, kstruct=c, n_loci=int(len(idx)),
                       constr_rank=rank, kax=kax, lambda_gif=float(lam),
                       n_cand=int(hit.sum()), n_blocks_cand=int(nblk),
                       min_q=float(q.min()), max_D2=float(D2.max()),
                       frac_praw_lt05=float((p_raw < 0.05).mean()),
                       frac_pgc_lt05=float((p_gc < 0.05).mean()))
            summary.append(row)
            print(f"  {tag}: rank={rank} kax={kax} lambda={lam:.3f} "
                  f"nCand(q<.05)={hit.sum()} nBlocks={nblk} minQ={q.min():.3g}",
                  flush=True)
        del Y, Yc, Yc_m
    sdf = pd.DataFrame(summary)
    sdf.to_csv(f"{OUTD}/summary.csv", index=False)
    with open(f"{OUTD}/summary.json", "w") as f:
        json.dump(summary, f, indent=2)
    print("\n=== SUMMARY ===")
    print(sdf.to_string(index=False))


if __name__ == "__main__":
    main()
