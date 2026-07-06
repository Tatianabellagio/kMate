#!/usr/bin/env python
"""Reduced-PC structure-corrected multivariate-climate regression (calibration check).

Companion to rda_prda.py. Instead of the parametric Mahalanobis-of-loadings null
(which comes out DEFLATED, lambda~0.8, on this 31-site data), this uses a proper
per-locus partial-F test of the multivariate climate model, then reports/tunes the
canonical Devlin-Roeder genomic-control lambda (chi^2_1 scale) toward ~1.0-1.2.

Per locus (over the 31 sites), F-test of the climate block controlling K structure PCs:
    y = Delta-p               (centered over sites)
    reduced design  = [K site-PCs]         (structure)
    full design     = [K site-PCs, climate q axes]
    partial F = (SS_climate / q) / (SS_resid / (n-1-K-q))
p_raw = F.sf(F, q, n-1-K-q).  Genomic control on the chi^2_1 transform of p_raw:
    lambda_gc = median(chi2_1^{-1}(p_raw)) / qchisq(0.5,1)
    p_gc = chi2_1.sf( chi2_1^{-1}(p_raw) / lambda_gc )
BH-FDR on p_gc. Reports lambda_gc for each K so K can be chosen for calibration
(target lambda ~1.0-1.2, NOT deflated), and the FDR-candidate clq0.9 block count.

Outputs under results/grenenet_gea/gea_newpanel/rda/ (pcreg_*).
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
CLIMSET = os.environ.get("RDA_CLIM", "bio5,bio6,bio12,bio15").split(",")
KSTRUCT = [int(k) for k in os.environ.get("RDA_KSTRUCT", "0,1,2,3,4,5").split(",")]
CLASSES = os.environ.get("RDA_CLASSES", "snp,nonsnp,sv").split(",")
P0_FOR = {"snp": "snp", "nonsnp": "nonsnp", "sv": "nonsnp"}
CHI2_1_MED = stats.chi2.ppf(0.5, 1)


def bh_fdr(p):
    p = np.asarray(p, float); n = p.size
    order = np.argsort(p)
    q = np.minimum.accumulate((p[order] * n / (np.arange(n) + 1))[::-1])[::-1]
    out = np.empty(n); out[order] = np.clip(q, 0, 1); return out


def load_class(cls):
    dims = open(f"{LFMM}/lfmm_{cls}_site_dims.txt").read().split()
    ns, nv = int(dims[0]), int(dims[1]); assert ns == NSITES
    Y = np.fromfile(f"{LFMM}/lfmm_{cls}_site_Y.f64", dtype=np.float64).reshape(ns, nv)
    rec = pd.read_csv(f"{RECD}/{cls}_gen9.records.csv"); assert len(rec) == nv
    p0 = np.load(f"{ROOT}/af_store/p0_{P0_FOR[cls]}.npy").astype(np.float64)
    col = rec["col"].to_numpy()
    site_af = Y + p0[col][None, :]
    maf = np.minimum(site_af.mean(0), 1.0 - site_af.mean(0))
    return Y, rec, maf


def load_climate():
    cols = [pd.read_csv(f"{ENVD}/env_site_{a}.csv").iloc[:, 0].to_numpy(float) for a in CLIMSET]
    X = np.column_stack(cols)
    X = X - X.mean(0, keepdims=True)
    X = X / X.std(0, ddof=0, keepdims=True)
    return X


def partial_F(Yc_m, X, k, Uall, SS_tot):
    """Per-locus partial-F of climate controlling k site-PCs. Returns F, p_raw, dfs.
    Uall = precomputed left singular vectors of Yc_m (structure basis); SS_tot precomputed."""
    n, p = Yc_m.shape
    q = X.shape[1]
    if k > 0:
        Z = Uall[:, :k]
        Xr = X - Z @ (Z.T @ X)
        ZtY = Z.T @ Yc_m
        SS_red = np.einsum("ki,ki->i", ZtY, ZtY)
    else:
        Z = np.zeros((n, 0)); Xr = X; SS_red = np.zeros(p)
    # orthonormal basis of climate after removing structure
    Qx, rq = np.linalg.qr(Xr)
    qeff = int(np.linalg.matrix_rank(Xr))
    B = Qx[:, :qeff].T @ Yc_m
    SS_clim = np.einsum("ki,ki->i", B, B)
    SS_res = np.clip(SS_tot - SS_red - SS_clim, 1e-300, None)
    df2 = n - 1 - k - qeff
    F = (SS_clim / qeff) / (SS_res / df2)
    p_raw = stats.f.sf(F, qeff, df2)
    return F, p_raw, qeff, df2


def main():
    genes = load_genes()
    X = load_climate()
    print(f"climate = {CLIMSET}  X{X.shape}", flush=True)
    summary = []
    for cls in CLASSES:
        t0 = time.time()
        Y, rec, maf = load_class(cls)
        keep = maf >= MAF_MIN
        blocks = assign_clq09_blocks(rec["chrom"].to_numpy(), rec["pos"].to_numpy())
        rec = rec.assign(site_maf=maf, block_clq09=blocks)
        idx = np.where(keep)[0]
        sub_base = rec.iloc[idx].reset_index(drop=True)
        Yc = Y - Y.mean(0, keepdims=True)
        Yc_m = np.ascontiguousarray(Yc[:, idx])
        Uall, _, _ = np.linalg.svd(Yc_m, full_matrices=False)   # structure basis, once
        SS_tot = np.einsum("ij,ij->j", Yc_m, Yc_m)
        print(f"[{cls}] maf>=.05 keeps {keep.sum():,}/{len(maf):,} ({time.time()-t0:.0f}s)", flush=True)
        for k in KSTRUCT:
            F, p_raw, qeff, df2 = partial_F(Yc_m, X, k, Uall, SS_tot)
            chi1 = stats.chi2.isf(np.clip(p_raw, 1e-300, 1.0), 1)
            lam = np.median(chi1) / CHI2_1_MED
            p_gc = stats.chi2.sf(chi1 / max(lam, 1e-9), 1)
            q = bh_fdr(p_gc)
            hit = q < 0.05
            nblk = pd.Series(sub_base["block_clq09"].to_numpy()[hit]).nunique() if hit.any() else 0
            if hit.any():
                sub = sub_base.iloc[np.where(hit)[0]].copy()
                sub["F"] = F[hit]; sub["p_raw"] = p_raw[hit]; sub["p_gc"] = p_gc[hit]; sub["fdr_q"] = q[hit]
                sub = annotate_svs(sub, flank=2000, genes=genes)
                sub.sort_values("fdr_q").to_csv(f"{OUTD}/pcreg_cand_{cls}_k{k}.csv", index=False)
            row = dict(cls=cls, kstruct=k, qeff=qeff, df2=df2, n_loci=int(len(idx)),
                       lambda_gc=float(lam), n_cand=int(hit.sum()), n_blocks_cand=int(nblk),
                       min_q=float(q.min()), frac_praw_lt05=float((p_raw < 0.05).mean()),
                       frac_pgc_lt05=float((p_gc < 0.05).mean()))
            summary.append(row)
            print(f"  {cls} k={k}: qeff={qeff} df2={df2} lambda_gc={lam:.3f} "
                  f"nCand(q<.05)={hit.sum()} nBlocks={nblk} minQ={q.min():.3g} "
                  f"frac_praw<.05={row['frac_praw_lt05']:.3f}", flush=True)
        del Y, Yc, Yc_m
    sdf = pd.DataFrame(summary)
    sdf.to_csv(f"{OUTD}/pcreg_summary.csv", index=False)
    with open(f"{OUTD}/pcreg_summary.json", "w") as f:
        json.dump(summary, f, indent=2)
    print("\n=== PCREG SUMMARY ===")
    print(sdf.to_string(index=False))


if __name__ == "__main__":
    main()
