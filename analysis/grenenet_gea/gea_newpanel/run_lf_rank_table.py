#!/usr/bin/env python
"""Per-axis significant-BLOCK table from the latent-factor-adjusted PARTIAL-RANK
(Spearman) climate-GEA, SNP vs non-SNP vs SV, with cross-class block overlaps.

This is the controlled (structure-adjusted, monotone/rank) analog of the LFMM
sig_blocks table: for every climate axis (bio1..19 + pc1) and every class, run the
K=16 latent-factor-adjusted partial Spearman (Agent-1 method, run_lf_rank.py),
GIF-calibrate the p per (class,axis) [deflate-only], collapse Bonferroni- and
FDR-significant records to clq0.9 haploblocks, and count blocks per class + the
SNP/non-SNP/SV overlaps. The per-axis GIF is reported and axes whose GIF stays
outside [0.85,1.25] even after calibration are flagged 'inflated'.

Efficiency: the latent-factor residual projector Mp = I - D(D'D)^-1 D' (D=[1|U]) is
the SAME for every axis, so we rank+residualize each class's AF ONCE, then every axis
is a single mat-vec. Latent factors U (K=16) are SNP-derived (common structure
reference for all three classes; 'structure is structure', matching the MSR run).

Out: results/grenenet_gea/gea_newpanel/lf_rank_table/
  sig_blocks_by_axis.csv   cls,axis,GIF,inflated,n_tested,n_bonf_blocks,n_fdr_blocks,min_p
  sig_blocks_overlap.csv   axis, per-class block counts + pairwise overlaps (bonf & fdr)
  + prints the axis x {SNP, non-SNP, SV, SNP∩SV, ...} pivot.
Run in kmate env (numpy/pandas/scipy/sklearn).
"""
from __future__ import annotations
import os, sys
import numpy as np
import pandas as pd
from scipy.stats import chi2
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.dirname(HERE))
import lib, blocks_clq09
from run_lf_rank import rank_avg, latent_factors, m_perp, pval_from_r, gif

CM = f"{lib.GEA}/phase1_replication/class_matrices"
OUT = f"{lib.GEA}/gea_newpanel/lf_rank_table"
AXES = [f"bio{i}" for i in range(1, 20)] + ["pc1"]
CLASSES = ["snp", "nonsnp", "sv"]
K = 16
GEN = 9
MAF_MIN = 0.05
CRED = (0.85, 1.25)


def bh(p):
    p = np.asarray(p, float); n = len(p); o = np.argsort(p); q = np.empty(n)
    q[o] = (p[o] * n) / (np.arange(n) + 1); q[o] = np.minimum.accumulate(q[o][::-1])[::-1]
    return np.clip(q, 0, 1)


def calibrate(p, g):
    """GIF deflate-only: divide chi2 by GIF when GIF>=1, else leave raw."""
    if not np.isfinite(g) or g < 1.0:
        return p
    x2 = chi2.isf(np.clip(p, 1e-300, 1.0), 1)
    return chi2.sf(x2 / g, 1)


def main():
    os.makedirs(OUT, exist_ok=True)
    pools = pd.read_csv(f"{CM}/gen{GEN}.pools.csv")
    n = len(pools)
    bios = [f"bio{i}" for i in range(1, 20)]
    B = pools[bios].to_numpy(float)
    Bz = (B - B.mean(0)) / B.std(0, ddof=0)
    # PC1 of the standardized 355x19 bioclim (plot-level composite axis)
    _, _, Vt = np.linalg.svd(Bz - Bz.mean(0), full_matrices=False)
    pc1 = Bz @ Vt[0]
    env = {ax: pools[ax].to_numpy(float) for ax in bios}
    env["pc1"] = pc1

    # latent factors from SNP AF (common structure reference), projector, env residuals
    print("[U] latent factors K=%d from SNP AF ..." % K, flush=True)
    snp_af = np.load(f"{CM}/snp_gen{GEN}_af.npy")
    U = latent_factors(snp_af, K)
    del snp_af
    Mp = m_perp(U, n)
    df = n - K - 2
    re_res, e_norm = {}, {}
    for ax in AXES:
        r = Mp @ rank_avg(env[ax])[0]
        re_res[ax] = r
        e_norm[ax] = float(np.linalg.norm(r))

    by_axis, sig = [], {}     # sig[(axis,thr)] -> {cls: set(blocks)}
    for cls in CLASSES:
        recs = pd.read_csv(f"{CM}/{cls}_gen{GEN}.records.csv")
        af = np.load(f"{CM}/{cls}_gen{GEN}_af.npy")
        M = af.shape[1]
        blocks = blocks_clq09.assign_clq09_blocks(recs.chrom.to_numpy(), recs.pos.to_numpy())
        keep = recs["maf"].to_numpy(float) >= MAF_MIN
        # rank + residualize each variant's AF ONCE; keep norms
        Ra_res = np.empty((M, n), np.float32)
        for a in range(0, M, 100000):
            b = min(a + 100000, M)
            Ra = rank_avg(af[:, a:b].T.astype(np.float64))
            Ra_res[a:b] = (Ra @ Mp).astype(np.float32)
        del af
        ra_norm = np.linalg.norm(Ra_res, axis=1)
        ra_norm[ra_norm == 0] = np.inf
        print(f"[{cls}] M={M:,} kept(MAF>={MAF_MIN})={int(keep.sum()):,}", flush=True)

        for ax in AXES:
            r = (Ra_res @ re_res[ax]) / (ra_norm * e_norm[ax])
            p = pval_from_r(r, df)
            g = gif(p)
            pf = calibrate(p, g)
            pk = pf[keep]; blk = blocks[keep]
            nb = int(np.isfinite(pk).sum())
            bthr = 0.05 / nb
            q = bh(pk)
            sb = set(blk[(pk < bthr) & (blk != "")])
            fb = set(blk[(q < 0.05) & (blk != "")])
            sig.setdefault((ax, "bonf"), {})[cls] = sb
            sig.setdefault((ax, "fdr"), {})[cls] = fb
            by_axis.append(dict(cls=cls, axis=ax, GIF=round(g, 2),
                inflated=not (CRED[0] <= g <= CRED[1]),
                n_tested=nb, n_bonf_blocks=len(sb), n_fdr_blocks=len(fb),
                min_p=float(np.nanmin(pf))))
        del Ra_res

    A = pd.DataFrame(by_axis)
    A.to_csv(f"{OUT}/sig_blocks_by_axis.csv", index=False)

    ov = []
    for (ax, thr), dd in sig.items():
        s, nn, v = dd.get("snp", set()), dd.get("nonsnp", set()), dd.get("sv", set())
        ov.append(dict(axis=ax, thr=thr, snp=len(s), nonsnp=len(nn), sv=len(v),
                       snp_nonsnp=len(s & nn), snp_sv=len(s & v), nonsnp_sv=len(nn & v)))
    O = pd.DataFrame(ov)
    O.to_csv(f"{OUT}/sig_blocks_overlap.csv", index=False)

    gif_snp = A[A.cls == "snp"].set_index("axis").GIF
    infl = A[A.cls == "snp"].set_index("axis").inflated
    pd.set_option("display.width", 200)
    for thr, label in [("bonf", "Bonferroni"), ("fdr", "FDR q<0.05")]:
        piv = O[O.thr == thr].set_index("axis")[["snp", "nonsnp", "sv",
                "snp_nonsnp", "snp_sv", "nonsnp_sv"]]
        piv = piv.rename(columns={"snp": "SNP", "nonsnp": "non-SNP", "sv": "SV",
                "snp_nonsnp": "SNP∩nonSNP", "snp_sv": "SNP∩SV", "nonsnp_sv": "nonSNP∩SV"})
        piv["GIF(snp)"] = gif_snp.round(2)
        piv["inflated"] = infl
        piv = piv.reindex(AXES)
        print(f"\n===== significant clq0.9 BLOCKS by axis ({label}, K=16 LF-adjusted "
              f"partial-Spearman, GIF-calibrated) =====")
        print(piv.sort_values(["SV", "non-SNP", "SNP"], ascending=False).to_string())
    print(f"\nwrote {OUT}/sig_blocks_by_axis.csv + sig_blocks_overlap.csv")


if __name__ == "__main__":
    main()
