#!/usr/bin/env python
"""Partial-RDA VARIANCE PARTITIONING by variant class (SNP / indel / SV / ALL):
how much of each class's among-site allele-frequency variance is explained by
climate (net of population structure), and does adding indels/SVs explain climate
variance the SNPs can't?

Design (locked with user):
  response  = site allele frequencies [31 sites x L], centered (no Hellinger, no
              scaling, no depth weight -- kMate AF is calibrated/equal-precision)
  predictor = climate PCs (PCA of 19 bioclim across 31 sites; broken-stick, <=5)
  covariate = common neutral structure = top SNP-PCA axes (SAME for every class)
  stat      = pure-climate adjusted R^2 (conditioned on structure) + permutation p
              (permute the 31 climate rows >=2000x); confounded = full - pure.
  MAF       = site-level >=0.05 primary + >=0.01 sensitivity (same filter all classes)
Efficiency: constrained SS = trace(Hx * G) with G = Yz Yz' (31x31) precomputed, so
each permutation is a 31x31 elementwise product -- 2000 perms are ~instant.
"""
from __future__ import annotations
import os, sys
import numpy as np, pandas as pd
sys.path.insert(0, "/global/scratch/users/tbellg/kmate/analysis/grenenet_gea")
import lib

CM = f"{lib.GEA}/phase1_replication/results/class_matrices"
OUT = f"{lib.GEA}/gea_newpanel/results/rda_varpart"; os.makedirs(OUT, exist_ok=True)
BIOS = [f"bio{i}" for i in range(1, 20)]
NPERM = 2000

pools = pd.read_csv(f"{CM}/gen9.pools.csv")
sites = np.sort(pools.site.unique()); n = len(sites)
w = pools.total_flowers.to_numpy(float); site_of = pools.site.to_numpy()


def brokenstick_k(evr):
    p = len(evr)
    bs = np.array([np.sum(1.0 / np.arange(k, p + 1)) for k in range(1, p + 1)]) / p
    above = evr > bs
    return int(np.argmax(~above)) if (~above).any() else p


def site_af(cls):
    af = np.load(f"{CM}/{cls}_gen9_af.npy")               # plots x L
    S = np.empty((n, af.shape[1]), dtype=np.float64)
    for i, s in enumerate(sites):
        m = site_of == s; ws = w[m] / w[m].sum()
        S[i] = ws @ af[m].astype(np.float64)
    del af
    return S


def maf_mask(S, thr):
    p = S.mean(0); mf = np.minimum(p, 1 - p)
    return (mf >= thr) & (S.std(0) > 1e-9)


def _prep(Y, Z):
    """Center Y, residualize on Z, return (G=31x31 Gram of residual, tot, Hz)."""
    Yc = Y - Y.mean(0)
    if Z is not None and Z.shape[1]:
        Hz = Z @ np.linalg.pinv(Z.T @ Z) @ Z.T
        Yz = Yc - Hz @ Yc
    else:
        Hz = None; Yz = Yc
    G = Yz @ Yz.T
    return G, float(np.trace(G)), Hz


def _con(X, Hz, G):
    Xz = X - (Hz @ X if Hz is not None else X.mean(0))
    Hx = Xz @ np.linalg.pinv(Xz.T @ Xz) @ Xz.T
    return float(np.sum(Hx * G))


def rda(Y, X, Z, perm=True):
    G, tot, Hz = _prep(Y, Z)
    obs = _con(X, Hz, G); R2 = obs / tot
    p = np.nan
    if perm:
        rng = np.random.RandomState(0); cnt = 1
        for _ in range(NPERM):
            if _con(X[rng.permutation(n)], Hz, G) >= obs:
                cnt += 1
        p = cnt / (NPERM + 1)
    return R2, p


def radj(R2, p, q):
    d = n - q - p - 1
    return 1 - (1 - R2) * (n - q - 1) / d if d > 0 else np.nan


def main():
    # climate PCs
    sp = pools.drop_duplicates("site").set_index("site").loc[sites]
    bcol = [b for b in BIOS if b in sp.columns]
    Bz = sp[bcol].to_numpy(float); Bz = (Bz - Bz.mean(0)) / Bz.std(0)
    Ub, Sb, _ = np.linalg.svd(Bz - Bz.mean(0), full_matrices=False)
    evr = Sb ** 2 / np.sum(Sb ** 2)
    kc = min(5, max(2, brokenstick_k(evr)))
    CLIM = Ub[:, :kc] * Sb[:kc]
    print(f"climate PCs kept: {kc} (broken-stick), cum var {evr[:kc].sum():.2f}", flush=True)

    # per-class site AF
    AF = {c: site_af(c) for c in ["snp", "smallindel", "sv"]}
    AF["all"] = np.hstack([AF["snp"], AF["smallindel"], AF["sv"]])
    for c in AF:
        print(f"  {c}: site-AF [{n} x {AF[c].shape[1]:,}]", flush=True)

    rows = []
    for thr in [0.05, 0.01]:
        # common neutral covariate = SNP-PCA axes (structure), same for all classes
        Ssnp = AF["snp"][:, maf_mask(AF["snp"], thr)]; Sc = Ssnp - Ssnp.mean(0)
        Us, ss, _ = np.linalg.svd(Sc, full_matrices=False)
        evs = ss ** 2 / np.sum(ss ** 2)
        kz = min(6, max(3, brokenstick_k(evs)))
        Zstruct = Us[:, :kz] * ss[:kz]
        print(f"\nMAF>={thr}: SNP-structure axes kz={kz} (cum var {evs[:kz].sum():.2f})", flush=True)
        for c in ["snp", "smallindel", "sv", "all"]:
            S = AF[c][:, maf_mask(AF[c], thr)]
            R2_pure, p_pure = rda(S, CLIM, Zstruct, perm=True)      # climate | structure
            R2_full, _ = rda(S, CLIM, None, perm=False)             # climate alone
            rows.append(dict(maf=thr, cls=c, L=int(S.shape[1]), kclim=kc, kstruct=kz,
                pure_climate_R2=round(R2_pure, 4), pure_climate_R2adj=round(radj(R2_pure, kc, kz), 4),
                perm_p=round(p_pure, 4), full_climate_R2=round(R2_full, 4),
                confounded=round(R2_full - R2_pure, 4)))
            print(f"  {c:11s} L={S.shape[1]:>9,}  pure-climate R2={R2_pure:.3f} "
                  f"(adj {radj(R2_pure,kc,kz):.3f}, p={p_pure:.4f})  full={R2_full:.3f} "
                  f"confounded={R2_full-R2_pure:.3f}", flush=True)
    df = pd.DataFrame(rows)
    df.to_csv(f"{OUT}/varpart_by_class.csv", index=False)
    print(f"\nwrote {OUT}/varpart_by_class.csv")
    # headline: does ALL beat SNP (added information)?
    print("\n=== added information: pure-climate adj R2, SNP vs +indel/SV (ALL) ===")
    for thr in [0.05, 0.01]:
        d = df[df.maf == thr].set_index("cls")
        print(f"  MAF>={thr}: SNP={d.loc['snp','pure_climate_R2adj']}  "
              f"indel={d.loc['smallindel','pure_climate_R2adj']}  "
              f"SV={d.loc['sv','pure_climate_R2adj']}  ALL={d.loc['all','pure_climate_R2adj']}  "
              f"| ALL-SNP={d.loc['all','pure_climate_R2adj']-d.loc['snp','pure_climate_R2adj']:+.4f}")


if __name__ == "__main__":
    main()
