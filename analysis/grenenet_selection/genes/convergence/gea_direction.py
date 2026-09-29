#!/usr/bin/env python
"""Direction of the GEA associations: do promoter deletions lean warm or cold?

The GEA's saved output carries p-values only. lfmm_signed/ re-ran the identical model (same
input, ridge K = 16, raw p) per axis and kept the effect size B and z-score; every axis
reproduced the saved p-values exactly (max |dlog10 p| <= 1.8e-15), so this is the GEA's own
sign, recovered, not a new statistic.

Sign. Y = delta p of the ALT allele, the axis is z-scored, so z > 0 means the ALT allele rose
more where the axis is higher. For a deletion ALT is the deletion. pc1 points to +bio1 and pc3 to
+bio10 (build_class_matrices.py), so on every temperature axis z > 0 = favoured WARM.

The comparison that matters is promoter deletions against INTERGENIC deletions -- the burden
model's contrast (notebooks/sv_climate_maf.ipynb). Both are deletions relative to Col-0, so any
reference-polarisation effect (non-Col-0 alleles rising in warm gardens because warm-favoured
lineages are more divergent from Col-0) is shared and cancels. Read two ways:

  hits    among the axis's GEA hits: share with z > 0, promoter vs intergenic deletions
  all     over every tested deletion: mean z, promoter minus intergenic -- the whole signed
          distribution, which is where a small shared shift (what a burden test sees) would
          show; stratified by size x MAF bin, 95% CI from resampling LD blocks

The sign comes from the pooled non-SNP scan; a record called only in the SV or small-indel scan
still gets its direction from the pooled one. Writes results/gea_direction.csv.
env: kmate. Compute node.
"""
import os, sys
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import gea_region_by_axis as GR                                  # noqa: E402

OUT = f"{HERE}/results"
SIGNED = f"{HERE}/lfmm_signed"
K4 = ["chrom", "pos", "ref_len", "alt_len"]
TEMP = ["bio1", "bio5", "bio6", "bio8", "bio9", "bio10", "bio11", "pc1", "pc3"]
REPS = 1000


def contrast_mean_z(d, a, b, rng):
    """mean z (a) - mean z (b), averaged over size x MAF strata; block bootstrap CI."""
    d = d[d.grp.isin([a, b])]
    st, si = np.unique(d.stratum.values, return_inverse=True)
    bl, bi = np.unique(d.block.values, return_inverse=True)
    g = (d.grp.values == a).astype(int)                       # 1 = a, 0 = b
    S = np.zeros((len(bl), len(st), 2)); N = np.zeros((len(bl), len(st), 2))
    np.add.at(S, (bi, si, g), d.z.values); np.add.at(N, (bi, si, g), 1)

    def est(Sx, Nx):
        ok = (Nx[..., 0] > 0) & (Nx[..., 1] > 0)
        diff = np.where(ok, Sx[..., 1] / np.where(ok, Nx[..., 1], 1) - Sx[..., 0] / np.where(ok, Nx[..., 0], 1), 0)
        w = np.where(ok, np.minimum(Nx[..., 0], Nx[..., 1]), 0)
        return (diff * w).sum(-1) / np.maximum(w.sum(-1), 1e-12)

    e = est(S.sum(0), N.sum(0))
    boot = []
    for _ in range(REPS // 100):
        W = np.stack([np.bincount(rng.integers(0, len(bl), len(bl)), minlength=len(bl)) for _ in range(100)]).astype(float)
        boot.append(est(np.einsum("rb,bsk->rsk", W, S), np.einsum("rb,bsk->rsk", W, N)))
    boot = np.concatenate(boot)
    return e, np.percentile(boot, 2.5), np.percentile(boot, 97.5), min(1.0, 2 * min((boot <= 0).mean(), (boot >= 0).mean()))


def main():
    U = pd.read_csv(f"{OUT}/gea_hit_landscape.csv.gz")
    U = U.sort_values(["hit", "min_p"], ascending=[False, True]).drop_duplicates(K4)
    U = U[U.kind.isin(["deletion", "insertion"])].copy()
    grp = {"3_promoter": "promoter", "7_proximal_intergenic": "intergenic", "8_gene_desert": "intergenic"}
    U["grp"] = U.tier_1kb.map(grp).fillna("other")
    U["stratum"] = (pd.cut(U["size"], GR.SIZE_BINS, labels=False).astype(str) + "|"
                    + pd.cut(U.MAF, GR.MAF_BINS, labels=False, include_lowest=True).astype(str))
    U["mafr"] = U.MAF.round(6)
    rng = np.random.default_rng(1); rows = []
    for ax in GR.AXES:
        Sg = pd.read_csv(f"{SIGNED}/lfmm_nonsnp_gen9_{ax}_signed.csv.gz")
        Sg["mafr"] = Sg.MAF.round(6)
        D = U.merge(Sg[K4 + ["mafr", "z"]].drop_duplicates(K4 + ["mafr"]), on=K4 + ["mafr"], how="left")
        D = D[D.z.notna()]
        D["axis_hit"] = D.hit_axes.fillna("").str.split(",").map(lambda s: ax in s)
        for cls in ("sv", "smallindel", "all"):
            for kind in ("deletion", "insertion"):
                d = D[(D.kind == kind) & ((D.cls == cls) if cls != "all" else True)]
                h = d[d.axis_hit]
                r = dict(axis=ax, cls=cls, kind=kind, tested=len(d), hits=len(h),
                         warm_share_all_promoter=(d[d.grp == "promoter"].z > 0).mean(),
                         warm_share_all_intergenic=(d[d.grp == "intergenic"].z > 0).mean(),
                         hits_promoter=int((h.grp == "promoter").sum()),
                         hits_intergenic=int((h.grp == "intergenic").sum()),
                         warm_share_hits_promoter=(h[h.grp == "promoter"].z > 0).mean() if (h.grp == "promoter").any() else np.nan,
                         warm_share_hits_intergenic=(h[h.grp == "intergenic"].z > 0).mean() if (h.grp == "intergenic").any() else np.nan)
                e, lo, hi, p = contrast_mean_z(d, "promoter", "intergenic", rng)
                r.update(dz_promoter_minus_intergenic=e, dz_lo=lo, dz_hi=hi, dz_p=p)
                rows.append(r)
        print(f"  {ax}", flush=True)
    R = pd.DataFrame(rows)
    R["temperature_axis"] = R.axis.isin(TEMP)
    R.to_csv(f"{OUT}/gea_direction.csv", index=False)
    pd.set_option("display.width", 250)
    c = ["axis", "tested", "hits", "hits_promoter", "warm_share_hits_promoter", "hits_intergenic",
         "warm_share_hits_intergenic", "dz_promoter_minus_intergenic", "dz_lo", "dz_hi", "dz_p"]
    for cls in ("sv", "all"):
        print(f"\n== {cls} DELETIONS, temperature axes (z > 0 = the deletion rises in warmer gardens)")
        print(R[(R.cls == cls) & (R.kind == "deletion") & R.temperature_axis][c].round(3).to_string(index=False))
    print("\n== all deletions, NON-temperature axes")
    print(R[(R.cls == "all") & (R.kind == "deletion") & ~R.temperature_axis][c].round(3).to_string(index=False))


if __name__ == "__main__":
    main()
