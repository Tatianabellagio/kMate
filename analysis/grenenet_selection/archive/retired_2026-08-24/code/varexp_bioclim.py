#!/usr/bin/env python
"""Design A -- Zhou et al. 2022 (Nature) analog on GrENE-net FOUNDERS.

Zhou showed a graph pangenome (SNP+indel+SV, incl. long-read-assembly SVs) captures more trait
heritability (0.41 vs 0.33) and higher genomic-prediction accuracy than a single linear reference
(SNP-dominated). Their design needs per-INDIVIDUAL genotypes + per-individual phenotypes + a
variant-discovery gap between the two genotype sources.

We replicate that at the founder level. Individuals = the 231 A. thaliana founders. Phenotype =
per-ecotype ORIGIN CLIMATE (bio1-19, worldclim), measured OUTSIDE the pool experiment -- so this
is NOT the pooled selection coefficient `s` (that was the null in varexp_selection.py); the genome
predicts provenance via local adaptation. Genotype "two ways" = the short-read->long-read contrast
INSIDE the arch3 panel:
    * SNP arm      = K_snp        (the layer short-read/linear calling recovers)
    * pangenome arm= K_all / +K_sv (adds the indel+SV layers that need the long-read graph)

Question (per bio-trait): does adding the SV layer raise heritability and out-of-sample
predictability of climate adaptation beyond SNPs, in a ~97%-selfing species?

Reuses the AI-REML / GBLUP-CV estimators from varexp_selection.py verbatim (same h2, joint-2GRM
LRT, repeated k-fold CV). Reads class_grms.npz + untagged_grms.npz + the bioclim CSV.
Writes bioclim_varexp.csv + bioclim_varexp_summary.json. Env: kmate. Light (231x231 from cache).
"""
from __future__ import annotations
import os, sys, json
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import lib
import varexp_selection as vs          # aireml, h2_1grm, joint_2grm, gblup_cv, rint

OUT = f"{lib.GEA}/r3_persite_gwas/results/varexp"
BIOCLIM_CSV = ("/global/scratch/users/tbellg/gea_grene-net/key_files/"
               "1001g_regmap_grenet_ecotype_info_corrected_bioclim_2024May16.csv")
BIOS = [f"bio{i}" for i in range(1, 20)]
# marginal single-GRM classes to report
MARGINAL = ["snp", "nonsnp", "indel", "sv", "all", "snp_matched", "untagged_r02"]


def zscore(x):
    x = np.asarray(x, float)
    return (x - x.mean()) / x.std()


def main():
    G = np.load(f"{OUT}/class_grms.npz", allow_pickle=True)
    U = np.load(f"{OUT}/untagged_grms.npz", allow_pickle=True)
    nmark = json.loads(str(G["n_markers"]))
    gf = G["founders"].astype(str)
    assert list(gf) == list(U["founders"].astype(str)), "GRM founder order mismatch"

    # ---- phenotype: per-ecotype origin bioclim, joined to GRM founders by ecotypeid ----
    bio = pd.read_csv(BIOCLIM_CSV).set_index("ecotypeid")
    bio.index = bio.index.astype(str)
    have = [i for i, f in enumerate(gf) if f in bio.index]
    matched = gf[have]
    Y = bio.loc[matched, BIOS].astype(float)
    ok = Y.notna().all(axis=1).values            # founders with complete bioclim
    have = list(np.array(have)[ok]); matched = matched[ok]; Y = Y.loc[matched]
    n = len(have)
    print(f"founders genotyped & with bioclim: {n} / {len(gf)}", flush=True)
    missing = [f for f in gf if f not in bio.index]
    if missing:
        print(f"  {len(missing)} founders lack bioclim (dropped): {missing[:12]}"
              + (" ..." if len(missing) > 12 else ""), flush=True)

    # subset all GRMs to the matched founders, same order
    Ks = {c: G[f"K_{c}"][np.ix_(have, have)] for c in
          ["snp", "nonsnp", "indel", "sv", "all", "snp_matched"]}
    Ks["untagged_r02"] = U["K_untagged_r02"][np.ix_(have, have)]
    tri = np.triu_indices(n, 1)
    print("corr(K_snp,K_sv)=%.3f  corr(K_snp,K_nonsnp)=%.3f  corr(K_snp,K_untagged_r02)=%.3f" % (
        np.corrcoef(Ks["snp"][tri], Ks["sv"][tri])[0, 1],
        np.corrcoef(Ks["snp"][tri], Ks["nonsnp"][tri])[0, 1],
        np.corrcoef(Ks["snp"][tri], Ks["untagged_r02"][tri])[0, 1]), flush=True)

    rows = []
    for b in BIOS:
        for transform in ("raw", "rint"):
            y0 = Y[b].values
            y = zscore(y0) if transform == "raw" else vs.rint(y0)
            y = np.asarray(y, float)

            # (1) marginal single-GRM h2 (+SE, +LRT p vs 0)
            for c in MARGINAL:
                h2, se, p = vs.h2_1grm(y, Ks[c])
                rows.append(dict(trait=b, transform=transform, method="marginal",
                                 cls=c, value=h2, se=se, p=p))

            # (2) joint 2-GRM partitions + LRT for the ADDED layer beyond SNP
            for k2, tag in [("sv", "sv"), ("nonsnp", "nonsnp"), ("untagged_r02", "untagged")]:
                frac, p_add = vs.joint_2grm(y, Ks["snp"], Ks[k2])
                for cl, v in zip([f"snp|{tag}", tag, f"resid|{tag}"], frac):
                    rows.append(dict(trait=b, transform=transform, method="joint2",
                                     cls=cl, value=v, se=np.nan,
                                     p=(p_add if cl == tag else np.nan)))

            # (3) GBLUP repeated k-fold CV -- predictive R2 & the gain from the added layer
            #     (Zhou composite form: extra variant class as a separate random effect)
            r_snp,  r2_snp  = vs.gblup_cv(y, [Ks["snp"]])
            r_ssv,  r2_ssv  = vs.gblup_cv(y, [Ks["snp"], Ks["sv"]])       # SNP + SV layer
            r_snn,  r2_snn  = vs.gblup_cv(y, [Ks["snp"], Ks["nonsnp"]])   # SNP + non-SNP (indel+SV)
            r_all,  r2_all  = vs.gblup_cv(y, [Ks["all"]])                 # full pangenome (single GRM)
            r_sv,   r2_sv   = vs.gblup_cv(y, [Ks["sv"]])                  # SV marginal
            r_ns,   r2_ns   = vs.gblup_cv(y, [Ks["nonsnp"]])              # non-SNP marginal
            r_m,    r2_m    = vs.gblup_cv(y, [Ks["snp_matched"]])         # marker-count control
            for cl, rr, r2 in [("snp", r_snp, r2_snp), ("snp+sv", r_ssv, r2_ssv),
                               ("snp+nonsnp", r_snn, r2_snn), ("all", r_all, r2_all),
                               ("sv", r_sv, r2_sv), ("nonsnp", r_ns, r2_ns),
                               ("snp_matched", r_m, r2_m)]:
                rows.append(dict(trait=b, transform=transform, method="predcv",
                                 cls=cl, value=r2, se=np.nan, p=rr))     # p col = pearson r
            rows.append(dict(trait=b, transform=transform, method="predcv",
                             cls="gain_sv", value=r2_ssv - r2_snp, se=np.nan, p=r_ssv - r_snp))
            rows.append(dict(trait=b, transform=transform, method="predcv",
                             cls="gain_nonsnp", value=r2_snn - r2_snp, se=np.nan, p=r_snn - r_snp))
            print(f"  {b:6s} {transform:4s}  r2_snp={r2_snp:+.3f}  "
                  f"gain_sv={r2_ssv-r2_snp:+.3f}  gain_nonsnp={r2_snn-r2_snp:+.3f}", flush=True)

    df = pd.DataFrame(rows)
    df.to_csv(f"{OUT}/bioclim_varexp.csv", index=False)

    # ---- summary: aggregate across the 19 bio-traits (raw transform) ----
    d = df[df.transform == "raw"]
    def marg(c):
        v = d[(d.method == "marginal") & (d.cls == c)].value
        return float(np.median(v)), float(np.mean(v))
    pcv = d[d.method == "predcv"]
    def cv(c):
        v = pcv[pcv.cls == c].value
        return float(np.median(v)), float(np.mean(v))
    def gainstat(cls):
        g = pcv[pcv.cls == cls].value.values
        return {"median": float(np.median(g)), "mean": float(np.mean(g)),
                "min": float(np.min(g)), "max": float(np.max(g)),
                "n_traits_positive": int((g > 0).sum()), "n_traits": int(len(g))}
    def lrt(tag):
        p = d[(d.method == "joint2") & (d.cls == tag)].p.values
        return {"median": float(np.median(p)), "n_sig_0.05": int((p < 0.05).sum()),
                "n_traits": int(len(p))}

    summary = {
        "n_founders": n, "n_bioclim_traits": len(BIOS), "n_markers": nmark,
        "corr_Ksnp_Ksv": float(np.corrcoef(Ks["snp"][tri], Ks["sv"][tri])[0, 1]),
        "corr_Ksnp_Knonsnp": float(np.corrcoef(Ks["snp"][tri], Ks["nonsnp"][tri])[0, 1]),
        "h2_marginal_median_mean": {c: marg(c) for c in MARGINAL},
        "predcv_r2_median_mean": {c: cv(c) for c in
                                  ["snp", "snp+sv", "snp+nonsnp", "all", "sv", "nonsnp", "snp_matched"]},
        "gain_sv_r2": gainstat("gain_sv"),
        "gain_nonsnp_r2": gainstat("gain_nonsnp"),
        "joint2_sv_LRT_p": lrt("sv"),
        "joint2_nonsnp_LRT_p": lrt("nonsnp"),
        "joint2_untagged_LRT_p": lrt("untagged"),
    }
    with open(f"{OUT}/bioclim_varexp_summary.json", "w") as fh:
        json.dump(summary, fh, indent=2)
    gsv, gns = pcv[pcv.cls == "gain_sv"].value.values, pcv[pcv.cls == "gain_nonsnp"].value.values
    psv = d[(d.method == "joint2") & (d.cls == "sv")].p.values
    pns = d[(d.method == "joint2") & (d.cls == "nonsnp")].p.values
    print("\n=== SUMMARY (raw, median over 19 bioclim traits) ===")
    print(f"n_founders={n}")
    print("marginal h2:  " + "  ".join(f"{c}={marg(c)[0]:.2f}" for c in MARGINAL))
    print("predcv   R2:  " + "  ".join(f"{c}={cv(c)[0]:+.3f}" for c in
                                       ["snp", "snp+sv", "snp+nonsnp", "all", "sv", "nonsnp"]))
    print(f"gain_sv     R2 (snp+sv - snp):     median={np.median(gsv):+.4f} "
          f"range=[{gsv.min():+.3f},{gsv.max():+.3f}] pos={int((gsv>0).sum())}/{len(gsv)}")
    print(f"gain_nonsnp R2 (snp+nonsnp - snp): median={np.median(gns):+.4f} "
          f"range=[{gns.min():+.3f},{gns.max():+.3f}] pos={int((gns>0).sum())}/{len(gns)}")
    print(f"joint2 LRT beyond SNP: SV median p={np.median(psv):.3f} sig={int((psv<0.05).sum())}/{len(psv)}"
          f"  |  nonSNP median p={np.median(pns):.3f} sig={int((pns<0.05).sum())}/{len(pns)}")
    print(f"\nwrote {OUT}/bioclim_varexp.csv + bioclim_varexp_summary.json")


if __name__ == "__main__":
    main()
