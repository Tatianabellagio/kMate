#!/usr/bin/env python
"""The SNPs that track the GPX6 haplotype as well as the deletion -- do any of them have a
mechanism?

audit_gpx6.py check 7 found that the 1,164 bp deletion is not the best MARKER at its
locus: 5 of 1,861 window SNPs match or beat its garden-level climate gradient. On
detection grounds the SV adds nothing, which is the standing "no kMate gain" result.

But detection is not the claim. The claim is causal candidacy: a 1.2 kb deletion that
removes open chromatin and destroys 15 TF binding sites is a mechanistic hypothesis, and a
SNP in perfect LD with it is not -- UNLESS that SNP also lands somewhere regulatory. So
this asks, for every SNP on the haplotype: where does it sit, does it overlap an ATAC
peak, and does it change a TF binding site under the same FIMO contrast the deletion went
through (tfbs_turnover.run_one, same pad, same q and low-complexity filters)?

If the rival SNPs are silent, LD-equivalence is irrelevant: the deletion is the only
variant on the haplotype that could do anything, and that is the whole argument for
carrying SVs through a scan that SNPs already detect.

env: kmate + the meme env for FIMO. Compute node, ~2 min.
Writes results/gpx6_rival_snps.csv.
"""
from __future__ import annotations
import os, sys
import numpy as np, pandas as pd, pysam

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(HERE, "..", "..")))
sys.path.insert(0, HERE)
import lib                                                       # noqa: E402
import axis_clusters as ac                                       # noqa: E402
import audit_gpx6 as AU                                          # noqa: E402
import tfbs_turnover as TF                                       # noqa: E402
import atac_overlap as AO                                        # noqa: E402

OUT = f"{HERE}/results"
CHROM, POS, R_SV, WIN, PAD = "Chr4", 7011705, -0.694903, 25_000, 400
R_MIN = -0.60          # SNPs at least this graded; the deletion itself is -0.695


def main():
    PM = pd.read_csv(ac.POOLMETA); site = PM.site.values
    A = ac.site_climate(); A.index = sorted(PM.site.astype(int).unique())
    pc1 = A["pc1"]

    Ms = np.load(f"{lib.GEA}/common/results/pool_matrices/pool_gen9_snp_af.npy", mmap_mode="r")
    p0s = np.load(f"{lib.AF_STORE}/p0_snp.npy")
    i = np.load(f"{lib.AF_STORE}/index_snp.npz", allow_pickle=True)
    ch, pos = i["chrom"].astype(str), i["pos"]
    sel = np.where((ch == CHROM) & (pos >= POS - WIN) & (pos <= POS + WIN))[0]
    sel = sel[np.minimum(p0s[sel], 1 - p0s[sel]) > 0.05]
    r = AU.garden_r(Ms, p0s, sel, site, pc1)
    T = pd.DataFrame({"pos": pos[sel], "r": r}).sort_values("r")
    T = T[T.r <= R_MIN].reset_index(drop=True)
    print(f"{len(T)} SNPs within +-{WIN//1000} kb with garden-level r <= {R_MIN} "
          f"(deletion r = {R_SV:+.3f})\n")

    # -- where do they sit? ----------------------------------------------------------------
    genes = lib.load_genes()
    g = genes[genes.gene == "AT4G11600"].iloc[0]          # minus strand: TSS = end
    peaks = AO.load_peaks()
    ov = AO.peak_overlap(peaks, [CHROM] * len(T), T.pos.values, np.ones(len(T), int))
    T = pd.concat([T, ov[["in_atac", "atac_n_tissues"]]], axis=1)
    T["dist_to_TSS"] = T.pos - g.end
    T["in_GPX6_promoter"] = T.dist_to_TSS.between(1, 1000)        # 1 kb upstream, minus strand
    T["in_GPX6_body"] = T.pos.between(g.start, g.end)
    T["inside_deletion"] = T.pos.between(POS, POS + 1225 - 1)

    # -- do they change a TF binding site? -------------------------------------------------
    fa_ref = pysam.FastaFile(TF.REF_FA)
    rows = []
    for t in T.itertuples():
        rec = pd.Series(dict(symbol=f"snp_{t.pos}", target_gene="AT4G11600", chrom=CHROM,
                             pos=int(t.pos), ref_len=1, alt_len=1, set="rival", ftier=""))
        res = TF.run_one(rec, PAD, 1e-4, TF.FIMO, fa_ref)
        rows.append(dict(pos=int(t.pos),
                         n_lost=res["n_lost"] if res else np.nan,
                         n_gained=res["n_gained"] if res else np.nan,
                         n_weak=res["n_weakened"] if res else np.nan,
                         n_strong=res["n_strengthened"] if res else np.nan,
                         fam_lost=res["families_lost"] if res else "",
                         fam_gained=res["families_gained"] if res else ""))
    T = T.merge(pd.DataFrame(rows), on="pos", how="left")
    T["any_tfbs_change"] = (T[["n_lost", "n_gained"]].fillna(0).sum(1) > 0)

    T.to_csv(f"{OUT}/gpx6_rival_snps.csv", index=False)
    print("\n" + "=" * 96)
    cols = ["pos", "r", "dist_to_TSS", "in_GPX6_promoter", "in_GPX6_body", "in_atac",
            "atac_n_tissues", "n_lost", "n_gained", "any_tfbs_change"]
    print(T[cols].to_string(index=False))

    print("\n" + "=" * 96)
    print(f"rival SNPs on the haplotype                : {len(T)}")
    print(f"  inside the deletion footprint            : {int(T.inside_deletion.sum())}")
    print(f"  in GPX6's 1 kb promoter                  : {int(T.in_GPX6_promoter.sum())}")
    print(f"  inside the GPX6 gene body                : {int(T.in_GPX6_body.sum())}")
    print(f"  overlapping an ATAC peak                 : {int(T.in_atac.sum())}")
    print(f"  changing any TF binding site (lost/gained): {int(T.any_tfbs_change.sum())}")
    print(f"\nthe deletion, for comparison: promoter (330 bp upstream of the TSS), "
          f"272 bp of ATAC\npeak removed across 2 peaks (1 outright), 15 TFBS lost / 6 gained.")


if __name__ == "__main__":
    main()
