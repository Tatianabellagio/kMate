#!/usr/bin/env python
"""Why are promoters already common among the variants the GEA tested? Two checks, written out
for notebooks/gea_hit_landscape.ipynb (section 2b).

A  Is the promoter call right? Every variant the classifier calls promoter
   (screen_sig_blocks.classify, PROM_BP = 1000) is checked against an INDEPENDENT signed
   distance to the nearest TSS (lib.load_genes, a different gene loader): a promoter call
   should sit 0-1 kb upstream of a TSS, and a non-genic variant NOT called promoter should not.

B  What share of each region would variants have by chance? A variant is an interval, not a
   base: a 2 kb deletion placed at random overlaps a gene body (which the classifier ranks
   first) far more often than a single base does, so a per-bp genome share is the wrong
   baseline for SVs. Here each class x type gets random positions across the five
   chromosomes with the SIZE distribution of the tested variants of that class x type, run
   through the SAME classify(). Uniform placement includes sequence where nothing can be
   called (centromeric repeats, gaps), which mainly affects the TE comparison.

Writes results/gea_region_checkA.csv and results/gea_region_baseline.csv.
env: kmate. Compute node.
"""
import os, sys
import numpy as np, pandas as pd, pysam

HERE = os.path.dirname(os.path.abspath(__file__))
GEA = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(GEA, "genes", "dissection")); sys.path.insert(0, GEA)
import screen_sig_blocks as S                                    # noqa: E402
import tfbs_turnover as TF                                       # noqa: E402

OUT = f"{HERE}/results"
K4 = ["chrom", "pos", "ref_len", "alt_len"]
GROUP = {"1_CDS": "gene body", "2_UTR": "gene body", "4_exon_noncoding": "gene body",
         "5_intron": "gene body", "3_promoter": "promoter", "6_TE": "TE",
         "7_proximal_intergenic": "intergenic", "8_gene_desert": "intergenic"}
REGIONS = ["gene body", "promoter", "TE", "intergenic"]
N_RANDOM = 20000


def main():
    U = pd.read_csv(f"{OUT}/gea_hit_landscape.csv.gz")
    U = U.sort_values(["hit", "min_p"], ascending=[False, True]).drop_duplicates(K4)
    U["g"] = U.tier_1kb.map(GROUP)

    # A -- promoter calls against an independent TSS distance
    ng = U[U.g != "gene body"]
    up = (ng.dist_tss < 0) & (ng.dist_tss >= -1000)
    A = pd.DataFrame([
        dict(set="called promoter", n=int((ng.g == "promoter").sum()),
             share_0_1kb_upstream=up[ng.g == "promoter"].mean()),
        dict(set="non-genic, not called promoter", n=int((ng.g != "promoter").sum()),
             share_0_1kb_upstream=up[ng.g != "promoter"].mean())])
    A.to_csv(f"{OUT}/gea_region_checkA.csv", index=False)
    print(A.to_string(index=False))

    # B -- size-matched random placement through the same classifier
    fa = pysam.FastaFile(TF.REF_FA)
    L = {c: fa.get_reference_length(c) for c in fa.references if c.startswith("Chr") and c[3:].isdigit()}
    chroms = np.array(list(L)); w = np.array([L[c] for c in chroms], float); w /= w.sum()
    rng = np.random.default_rng(1); rows = []
    for cls in ("sv", "smallindel"):
        for k in ("deletion", "insertion"):
            d = U[(U.cls == cls) & (U.kind == k)]
            ref = rng.choice(d.ref_len.values, N_RANDOM, replace=True)
            ch = rng.choice(chroms, N_RANDOM, p=w)
            pos = np.array([rng.integers(1, L[c] - r - 1) for c, r in zip(ch, ref)])
            C, _ = S.classify(pd.DataFrame(dict(chrom=ch, pos=pos, ref_len=ref, alt_len=ref)))
            rs = C.tier.map(GROUP).value_counts(normalize=True)
            ts = d.g.value_counts(normalize=True)
            for r in REGIONS:
                rows.append(dict(cls=cls, kind=k, region=r, random_share=rs.get(r, 0.0),
                                 tested_share=ts.get(r, 0.0)))
            print(f"  {cls} {k}: done", flush=True)
    B = pd.DataFrame(rows)
    B["ratio"] = B.tested_share / B.random_share
    body = B[B.region == "gene body"].set_index(["cls", "kind"])
    nr = 1 - B.set_index(["cls", "kind"]).index.map(body.random_share)
    nt = 1 - B.set_index(["cls", "kind"]).index.map(body.tested_share)
    B["ratio_outside_genes"] = np.where(B.region == "gene body", np.nan,
                                        (B.tested_share / nt) / (B.random_share / nr))
    B.to_csv(f"{OUT}/gea_region_baseline.csv", index=False)
    print(B.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
