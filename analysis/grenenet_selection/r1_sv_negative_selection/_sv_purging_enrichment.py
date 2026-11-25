#!/usr/bin/env python
"""Are the candidate purged insertions enriched in gene-related sequence?

Candidates = cold-specific insertions with a global decline (sv_purging_candidates.csv). Their
genomic context (CDS > UTR > intron > TE > promoter <=2 kb > intergenic, first match wins) is
compared against three backgrounds:
  (a) the OTHER cold-specific insertions (same lineage group, not called global) -- the cleanest test
  (b) all tested insertions
  (c) the GENOME: fraction of bp in each class under the same priority rules
Fisher exact per class and for pooled "gene-related" = CDS+UTR+intron+promoter. Candidates are also
collapsed to one per chrom:pos (two alleles of one event were counted twice) as a sensitivity check.

Env: basic. Writes results/sv_adaptive/sv_purging_enrichment.csv.
"""
import os, sys
import numpy as np, pandas as pd
from scipy import stats
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import lib

OUT = f"{lib.GEA}/r1_sv_negative_selection/results/sv_adaptive"
CTX = ["CDS", "UTR", "intron", "TE", "promoter", "intergenic"]
GENIC = ["CDS", "UTR", "intron", "promoter"]
CHRLEN = {"Chr1": 30427671, "Chr2": 19698289, "Chr3": 23459830, "Chr4": 18585056, "Chr5": 26975502}


def genome_fraction():
    g = pd.read_csv(lib.TAIR10_GENES_TE, sep="\t", header=None, comment="#",
                    names=["chrom", "src", "feat", "start", "end", "score", "strand", "frame", "attr"])
    g = g[g.chrom.isin(CHRLEN)]
    genes = g[g.feat == "gene"]
    prom = genes.copy()
    prom["start"] = np.where(genes.strand == "+", genes.start - 2000, genes.end + 1)
    prom["end"] = np.where(genes.strand == "+", genes.start - 1, genes.end + 2000)
    layers = [("promoter", prom), ("TE", g[g.feat == "transposable_element"]),
              ("gene", genes), ("exon", g[g.feat == "exon"]),
              ("UTR", g[g.feat.isin(["five_prime_UTR", "three_prime_UTR"])]), ("CDS", g[g.feat == "CDS"])]
    tot = {c: 0 for c in CTX}
    for ch, L in CHRLEN.items():
        m = {}
        for k, d in layers:
            a = np.zeros(L + 2, bool)
            for s, e in zip(d[d.chrom == ch].start.clip(1, L), d[d.chrom == ch].end.clip(1, L)):
                a[s:e + 1] = True
            m[k] = a[1:L + 1]
        cat = np.full(L, 5, np.int8)                       # intergenic
        cat[m["promoter"]] = 4; cat[m["TE"]] = 3
        cat[m["gene"] & ~m["exon"]] = 2; cat[m["UTR"]] = 1; cat[m["CDS"]] = 0
        for i, c in enumerate(CTX):
            tot[c] += int((cat == i).sum())
    G = sum(tot.values())
    return {c: tot[c] / G for c in CTX}


def fisher(k_in, n_in, k_bg, n_bg):
    t = [[k_in, n_in - k_in], [k_bg, n_bg - k_bg]]
    o, p = stats.fisher_exact(t)
    return o, p


def main():
    A = pd.read_csv(f"{OUT}/sv_purging_annotated.csv")
    I = A[A.cls == "insertion"].copy()
    cold = I[I.origin == "cold-specific"]
    cand = cold[cold.is_global]; rest = cold[~cold.is_global]
    cand_u = cand.drop_duplicates(["chrom", "pos"])
    gf = genome_fraction()
    rows = []
    for lab, fg in (("candidates (alleles)", cand), ("candidates (1 per position)", cand_u)):
        for c in CTX + ["gene-related"]:
            sel = (lambda d: d.context.isin(GENIC)) if c == "gene-related" else (lambda d, c=c: d.context == c)
            k, n = int(sel(fg).sum()), len(fg)
            o1, p1 = fisher(k, n, int(sel(rest).sum()), len(rest))
            o2, p2 = fisher(k, n, int(sel(I).sum()), len(I))
            gfrac = sum(gf[x] for x in GENIC) if c == "gene-related" else gf[c]
            pg = stats.binomtest(k, n, gfrac).pvalue
            rows.append(dict(set=lab, context=c, k=k, n=n, frac=k / n,
                             frac_other_cold=sel(rest).mean(), OR_vs_other_cold=o1, p_vs_other_cold=p1,
                             frac_all_ins=sel(I).mean(), OR_vs_all_ins=o2, p_vs_all_ins=p2,
                             frac_genome=gfrac, fold_vs_genome=(k / n) / gfrac, p_vs_genome=pg))
    R = pd.DataFrame(rows); R.to_csv(f"{OUT}/sv_purging_enrichment.csv", index=False)
    pd.set_option("display.width", 250)
    for lab in R.set.unique():
        d = R[R.set == lab]
        print(f"\n=== {lab}: n = {d.n.iloc[0]}")
        print(f"{'context':<14}{'cand':>12}{'other cold':>12}{'OR':>7}{'p':>9}{'all ins':>10}{'OR':>7}{'p':>9}"
              f"{'genome':>9}{'fold':>7}{'p':>10}")
        for _, r in d.iterrows():
            print(f"{r.context:<14}{f'{r.k} ({100*r.frac:.0f}%)':>12}{100*r.frac_other_cold:>11.0f}%{r.OR_vs_other_cold:>7.2f}"
                  f"{r.p_vs_other_cold:>9.2g}{100*r.frac_all_ins:>9.0f}%{r.OR_vs_all_ins:>7.2f}{r.p_vs_all_ins:>9.2g}"
                  f"{100*r.frac_genome:>8.0f}%{r.fold_vs_genome:>7.2f}{r.p_vs_genome:>10.2g}")
    print(f"\n[wrote] {OUT}/sv_purging_enrichment.csv")


if __name__ == "__main__":
    main()
