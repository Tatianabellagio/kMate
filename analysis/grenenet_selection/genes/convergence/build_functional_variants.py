#!/usr/bin/env python
"""SV / small-indel candidates with a plausible functional handle on a gene.

The question this answers
-------------------------
Of the pooled Bonferroni non-SNP candidates, which ones physically land **inside a
gene, or immediately before or after it** -- i.e. where a length-changing variant has a
mechanism, not just a position? An indel or SV in a gene desert may still be real, but
there is nothing to follow up; one that removes 400 bp of a promoter, or shifts a
reading frame, is a testable hypothesis.

Restricted to `size > 0` -- SVs and small indels only. MNPs (ref_len==alt_len>=2) are
dropped: they are non-SNP, but they are substitutions, not length changes, so the
"structural variant with a functional effect" framing does not apply to them.

What this adds over `variants_classified.csv`
---------------------------------------------
1. **A downstream flank category.** `dissection/screen_sig_blocks.py` only has a
   strand-aware 1 kb *upstream* promoter tier, so a variant sitting just past a gene's
   3' end -- terminators, 3' regulatory elements -- was landing in generic
   `7_proximal_intergenic`. "Right after the gene" now has its own tier, symmetric with
   the promoter one, and it is worth 280 variants that were previously invisible.
2. **Frameshift status for coding indels.** A CDS indel whose length is not a multiple
   of 3 shifts the reading frame for the rest of the transcript -- by far the strongest
   functional prior available without expression data. `size % 3 != 0` for indels; SVs
   overlapping CDS get flagged separately since a large deletion is more likely to
   remove whole exons than to shift a frame.
3. **The gene's evidence** from `candidate_genes.csv`, so functional plausibility and
   statistical support can be read on one row.

Priority is by mechanism, not by p-value:
  F1_CDS_frameshift > F2_CDS_inframe > F3_UTR > F4_promoter > F5_downstream > F6_intron

Outputs -> results/
  functional_variants.csv        every SV/indel with a gene handle, ranked
  functional_shortlist.csv       the above, restricted to genes carrying >=1 gated
                                 line of evidence (see build_convergence.py)

env: kmate.  Run after build_convergence.py, on a compute node.
"""
from __future__ import annotations
import os
import sys
import importlib.util as ilu
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(HERE, "..", "..")))
import lib                                                    # noqa: E402

OUT = f"{HERE}/results"
FLANK_BP = 1000          # "right after the gene", symmetric with the 1 kb promoter

# mechanism tiers, strongest functional prior first
PRIORITY = ["F1_CDS_frameshift", "F2_CDS_inframe", "F3_UTR", "F4_promoter",
            "F5_downstream", "F6_intron"]


def downstream_flank(V: pd.DataFrame) -> pd.DataFrame:
    """Strand-aware: is a non-genic variant within FLANK_BP past a gene's 3' end?

    The promoter tier already claims 1 kb upstream of the TSS, so this is the mirror
    case. A variant can be downstream of one gene and upstream of another; the
    promoter assignment wins, because a 5' regulatory hypothesis is the more specific
    one. Only variants the classifier left non-genic are considered here.
    """
    G = lib.load_genes()
    idx = {}
    for ch, d in G.groupby("chrom"):
        idx[ch] = (d.start.to_numpy(np.int64), d.end.to_numpy(np.int64),
                   d.strand.to_numpy(str), d.gene.to_numpy(str))

    genes, dists = [], []
    for r in V.itertuples():
        gene, dist = "", np.nan
        ch = r.chrom
        if ch in idx:
            st, en, sd, gn = idx[ch]
            vs = int(r.pos)
            ve = vs + int(r.ref_len) - 1
            # 3' end is `end` on the + strand and `start` on the - strand
            lo = np.where(sd == "+", en + 1, st - FLANK_BP)
            hi = np.where(sd == "+", en + FLANK_BP, st - 1)
            m = (vs <= hi) & (ve >= lo)
            if m.any():
                k = int(np.where(m)[0][0])
                gene = gn[k]
                dist = int(vs - en[k]) if sd[k] == "+" else int(st[k] - ve)
        genes.append(gene)
        dists.append(dist)
    V = V.copy()
    V["downstream_gene"] = genes
    V["downstream_dist"] = dists
    return V


def mechanism(V: pd.DataFrame) -> pd.DataFrame:
    """Assign the functional tier and the frameshift flag."""
    V = V.copy()
    # frameshift only means anything for an indel inside coding sequence
    is_cds = V.tier.eq("1_CDS")
    V["frameshift"] = is_cds & V.vclass.eq("smallindel") & V["size"].mod(3).ne(0)
    V["cds_sv"] = is_cds & V.vclass.eq("sv")

    f = pd.Series("", index=V.index, dtype=object)
    f[V.tier.eq("5_intron")] = "F6_intron"
    f[V.downstream_gene.ne("") & f.eq("")] = "F5_downstream"
    f[V.tier.eq("3_promoter")] = "F4_promoter"
    f[V.tier.isin(["2_UTR", "4_exon_noncoding"])] = "F3_UTR"
    f[is_cds] = "F2_CDS_inframe"
    f[V.frameshift] = "F1_CDS_frameshift"
    V["ftier"] = f

    # the gene the functional call is about
    V["target_gene"] = np.where(V.ftier.eq("F5_downstream"),
                                V.downstream_gene, V.gene)
    return V[V.ftier.ne("")]


def main():
    V = pd.read_csv(f"{OUT}/variants_classified.csv")
    A = pd.read_csv(f"{OUT}/candidate_genes.csv")

    n0 = len(V)
    V = V[V["size"] > 0].copy()             # SV + small indel only, no MNPs
    print(f"{n0:,} pooled variants -> {len(V):,} length-changing (SV/indel)")

    nong = V.tier.isin(["6_TE", "7_proximal_intergenic", "8_gene_desert"])
    V.loc[:, "downstream_gene"] = ""
    V.loc[:, "downstream_dist"] = np.nan
    flk = downstream_flank(V[nong])
    V.loc[flk.index, "downstream_gene"] = flk.downstream_gene
    V.loc[flk.index, "downstream_dist"] = flk.downstream_dist
    print(f"  recovered {int(flk.downstream_gene.ne('').sum()):,} variants within "
          f"{FLANK_BP} bp past a gene's 3' end (previously untiered)")

    F = mechanism(V)
    print(f"  {len(F):,} have a gene handle "
          f"({len(V) - len(F):,} are TE/intergenic with none)")

    ev = A.set_index("gene")
    for c in ("n_independent", "L_gea", "L_gwas", "L_cross", "gea_n_clusters",
              "gwas_n_gardens", "locus", "locus_n_genes"):
        F[c] = F.target_gene.map(ev[c]) if c in ev.columns else np.nan
    # genes reached only via the new downstream tier are not in candidate_genes.csv
    F["n_independent"] = F.n_independent.fillna(0).astype(int)

    F["size_inferred"] = F.get("size_inferred", pd.Series(False, index=F.index))
    F["size_inferred"] = F.size_inferred.fillna(False)
    F["ftier_rank"] = F.ftier.map({t: i for i, t in enumerate(PRIORITY)})
    F = F.sort_values(["ftier_rank", "n_independent", "size"],
                      ascending=[True, False, False]).reset_index(drop=True)

    cols = ["chrom", "pos", "ref_len", "alt_len", "size", "vclass", "ftier",
            "frameshift", "cds_sv", "target_gene", "region", "tier",
            "downstream_dist", "in_gea", "in_gwas", "gea_nlp", "gea_n_clusters",
            "gea_sig_clusters", "gea_min_lam", "gwas_nlp", "gwas_n_gardens",
            "gwas_gardens", "gwas_mac", "n_independent", "L_gea", "L_gwas",
            "L_cross", "locus", "locus_n_genes", "te_overlap", "size_inferred",
            "gwas_pos_multiallelic"]
    cols = [c for c in cols if c in F.columns]
    F[cols].to_csv(f"{OUT}/functional_variants.csv", index=False)

    S = F[F.n_independent >= 1]
    S[cols].to_csv(f"{OUT}/functional_shortlist.csv", index=False)

    print("\nby mechanism (all pooled SV/indels with a gene handle):")
    print(pd.crosstab(F.ftier, F.vclass).reindex(PRIORITY).fillna(0)
          .astype(int).to_string())
    print(f"\nCDS frameshift indels: {int(F.frameshift.sum())}   "
          f"SVs overlapping CDS: {int(F.cds_sv.sum())}")
    amb = F[F.size_inferred]
    print(f"  !! {len(amb)} rows have an INFERRED size (multiallelic position, allele "
          f"guessed as the largest record) -- vclass/ftier/frameshift unreliable there; "
          f"{int(amb.frameshift.sum())} of them carry a frameshift call")
    print(f"\nrestricted to genes with >=1 gated line of evidence: {len(S)} variants "
          f"in {S.target_gene.nunique()} genes / {S.locus.nunique()} loci")
    print(pd.crosstab(S.ftier, S.vclass).reindex(PRIORITY).fillna(0)
          .astype(int).to_string())
    print(f"\nwrote {OUT}/functional_variants.csv, {OUT}/functional_shortlist.csv")


if __name__ == "__main__":
    main()
