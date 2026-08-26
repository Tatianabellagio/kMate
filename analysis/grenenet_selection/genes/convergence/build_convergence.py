#!/usr/bin/env python
"""Cross the GEA and GWAS non-SNP candidate pools and count lines of evidence.

The premise
-----------
Neither upstream scan is calibrated well enough to trust on p-value alone:

  * the GEA is raw uncalibrated LFMM (lambda ~1.7 median, up to 3.13 on pc3) --
    memory `raw-lfmm-over-wza-decision`;
  * the per-garden GWAS has a tail inflated 5-60x, worst exactly where its hits live --
    memory `persite-gwas-low-mac-tail-inflation`.

Both inflations are *structure* artifacts, but they act on different traits (climate at
the site vs. founder selection coefficient within a garden), different units (352 pools
vs. 231 founders), and different software (LFMM vs. GEMMA). So a locus that shows up in
both is not explained by either scan's failure mode, and recurrence WITHIN a scan --
across uncorrelated climate clusters, or across independent gardens -- is likewise not
what a single inflated axis or a single badly-behaved garden produces.

Hence: rank by independent lines of evidence, not by p.

  E1  gea_multi_cluster   significant on >=2 of the 7 empirical climate clusters
  E2  gea_clean_axis      significant on >=1 axis with lambda < 2.0
  E3  gwas_multi_garden   Bonferroni in >=2 of the 30 gardens
  E4  both_scans          the gene is hit by BOTH the GEA and the GWAS

Every variant is assigned to a gene AT ITS OWN POSITION via the TAIR10 GFF (the CARK
lesson: a block's attributed gene can be one the lead variant never touches). The
classification code is imported from `dissection/screen_sig_blocks.py` rather than
reimplemented.

Outputs -> results/
  variants_classified.csv  every pooled variant with its genomic context
  candidate_genes.csv      one row per gene, all evidence columns, both scans
  convergent_genes.csv     genes carrying >=2 independent lines of evidence

env: kmate.  Run on a compute node (parses the full TAIR10 GFF).
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
GENIC = ["1_CDS", "2_UTR", "3_promoter", "4_exon_noncoding", "5_intron"]
LOCUS_MERGE_BP = 200_000        # genes closer than this are one locus, not two


def _screen_mod():
    """Import dissection/screen_sig_blocks.py for its GFF classifier."""
    p = os.path.join(os.path.dirname(HERE), "dissection", "screen_sig_blocks.py")
    spec = ilu.spec_from_file_location("screen_sig_blocks", p)
    mod = ilu.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def norm_chrom(s: pd.Series) -> pd.Series:
    """'chr4' / 'Chr4' -> 'Chr4' (the GFF and lib.CHROMS convention)."""
    return s.str.replace(r"^chr", "Chr", case=False, regex=True)


def load_pools() -> pd.DataFrame:
    G = pd.read_csv(f"{OUT}/gea_pool.csv")
    W = pd.read_csv(f"{OUT}/gwas_pool.csv")
    G["chrom"] = norm_chrom(G.chrom)
    W["chrom"] = norm_chrom(W.chrom)

    G = G.rename(columns={"best_nlp": "gea_nlp", "best_axis": "gea_best_axis",
                          "min_lam": "gea_min_lam", "n_axes": "gea_n_axes",
                          "n_clusters": "gea_n_clusters",
                          "sig_axes": "gea_sig_axes",
                          "sig_clusters": "gea_sig_clusters",
                          "cleanest_axis": "gea_cleanest_axis",
                          "has_clean_axis": "gea_clean_axis",
                          "sig_classes": "gea_classes", "MAF": "gea_MAF"})
    W["size_inferred"] = W.get("size_inferred", False)
    W = W.rename(columns={"best_nlp": "gwas_nlp", "best_garden": "gwas_best_garden",
                          "n_gardens": "gwas_n_gardens", "gardens": "gwas_gardens",
                          "mac": "gwas_mac", "maf": "gwas_maf",
                          "bio1_min": "gwas_bio1_min", "bio1_max": "gwas_bio1_max",
                          "sig_classes": "gwas_classes",
                          "pos_multiallelic": "gwas_pos_multiallelic"})

    gcols = ["chrom", "pos", "ref_len", "alt_len", "size", "vclass", "gea_MAF",
             "gea_best_axis", "gea_nlp", "gea_min_lam", "gea_n_axes",
             "gea_n_clusters", "gea_sig_axes", "gea_sig_clusters",
             "gea_cleanest_axis", "gea_clean_axis", "gea_classes", "block"]
    wcols = ["chrom", "pos", "ref_len", "alt_len", "size", "vclass", "gwas_nlp",
             "gwas_best_garden", "gwas_n_gardens", "gwas_gardens", "gwas_mac",
             "gwas_maf", "gwas_bio1_min", "gwas_bio1_max", "gwas_classes",
             "gwas_pos_multiallelic", "size_inferred"]

    # union on position (the GWAS arrays carry no allele, so position is the only
    # key available for the cross-scan merge -- see panel-multiallelic-pos-key-trap;
    # gwas_pos_multiallelic marks where that key is ambiguous)
    U = pd.merge(G[gcols], W[wcols], on=["chrom", "pos"], how="outer",
                 suffixes=("", "_w"))
    for c in ("ref_len", "alt_len", "size", "vclass"):
        U[c] = U[c].where(U[c].notna(), U[f"{c}_w"])
        U = U.drop(columns=[f"{c}_w"])
    U["in_gea"] = U.gea_nlp.notna()
    U["in_gwas"] = U.gwas_nlp.notna()
    U["ref_len"] = U.ref_len.fillna(1).astype(int)
    U["alt_len"] = U.alt_len.fillna(1).astype(int)
    return U


def add_loci(A: pd.DataFrame, merge_bp: int = LOCUS_MERGE_BP) -> pd.DataFrame:
    """Group candidate genes into physical loci.

    A gene list is NOT a locus list here. Neighbouring genes in one LD block are hit
    by the same haplotype and must not be counted as independent evidence -- the Chr4
    ~6.98-7.21 Mb region alone contributes 8 convergent genes, and it is the cluster
    already flagged as "possibly one adaptive haplotype" in `genes-expl-candidate-
    pipeline`. This is the LD-scale version of the block-merge artifact in
    `sv-block-merge-artifact`. Genes within `merge_bp` are collapsed to one locus id,
    and `locus_n_genes` says how many genes a locus contributes.
    """
    G = lib.load_genes().set_index("gene")
    A = A.copy()
    A["gene_start"] = A.gene.map(G.start)
    A["gene_end"] = A.gene.map(G.end)
    A = A.sort_values(["chrom", "gene_start"])
    ids, cur, prev = [], -1, None
    for r in A.itertuples():
        if prev is None or r.chrom != prev[0] or (r.gene_start - prev[1]) > merge_bp:
            cur += 1
        ids.append(f"L{cur:04d}")
        prev = (r.chrom, r.gene_end)
    A["locus"] = ids
    A["locus_n_genes"] = A.groupby("locus")["gene"].transform("size")
    return A


def annotate(genes: list[str]) -> pd.DataFrame:
    """TAIR GO + UniProt, via the annotator mandated in CLAUDE.md."""
    p = os.path.join(os.path.dirname(HERE), "..", "r2_gea_nonsnp",
                     "phase1_replication", "annotate_genes_tair_uniprot.py")
    spec = ilu.spec_from_file_location("ann_tu", os.path.abspath(p))
    mod = ilu.module_from_spec(spec)
    spec.loader.exec_module(mod)
    print(f"annotating {len(genes)} genes (TAIR GO + UniProt) ...", flush=True)
    return mod.annotate(sorted(genes))


def main():
    os.makedirs(OUT, exist_ok=True)
    U = load_pools()
    print(f"pooled variants: {len(U):,}  "
          f"(GEA {int(U.in_gea.sum()):,}, GWAS {int(U.in_gwas.sum()):,}, "
          f"both {int((U.in_gea & U.in_gwas).sum()):,})")

    sm = _screen_mod()
    print("classifying against TAIR10 GFF ...", flush=True)
    V, _G = sm.classify(U)
    V.to_csv(f"{OUT}/variants_classified.csv", index=False)
    print(V.tier.value_counts().sort_index().to_string())

    # ---- gene level: the variant's OWN gene (genic/promoter only) -------------
    g = V[V.tier.isin(GENIC) & V.gene.ne("")].copy()
    print(f"\n{len(g):,} variants land in a gene or its promoter "
          f"-> {g.gene.nunique():,} genes")

    def agg(d: pd.DataFrame) -> pd.Series:
        gea, gw = d[d.in_gea], d[d.in_gwas]
        return pd.Series({
            "chrom": d.chrom.iat[0],
            "n_variants": len(d),
            "regions": ",".join(sorted(set(d.region))),
            "vclasses": ",".join(sorted(set(d.vclass.dropna()))),
            "max_size": int(d["size"].max()),
            "in_gea": bool(len(gea)),
            "in_gwas": bool(len(gw)),
            "gea_nlp": gea.gea_nlp.max() if len(gea) else np.nan,
            "gea_n_clusters": int(gea.gea_n_clusters.max()) if len(gea) else 0,
            "gea_n_axes": int(gea.gea_n_axes.max()) if len(gea) else 0,
            "gea_min_lam": gea.gea_min_lam.min() if len(gea) else np.nan,
            "gea_sig_clusters": (",".join(sorted({c for s in gea.gea_sig_clusters
                                                  for c in str(s).split(",")}))
                                 if len(gea) else ""),
            "gwas_nlp": gw.gwas_nlp.max() if len(gw) else np.nan,
            "gwas_n_gardens": int(gw.gwas_n_gardens.max()) if len(gw) else 0,
            "gwas_gardens": (",".join(sorted({x for s in gw.gwas_gardens
                                              for x in str(s).split(",")}))
                             if len(gw) else ""),
            "gwas_min_mac": gw.gwas_mac.min() if len(gw) else np.nan,
            "any_block_mismatch": np.nan,      # filled below, GEA variants only
        })

    A = g.groupby("gene").apply(agg, include_groups=False).reset_index()
    A = A.drop(columns=["any_block_mismatch"])
    A = add_loci(A)

    # ---- lines of evidence ----------------------------------------------------
    A["E1_gea_multi_cluster"] = A.gea_n_clusters >= 2
    A["E2_gea_clean_axis"] = A.gea_min_lam < 2.0
    A["E3_gwas_multi_garden"] = A.gwas_n_gardens >= 2
    A["E4_both_scans"] = A.in_gea & A.in_gwas
    ecols = ["E1_gea_multi_cluster", "E2_gea_clean_axis",
             "E3_gwas_multi_garden", "E4_both_scans"]

    # E1 and E2 are NOT independent -- 48% of multi-cluster genes have a clean axis
    # against 10% of the rest, because recurring across uncorrelated climate clusters
    # is most of what it takes to land on a well-behaved axis. Summing the four flags
    # would score a purely GEA-internal result as "2 lines of evidence". So E1 and E2
    # are combined into ONE gated GEA line (recurrent AND not resting on an inflated
    # axis), and independence is counted over the three genuinely separate sources:
    # the GEA scan, the GWAS scan, and agreement between them.
    A["L_gea"] = A.E1_gea_multi_cluster & A.E2_gea_clean_axis
    A["L_gwas"] = A.E3_gwas_multi_garden
    A["L_cross"] = A.E4_both_scans
    A["n_evidence"] = A[ecols].sum(axis=1)          # kept for auditing, do not rank on
    A["n_independent"] = A[["L_gea", "L_gwas", "L_cross"]].sum(axis=1)
    A = A.sort_values(["n_independent", "gea_n_clusters", "gwas_n_gardens",
                       "gea_nlp"], ascending=False).reset_index(drop=True)
    A.to_csv(f"{OUT}/candidate_genes.csv", index=False)

    # the shortlist: at least one FULLY gated line of evidence
    C = A[A.n_independent >= 1].copy()
    # locus ids recomputed within the convergent set, so locus_n_genes counts
    # convergent genes per locus rather than pool genes per locus
    C = add_loci(C.drop(columns=["locus", "locus_n_genes"]))
    ann = annotate(C.gene.tolist()).set_index("gene")
    for col in ("symbol", "protein_name", "uniprot_function", "uniprot_keywords",
                "go_bp", "categories"):
        C[col] = C.gene.map(ann[col]) if col in ann.columns else ""
    C = C.sort_values(["n_independent", "gea_n_clusters", "gwas_n_gardens",
                       "gea_nlp"], ascending=False).reset_index(drop=True)
    C.to_csv(f"{OUT}/convergent_genes.csv", index=False)

    print("\nraw flags (correlated -- for auditing only):")
    for c in ecols:
        print(f"  {c:22s} {int(A[c].sum()):>5,}")
    print("\nindependent lines of evidence:")
    for c in ("L_gea", "L_gwas", "L_cross"):
        print(f"  {c:22s} {int(A[c].sum()):>5,}")
    print(A.n_independent.value_counts().sort_index().to_string())
    print(f"\nshortlist (>=1 gated line): {len(C):,} genes in "
          f"{C.locus.nunique():,} physical loci ({LOCUS_MERGE_BP//1000} kb merge)")
    two = C[C.n_independent >= 2]
    print(f"  of which >=2 INDEPENDENT lines: {len(two)} "
          f"({', '.join(two.gene)})")
    top = C.groupby("locus").agg(chrom=("chrom", "first"), n=("gene", "size"),
                                 lo=("gene_start", "min"), hi=("gene_end", "max"))
    top = top.sort_values("n", ascending=False).head(5)
    for r in top.itertuples():
        print(f"    {r.chrom}:{r.lo/1e6:.2f}-{r.hi/1e6:.2f} Mb  {r.n} genes")
    print(f"\nwrote {OUT}/candidate_genes.csv, {OUT}/convergent_genes.csv, "
          f"{OUT}/variants_classified.csv")


if __name__ == "__main__":
    main()
