#!/usr/bin/env python
"""v2 step 2 -- gene-level candidate table for the NEW list: GEA only (events + call-rate 0.50), no GWAS arm.

Same logic as `build_convergence.py` (classify every pooled variant at its OWN position via the TAIR10 GFF (promoter BEFORE TE for this list); assign it to the
gene it lands in or its promoter; aggregate per gene; E1/E2 -> one gated GEA line; physical loci) with the GWAS pool, the cross-scan
union and the `E3/E4 / L_gwas / L_cross` lines REMOVED -- not just left empty. There is no `gwas_*` column anywhere in the output.

Input : results/candidate_list_v2_events_cr50/gea_pool.csv   (step 1; one row per UNIT: short-indel/MNP record or SV event)
Output (all in results/candidate_list_v2_events_cr50/, `v2_` prefixed so a file name never collides with the old list):
  v2_variants_classified.csv  every pooled unit with its genomic context
  v2_candidate_genes.csv      one row per gene a pooled unit lands in or beside; GEA evidence only
  v2_gea_gated_genes.csv      genes with the gated GEA line (recurrent over >= 2 climate clusters AND a clean axis), annotated
Lines of evidence (GEA only):  E1 gea_multi_cluster (>= 2 of the 7 climate clusters), E2 gea_clean_axis (some axis with lambda < 2.0),
L_gea = E1 & E2.  Raw uncalibrated LFMM p; this is a candidate net, not a calibrated hit list.
env: kmate. Compute node (parses the TAIR10 GFF; the annotator needs outbound HTTPS).
"""
from __future__ import annotations
import os, sys
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__)); CONV = os.path.dirname(HERE)
sys.path.insert(0, CONV); sys.path.insert(0, os.path.abspath(os.path.join(CONV, "..", "..")))
import build_convergence as bc                                  # noqa: E402  (helpers only: classify module, loci, annotation)
V2 = os.path.join(CONV, "results", "candidate_list_v2_events_cr50")


def load_gea_pool() -> pd.DataFrame:
    G = pd.read_csv(f"{V2}/gea_pool.csv")
    G["chrom"] = bc.norm_chrom(G.chrom)
    G = G.rename(columns={"best_nlp": "gea_nlp", "best_axis": "gea_best_axis", "min_lam": "gea_min_lam",
                          "n_axes": "gea_n_axes", "n_clusters": "gea_n_clusters", "sig_axes": "gea_sig_axes",
                          "sig_clusters": "gea_sig_clusters", "cleanest_axis": "gea_cleanest_axis",
                          "has_clean_axis": "gea_clean_axis", "sig_classes": "gea_classes", "MAF": "gea_MAF"})
    cols = ["unit", "chrom", "pos", "ref_len", "alt_len", "size", "vclass", "gea_MAF", "gea_best_axis", "gea_nlp",
            "gea_min_lam", "gea_n_axes", "gea_n_clusters", "gea_sig_axes", "gea_sig_clusters", "gea_cleanest_axis",
            "gea_clean_axis", "gea_classes", "block", "n_records"]
    U = G[cols].copy()
    assert U.unit.is_unique
    return U


def main():
    U = load_gea_pool()
    print(f"pooled GEA units: {len(U):,}  ({U.vclass.value_counts().to_dict()})")
    sm = bc._screen_mod()
    sm.TE_ABOVE_PROMOTER = False        # this list only (user, 2026-10-02): a promoter call wins over a TE overlap; te_overlap is kept as a column
    print("classifying against TAIR10 GFF (promoter <= 1 kb takes precedence over TE) ...", flush=True)
    V, _ = sm.classify(U)
    assert "unit" in V.columns and V.unit.is_unique and len(V) == len(U), "classification must keep one row per unit"
    V.to_csv(f"{V2}/v2_variants_classified.csv", index=False)
    print(V.tier.value_counts().sort_index().to_string())

    g = V[V.tier.isin(bc.GENIC) & V.gene.ne("")].copy()
    print(f"\n{len(g):,} units land in a gene or its promoter -> {g.gene.nunique():,} genes")

    def agg(d: pd.DataFrame) -> pd.Series:
        return pd.Series({
            "chrom": d.chrom.iat[0], "n_units": len(d),
            "n_sv_events": int((d.vclass == "sv").sum()), "n_smallindel": int((d.vclass == "smallindel").sum()),
            "regions": ",".join(sorted(set(d.region))), "vclasses": ",".join(sorted(set(d.vclass.dropna()))),
            "max_size": int(d["size"].max()),
            "gea_nlp": d.gea_nlp.max(), "gea_n_clusters": int(d.gea_n_clusters.max()),
            "gea_n_axes": int(d.gea_n_axes.max()), "gea_min_lam": d.gea_min_lam.min(),
            "gea_sig_clusters": ",".join(sorted({c for s in d.gea_sig_clusters for c in str(s).split(",")})),
            "units": ",".join(map(str, sorted(d.unit))),
        })

    A = g.groupby("gene").apply(agg, include_groups=False).reset_index()
    A = bc.add_loci(A)
    A["E1_gea_multi_cluster"] = A.gea_n_clusters >= 2
    A["E2_gea_clean_axis"] = A.gea_min_lam < 2.0
    A["L_gea"] = A.E1_gea_multi_cluster & A.E2_gea_clean_axis       # one gated line: E1 and E2 are not independent
    A = A.sort_values(["L_gea", "gea_n_clusters", "gea_nlp"], ascending=False).reset_index(drop=True)
    A.to_csv(f"{V2}/v2_candidate_genes.csv", index=False)

    C = A[A.L_gea].copy()
    C = bc.add_loci(C.drop(columns=["locus", "locus_n_genes"]))     # loci recomputed within the gated set
    ann = bc.annotate(C.gene.tolist()).set_index("gene")
    for col in ("symbol", "protein_name", "uniprot_function", "uniprot_keywords", "go_bp", "categories"):
        C[col] = C.gene.map(ann[col]) if col in ann.columns else ""
    C = C.sort_values(["gea_n_clusters", "gea_nlp"], ascending=False).reset_index(drop=True)
    C.to_csv(f"{V2}/v2_gea_gated_genes.csv", index=False)

    print(f"\ngenes with any pooled unit: {len(A):,} | E1 {int(A.E1_gea_multi_cluster.sum()):,} | E2 {int(A.E2_gea_clean_axis.sum()):,} "
          f"| L_gea (E1&E2) {int(A.L_gea.sum()):,}")
    print(f"gated set: {len(C):,} genes in {C.locus.nunique():,} physical loci ({bc.LOCUS_MERGE_BP // 1000} kb merge)")
    print(f"wrote {V2}/v2_variants_classified.csv, v2_candidate_genes.csv, v2_gea_gated_genes.csv")


if __name__ == "__main__":
    main()
