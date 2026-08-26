#!/usr/bin/env python
"""The master candidate-gene table: one row per gene, everything we know about it.

Consolidates every layer this section produces into a single wide table, so a candidate
can be judged without opening five CSVs:

  identity      gene, symbol, protein name, coordinates, strand
  locus         physical locus id + how many candidate genes share it (pileup warning)
  discovery     which scan found it (LFMM GEA / GEMMA per-garden GWAS / both), which
                class scans, which climate axes and clusters, which gardens
  blocks        the clq0.9 tiling block(s) its variants sit in, and the CARK-type
                block/gene mismatch flag
  variants      how many, of what class, how big, where they land, lead variant
  mechanism     strongest functional tier, frameshift / CDS-SV flags
  snp visibility  statistical shadowing AND founder-panel tagging r2 (different things)
  evidence      the three independent lines and their total
  annotation    TAIR GO + UniProt for EVERY candidate gene, not just the shortlist
  caveats       inferred-allele flag, untestable-tagging flag

Gene set = every gene a pooled variant lands **in or beside**: genic/promoter from
`variants_classified.csv`, plus the genes reached only through the downstream-flank tier
added by `build_functional_variants.py`.

Blocks. The GEA records carry their clq0.9 tiling block; GWAS-only variants do not, so
their block is recomputed from position with the same partition. `block_gene_mismatch`
reproduces the CARK check -- the block's attributed genes do not include the gene the
variant physically sits in, i.e. block-level attribution would have named the wrong
gene. It is a property of the *block*, so it is reported per gene as "any of this gene's
variants sits in a block that would have been attributed elsewhere".

Outputs -> results/
  master_candidate_genes.csv   every candidate gene, all columns
  master_candidate_genes.md    the shortlist (>=1 gated line), human-readable

env: kmate.  Run last, after the other four builders. Needs outbound HTTPS.
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
import build_convergence as bc                                # noqa: E402

OUT = f"{HERE}/results"
GENIC = bc.GENIC
FTIER_ORDER = ["F1_CDS_frameshift", "F2_CDS_inframe", "F3_UTR", "F4_promoter",
               "F5_downstream", "F6_intron"]


def assign_blocks(V: pd.DataFrame) -> pd.DataFrame:
    """Fill the clq0.9 tiling block for variants that have none (GWAS-only rows)."""
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(HERE)), "blocks"))
    import blocks_tiling as bt                                # noqa: E402
    V = V.copy()
    need = V.block.isna() | V.block.eq("")
    if not need.any():
        return V
    for ch in lib.CHROMS:
        m = need & V.chrom.eq(ch)
        if not m.any():
            continue
        b = bt.load_blocks(0.9, ch)
        ends = b.end_pos.to_numpy(np.int64)
        # block i spans (ends[i-1], ends[i]]; searchsorted gives the index directly
        idx = np.searchsorted(ends, V.loc[m, "pos"].to_numpy(np.int64))
        idx = np.clip(idx, 0, len(ends) - 1)
        V.loc[m, "block"] = [f"{ch}_{i}" for i in idx]
    return V


def per_gene(V: pd.DataFrame, F: pd.DataFrame) -> pd.DataFrame:
    """Collapse the variant table to one row per gene."""
    # V already carries ftier/frameshift/cds_sv (merged in main); do not re-merge
    fk = ["chrom", "pos", "size"]

    # a variant contributes to its own gene, and (for downstream-flank rows) to the
    # gene it sits behind -- so build the gene-variant edge list explicitly
    # V already carries every column, so the downstream edge is just V rows relabelled
    # with the gene they sit behind -- no merge, hence no column loss.
    g1 = V[V.tier.isin(GENIC) & V.gene.ne("")].copy()
    g1["target_gene"] = g1.gene

    dmap = (F[F.ftier.eq("F5_downstream")]
            .drop_duplicates(subset=fk)
            .set_index(fk)["target_gene"])
    key = pd.MultiIndex.from_frame(V[fk])
    down = V[key.isin(dmap.index)].copy()
    down["target_gene"] = pd.MultiIndex.from_frame(down[fk]).map(dmap)
    down = down[down.target_gene.notna() & down.target_gene.ne("")]

    E = pd.concat([g1, down], ignore_index=True)
    E["ftier_rank"] = E.ftier.map({t: i for i, t in enumerate(FTIER_ORDER)})

    def one(d: pd.DataFrame) -> pd.Series:
        gea, gw = d[d.in_gea], d[d.in_gwas]
        lead = d.sort_values(["gea_nlp", "gwas_nlp"], ascending=False).iloc[0]
        best_f = d.ftier.dropna()
        best_f = (d.loc[d.ftier_rank.idxmin(), "ftier"]
                  if d.ftier_rank.notna().any() else "")
        return pd.Series({
            "chrom": d.chrom.iat[0],
            # ---- discovery ------------------------------------------------
            "found_by": ("GEA+GWAS" if len(gea) and len(gw)
                         else ("GEA" if len(gea) else "GWAS")),
            "gea_classes": ",".join(sorted({c for s in gea.gea_classes.dropna()
                                            for c in str(s).split(",")})),
            "gwas_classes": ",".join(sorted({c for s in gw.gwas_classes.dropna()
                                             for c in str(s).split(",")})),
            "gea_nlp": gea.gea_nlp.max() if len(gea) else np.nan,
            "gea_best_axis": (gea.sort_values("gea_nlp", ascending=False)
                              .gea_best_axis.iat[0] if len(gea) else ""),
            "gea_n_axes": int(gea.gea_n_axes.max()) if len(gea) else 0,
            "gea_sig_axes": ",".join(sorted({a for s in gea.gea_sig_axes.dropna()
                                             for a in str(s).split(",")})),
            "gea_n_clusters": int(gea.gea_n_clusters.max()) if len(gea) else 0,
            "gea_sig_clusters": ",".join(sorted({a for s in
                                                 gea.gea_sig_clusters.dropna()
                                                 for a in str(s).split(",")})),
            "gea_min_lam": gea.gea_min_lam.min() if len(gea) else np.nan,
            "gwas_nlp": gw.gwas_nlp.max() if len(gw) else np.nan,
            "gwas_n_gardens": int(gw.gwas_n_gardens.max()) if len(gw) else 0,
            "gwas_gardens": ",".join(sorted({a for s in gw.gwas_gardens.dropna()
                                             for a in str(s).split(",")}, key=str)),
            "gwas_min_mac": gw.gwas_mac.min() if len(gw) else np.nan,
            "gwas_bio1_min": gw.gwas_bio1_min.min() if len(gw) else np.nan,
            "gwas_bio1_max": gw.gwas_bio1_max.max() if len(gw) else np.nan,
            # ---- blocks ---------------------------------------------------
            "blocks": ",".join(sorted(set(d.block.dropna().astype(str)) - {""})),
            "n_blocks": d.block.replace("", np.nan).nunique(),
            "block_gene_mismatch": bool(d.block_gene_mismatch.any())
            if "block_gene_mismatch" in d else False,
            # ---- variants -------------------------------------------------
            "n_variants": len(d),
            "n_sv": int((d.vclass == "sv").sum()),
            "n_smallindel": int((d.vclass == "smallindel").sum()),
            "n_mnp": int((d.vclass == "mnp").sum()),
            "max_size": int(d["size"].max()),
            "regions": ",".join(sorted(set(d.region.dropna()))),
            "lead_variant": f"{lead.chrom}:{int(lead.pos)}",
            "lead_size": int(lead["size"]),
            "lead_vclass": lead.vclass,
            # ---- mechanism ------------------------------------------------
            "best_ftier": best_f,
            "has_frameshift": bool(d.frameshift.fillna(False).any()),
            "has_cds_sv": bool(d.cds_sv.fillna(False).any()),
            # ---- SNP visibility -------------------------------------------
            "snp_cosig_any": bool(d.snp_cosig_2kb.fillna(False).any())
            if "snp_cosig_2kb" in d else False,
            "snp_cosig_min_dist": (d.snp_cosig_dist.min()
                                   if "snp_cosig_dist" in d else np.nan),
            "best_r2_snp": d.best_r2_snp.max() if "best_r2_snp" in d else np.nan,
            "snp_blind_any": bool((d.best_r2_snp < 0.2).any())
            if "best_r2_snp" in d else False,
            "tag_untestable_all": bool(d.best_r2_snp.isna().all())
            if "best_r2_snp" in d else True,
            # ---- caveats --------------------------------------------------
            "size_inferred_any": bool(d.size_inferred.fillna(False).any())
            if "size_inferred" in d else False,
        })

    return E.groupby("target_gene").apply(one, include_groups=False).reset_index() \
            .rename(columns={"target_gene": "gene"})


def main():
    V = pd.read_csv(f"{OUT}/variants_classified.csv")
    F = pd.read_csv(f"{OUT}/functional_variants.csv")
    A = pd.read_csv(f"{OUT}/candidate_genes.csv")

    # bring the SNP-visibility + inferred-size columns onto the variant table
    keep = [c for c in ("snp_cosig_2kb", "snp_cosig_dist", "snp_cosig_axes",
                        "best_r2_snp", "size_inferred", "ftier", "frameshift",
                        "cds_sv") if c in F.columns]
    # one row per variant key -- F can hold several rows for one variant (a variant in
    # two genes' promoters), and merging without collapsing would fan V out
    Fk = (F[["chrom", "pos", "size"] + keep]
          .sort_values("ftier", na_position="last")
          .drop_duplicates(subset=["chrom", "pos", "size"], keep="first"))
    n_before = len(V)
    V = V.merge(Fk, on=["chrom", "pos", "size"], how="left")
    assert len(V) == n_before, f"merge fanned out {n_before} -> {len(V)}"
    if "block" not in V.columns:
        V["block"] = ""
    V = assign_blocks(V)

    # CARK-type block/gene mismatch, reusing the dissection implementation
    sm = bc._screen_mod()
    G = lib.load_genes()
    print("computing block spans + block/gene mismatch ...", flush=True)
    V = sm.block_span_genes(V, G)
    print(f"  {int(V.block_gene_mismatch.sum()):,} of {len(V):,} variants sit in a "
          f"block whose attributed genes exclude their own gene")

    M = per_gene(V, F)
    print(f"\n{len(M):,} candidate genes")

    # evidence + locus from the convergence layer
    ev = A.set_index("gene")
    for c in ("E1_gea_multi_cluster", "E2_gea_clean_axis", "E3_gwas_multi_garden",
              "E4_both_scans", "L_gea", "L_gwas", "L_cross", "n_independent",
              "n_evidence"):
        M[c] = M.gene.map(ev[c]) if c in ev.columns else False
    M["n_independent"] = M.n_independent.fillna(0).astype(int)

    # coordinates, then loci over the FULL candidate set
    Gi = G.set_index("gene")
    M["gene_start"] = M.gene.map(Gi.start)
    M["gene_end"] = M.gene.map(Gi.end)
    M["strand"] = M.gene.map(Gi.strand)
    M = bc.add_loci(M.assign(chrom=M.chrom))
    M["locus_span_kb"] = (M.groupby("locus").gene_end.transform("max")
                          - M.groupby("locus").gene_start.transform("min")) / 1000.0

    print(f"annotating all {len(M):,} genes (TAIR GO + UniProt) ...", flush=True)
    ann = bc.annotate(M.gene.tolist()).set_index("gene")
    for col in ("symbol", "protein_name", "uniprot_function", "uniprot_keywords",
                "uniprot_reviewed", "go_bp", "categories"):
        M[col] = M.gene.map(ann[col]) if col in ann.columns else ""

    order = ["gene", "symbol", "protein_name", "chrom", "gene_start", "gene_end",
             "strand", "locus", "locus_n_genes", "locus_span_kb",
             "n_independent", "L_gea", "L_gwas", "L_cross",
             "E1_gea_multi_cluster", "E2_gea_clean_axis", "E3_gwas_multi_garden",
             "E4_both_scans", "n_evidence",
             "found_by", "gea_classes", "gea_nlp", "gea_best_axis", "gea_n_axes",
             "gea_sig_axes", "gea_n_clusters", "gea_sig_clusters", "gea_min_lam",
             "gwas_classes", "gwas_nlp", "gwas_n_gardens", "gwas_gardens",
             "gwas_min_mac", "gwas_bio1_min", "gwas_bio1_max",
             "blocks", "n_blocks", "block_gene_mismatch",
             "n_variants", "n_sv", "n_smallindel", "n_mnp", "max_size", "regions",
             "lead_variant", "lead_size", "lead_vclass",
             "best_ftier", "has_frameshift", "has_cds_sv",
             "snp_cosig_any", "snp_cosig_min_dist", "best_r2_snp", "snp_blind_any",
             "tag_untestable_all", "size_inferred_any",
             "uniprot_function", "uniprot_keywords", "uniprot_reviewed", "go_bp",
             "categories"]
    order = [c for c in order if c in M.columns]
    M["_ft"] = M.best_ftier.map({t: i for i, t in enumerate(FTIER_ORDER)}).fillna(9)
    M = M.sort_values(["n_independent", "gea_n_clusters", "gwas_n_gardens", "_ft"],
                      ascending=[False, False, False, True])
    M[order].to_csv(f"{OUT}/master_candidate_genes.csv", index=False)

    print(f"\nfound_by:"); print(M.found_by.value_counts().to_string())
    print(f"\nblock/gene mismatch (block attribution would name another gene): "
          f"{int(M.block_gene_mismatch.sum()):,}")
    print(f"independent lines:"); print(M.n_independent.value_counts()
                                        .sort_index().to_string())
    print(f"\nwrote {OUT}/master_candidate_genes.csv  "
          f"({len(M):,} genes x {len(order)} columns)")

    # human-readable shortlist
    S = M[M.n_independent >= 1]
    with open(f"{OUT}/master_candidate_genes.md", "w") as fh:
        fh.write("# Candidate genes carrying at least one gated line of evidence\n\n")
        fh.write(f"{len(S)} genes in {S.locus.nunique()} loci, of "
                 f"{len(M)} candidates total. Full table: "
                 f"`master_candidate_genes.csv`.\n\n")
        fh.write("| gene | sym | locus (n genes) | found by | GEA clusters/axes | "
                 "gardens | mechanism | lead variant | r2 SNP | lines |\n")
        fh.write("|---|---|---|---|---|---|---|---|---|---|\n")
        for r in S.itertuples():
            sym = r.symbol if isinstance(r.symbol, str) and r.symbol else "-"
            r2 = "-" if pd.isna(r.best_r2_snp) else f"{r.best_r2_snp:.2f}"
            fh.write(f"| {r.gene} | {sym} | {r.locus} ({r.locus_n_genes}) | "
                     f"{r.found_by} | {r.gea_n_clusters}/{r.gea_n_axes} | "
                     f"{r.gwas_n_gardens} | {r.best_ftier} | {r.lead_variant} "
                     f"({r.lead_vclass} {r.lead_size}bp) | {r2} | "
                     f"{r.n_independent} |\n")
    print(f"wrote {OUT}/master_candidate_genes.md ({len(S)} shortlist genes)")


if __name__ == "__main__":
    main()
