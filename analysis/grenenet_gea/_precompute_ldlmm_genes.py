#!/usr/bin/env python
"""Precompute the gene-annotation table for the FDR-significant (q<0.05) temporal LD-LMM
blocks, so the notebook stays load-only (no network at render).

Takes the q<0.05 subset of site{SITE}_block_ld_lmm_temporal.csv, queries Ensembl Plants
(via genes_from_regions.genes_in_region) for every overlapping gene, and writes one row per
significant block: chrom,start,end,q,resid(effect),n_genes,symbols, and the gene
symbols+descriptions joined for the notebook table. Run in the `basic` env on a compute node
(needs outbound HTTPS).
"""
import os
import sys
import pandas as pd

ROOT = "/global/scratch/users/tbellg/kmate"
sys.path.insert(0, f"{ROOT}/analysis/grenenet_gea")
from genes_from_regions import genes_in_region, label  # noqa: E402

SITE = int(os.environ.get("SITE", 4))
H = f"{ROOT}/results/grenenet_gea/hapfreq"
OUT = f"{H}/site{SITE}_ldlmm_temporal_genes.csv"

d = pd.read_csv(f"{H}/site{SITE}_block_ld_lmm_temporal.csv")
bonf = 0.05 / len(d)                                   # Bonferroni threshold on p
sig = d[d.q < 0.05].sort_values("p").reset_index(drop=True)
print(f"site {SITE}: {len(sig)} FDR-significant (q<0.05) blocks of {len(d)} testable; "
      f"{int((sig.p < bonf).sum())} also Bonferroni-significant (p<{bonf:.1e})")

rows = []
for _, r in sig.iterrows():
    g = genes_in_region(r.chrom, int(r.unit_start), int(r.unit_end))
    named = g[g.symbol != ""] if len(g) else g
    # symbol-or-id : description pairs (factual Ensembl descriptions only)
    gene_descs = "; ".join(
        f"{(row.symbol or row.gene_id)}: {row.description}".strip(": ")
        for _, row in g.iterrows()
    ) if len(g) else ""
    rows.append({
        "bonferroni": "★" if r.p < bonf else "",   # most-significant (p < 0.05/M)
        "chrom": r.chrom,
        "start": int(r.unit_start),
        "end": int(r.unit_end),
        "p": float(r.p),
        "q": float(r.q),
        "resid": float(r.resid),          # block-specific selection (effect)
        "s": float(r.s),
        "n_genes": int(len(g)),
        "symbols": ";".join(named.symbol.tolist()) if len(named) else "",
        "genes": label(g) if len(g) else "",
        "gene_descriptions": gene_descs,
    })

out = pd.DataFrame(rows)
out.to_csv(OUT, index=False)
print(f"[done] {len(out)} significant blocks ({int((out.n_genes > 0).sum())} overlap >=1 gene) -> {OUT}")
