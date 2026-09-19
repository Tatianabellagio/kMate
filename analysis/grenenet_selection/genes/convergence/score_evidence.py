#!/usr/bin/env python
"""Merge the evidence-matrix stages and count lines of evidence per variant.

Lines (each a column, each 0/1; the count is how many are satisfied -- never a p-value):
  climate      selection-aware garden gradient >= 97th pct, linear (pct_sel) OR non-linear
               (pct_quad: humps / U-shapes, so a temperate optimum is not missed)
  specific     locus-specific (genome-wide profile removed) >= 90th pct, linear or
               non-linear, AND the variant is the top-2 non-SNP record at its locus
  gwas         Bonferroni in >= 2 per-garden GWAS gardens (a different trait: founder `s`)
  mechanism    2 = strong (gene loss, truncation, promoter/5'UTR SV >= 50 bp within 2 kb of
               the TSS), 1 = moderate (smaller promoter/5'UTR indel, splice, 3'UTR), 0 else;
               the line counts at >= 1, the grade is kept
  chromatin    REF footprint removes ATAC peak sequence (regulatory positions only)
  expression   carriers vs non-carriers FDR < 0.10 AND beats the lineage control (p_emp <= 0.05)

Flags (caveats, never filters): af_shared_vector, low_support (called < 100 of 231),
few_carriers (< 5), gene_silent (median raw < 1 in 1001T), u_shape_only, mnp. Plus
snp_rivals == 0 reported as `snp_blind_locus` (SNPs do not mark it as well), which is about
kMate's added value, not about function.

Writes results/evidence_matrix.csv.
"""
import os, glob
import numpy as np, pandas as pd
from statsmodels.stats.multitest import multipletests
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = f"{HERE}/results"

A = pd.read_csv(f"{OUT}/evidence_matrix_A.csv")
B = pd.concat([pd.read_csv(f) for f in sorted(glob.glob(f"{OUT}/evidence_matrix_B_*.csv"))])
C = pd.read_csv(f"{OUT}/evidence_matrix_C.csv")
D = A.merge(B, on="store_row", how="left").merge(C, on="store_row", how="left")

t = D.p_expr.notna()
D.loc[t, "q_expr"] = multipletests(D.loc[t, "p_expr"], method="fdr_bh")[1]
D["L_climate"] = ((D.pct_sel >= 97) | (D.pct_quad >= 97)).astype(int)
D["L_specific"] = (((D.pct_resid >= 90) | (D.pct_quad_resid >= 90))
                   & (D.local_rank <= 2)).astype(int)
D["L_gwas"] = (D.gwas_n_gardens.fillna(0) >= 2).astype(int)
size = (D.alt_len - D.ref_len).abs()
reg = D["mode"].isin(["promoter", "5'UTR"])
near = D.dist_tss.abs() <= 2000
D["mech_grade"] = np.select(
    [D["mode"].isin(["gene loss", "truncation"]) | (reg & (size >= 50) & near),
     (reg & near) | D["mode"].isin(["splice", "3'UTR", "in-frame / minor coding"])],
    [2, 1], 0)
D["L_mechanism"] = (D.mech_grade >= 1).astype(int)
regulatory = D["mode"].isin(["promoter", "5'UTR", "3'UTR", "intron", "intergenic", "splice"])
D["L_chromatin"] = ((D.atac_bp > 0) & regulatory).astype(int)
D["L_expression"] = ((D.q_expr < 0.10) & (D.p_expr_emp <= 0.05)).astype(int)
L = ["L_climate", "L_specific", "L_gwas", "L_mechanism", "L_chromatin", "L_expression"]
D["n_lines"] = D[L].sum(1)

D["low_support"] = D.called < 100
D["few_carriers"] = D.carriers < 5
D["gene_silent"] = D.gene_median_raw < 1
D["u_shape_only"] = (D["shape"] == "U") & (D.pct_sel < 97)
D["mnp"] = D.vclass == "mnp"
D["snp_blind_locus"] = D.snp_rivals == 0
D.to_csv(f"{OUT}/evidence_matrix.csv", index=False)

print(f"{len(D)} variants, {D.target_gene.nunique()} genes")
print("variants per line:", {c[2:]: int(D[c].sum()) for c in L})
print("distribution of n_lines:", D.n_lines.value_counts().sort_index().to_dict())
