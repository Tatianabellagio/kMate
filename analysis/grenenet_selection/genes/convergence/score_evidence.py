#!/usr/bin/env python
"""Merge the evidence-matrix stages and count lines of evidence per variant.

Lines (each a column, each 0/1; the count is how many are satisfied -- never a p-value):
  climate      selection-aware garden gradient >= 97th pct, linear (pct_sel) OR non-linear
               (pct_quad: humps / U-shapes, so a temperate optimum is not missed)
  locus        the variant is the top-2 non-SNP record within +-25 kb on its best axis
  gwas         per-garden GWAS, judged directionally: allele moves the way GEMMA's Z says
               in >= 60% of its significant gardens' pools, more than elsewhere
  mechanism    2 = strong (gene loss, truncation, promoter/5'UTR SV >= 50 bp within 2 kb of
               the TSS), 1 = moderate (smaller promoter/5'UTR indel, splice, 3'UTR), 0 else;
               the line counts at >= 1, the grade is kept
  chromatin    REF footprint removes ATAC peak sequence (regulatory positions only)
  motif        >= 1 TF binding site lost or gained (FIMO, repeat-tract matches excluded;
               promoter / 5'UTR variants only -- others are not scanned)
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
# motif turnover (tfbs_pool.py) and the directional per-garden GWAS evaluation
if os.path.exists(f"{OUT}/tfbs_pool_summary.csv"):
    D = D.merge(pd.read_csv(f"{OUT}/tfbs_pool_summary.csv"), on="store_row", how="left")
W = pd.read_csv(f"{OUT}/gwas_evaluation.csv")[["store_row", "consistency_dir", "excess",
                                                 "n_sig_gardens"]]
W = W.sort_values("consistency_dir", ascending=False).drop_duplicates("store_row")
D = D.merge(W, on="store_row", how="left")

t = D.p_expr.notna()
D.loc[t, "q_expr"] = multipletests(D.loc[t, "p_expr"], method="fdr_bh")[1]
D["L_climate"] = ((D.pct_sel >= 97) | (D.pct_quad >= 97)).astype(int)
# LOCUS, not "locus-specific": removing the genome-wide garden components over-corrects.
# Unlinked candidates on one axis share garden profiles, but their founder carrier sets
# overlap at chance -- a parallel response of different founders' alleles to the same
# climate, i.e. real polygenic signal (memory raw-lfmm-over-wza-decision: do not correct
# it away). What attributes a signal to THIS locus is being its best marker locally.
# pct_resid / pct_quad_resid stay in the table as information.
D["L_locus"] = (D.local_rank <= 2).astype(int)
# GWAS judged directionally on its own criterion (gwas_evaluation.py): the allele moves
# the way GEMMA's Z says, in most pools of its significant garden(s), more than elsewhere
D["L_gwas"] = ((D.consistency_dir >= 0.6) & (D.excess > 0)).astype(int)
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
# a motif call that is repeat-driven (same motif matched at many offsets, or an indel inside a
# microsatellite -- repeat_context.py) is not counted as a line; tfbs_repeat stays as a flag
turnover = (D.get("tfbs_lost", 0).fillna(0) + D.get("tfbs_gained", 0).fillna(0)) > 0
D["tfbs_repeat"] = D.get("tfbs_repeat", False).fillna(False).astype(bool)
D["L_motif"] = (turnover & ~D.tfbs_repeat).astype(int)
L = ["L_climate", "L_locus", "L_gwas", "L_mechanism", "L_chromatin", "L_motif", "L_expression"]
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
