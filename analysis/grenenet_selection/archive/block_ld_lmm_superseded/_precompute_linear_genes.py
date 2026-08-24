#!/usr/bin/env python
"""Precompute Ensembl-Plants gene annotation for the rose-MORE-than-ecotype blocks of the
LINEAR + sampling-floor + drift-null model (site{SITE}). Significant set = boundary-guarded
(panel_freq in [0.10,0.90]) blocks with drift-null one-sided q_pos_drift < QCUT, joined to the
estimator effects (s, linked_bg, resid) on the unique haplotype id gid. Run in `basic` env
(outbound HTTPS, urllib). Writes site{SITE}_linear_rosemore_genes.csv (load-only for the notebook).
"""
import os, sys
import pandas as pd

ROOT = "/global/scratch/users/tbellg/kmate"
sys.path.insert(0, f"{ROOT}/analysis/grenenet_gea")
from genes_from_regions import genes_in_region, label  # noqa: E402

SITE = int(os.environ.get("SITE", 4))
QCUT = float(os.environ.get("QCUT", 0.05))
H = f"{ROOT}/results/grenenet_gea/hapfreq"
OUT = f"{H}/site{SITE}_linear_rosemore_genes.csv"

dn = pd.read_csv(f"{H}/site{SITE}_block_ld_lmm_driftnull_linear.csv")
mod = pd.read_csv(f"{H}/site{SITE}_block_ld_lmm_linear_sampvar.csv")
d = dn.merge(mod[["gid", "s", "linked_bg", "resid"]], on="gid", how="left")   # unique-key join
sig = d[(d.guard == True) & (d.q_pos_drift < QCUT)].sort_values("p_pos_marg").reset_index(drop=True)
print(f"site {SITE}: {len(sig)} rose-more blocks (guarded, q_pos_drift<{QCUT}); "
      f"{int((sig.p_pos_fwer < 0.05).sum())} also genome-wide FWER (p_pos_fwer<0.05)")

rows = []
for _, r in sig.iterrows():
    g = genes_in_region(r.chrom, int(r.unit_start), int(r.unit_end))
    named = g[g.symbol != ""] if len(g) else g
    gene_descs = "; ".join(
        f"{(row.symbol or row.gene_id)}: {row.description}".strip(": ")
        for _, row in g.iterrows()) if len(g) else ""
    rows.append({
        "fwer": "★" if r.p_pos_fwer < 0.05 else "",
        "chrom": r.chrom, "start": int(r.unit_start), "end": int(r.unit_end),
        "panel_freq": float(r.panel_freq), "s": float(r.s), "linked_bg": float(r.linked_bg),
        "resid": float(r.resid), "z": float(r.z),
        "p_pos_fwer": float(r.p_pos_fwer), "q_pos_drift": float(r.q_pos_drift),
        "n_genes": int(len(g)),
        "symbols": ";".join(named.symbol.tolist()) if len(named) else "",
        "genes": label(g) if len(g) else "",
        "gene_descriptions": gene_descs,
    })

out = pd.DataFrame(rows)
out.to_csv(OUT, index=False)
print(f"[done] {len(out)} blocks ({int((out.n_genes > 0).sum())} overlap >=1 gene) -> {OUT}")
