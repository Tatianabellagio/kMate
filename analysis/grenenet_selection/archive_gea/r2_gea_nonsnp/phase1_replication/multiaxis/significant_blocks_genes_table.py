#!/usr/bin/env python
"""Master significant-blocks-and-genes table: every (model x climate axis x class)
combo, current production regime (isotonic WZA, 3-class snp/sv/smallindel + the
pooled nonsnp for reference), at BOTH Bonferroni and BH-FDR.

One row per (model, axis, class, block): block span, block p (Z_pVal), BH q,
and the overlapping TAIR10 genes (semicolon-joined). This is the base table
everything else (class-specific / overlap / candidate-gene lists) is derived
from -- kept separate from those derived tables so it can be re-sliced any way
(by axis, by model, by class) without re-reading 240 WZA files each time.

Output:
  multiaxis/significant_blocks_bonferroni.csv
  multiaxis/significant_blocks_fdr.csv

  PY=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python
  $PY significant_blocks_genes_table.py
"""
from __future__ import annotations
import os, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))))
import lib

MA = f"{lib.GEA}/r2_gea_nonsnp/phase1_replication/results/multiaxis"
WD = f"{MA}/wza"
AXES = [f"bio{i}" for i in range(1, 20)] + ["pc1"]
MODELS = ["kendall", "lfmm", "binomial"]
CLASSES = ["snp", "sv", "smallindel", "nonsnp"]


def block_spans(r2=0.9):
    tag = f"clq{r2}"; rows = []
    for ci in range(1, 6):
        g = pd.read_csv(f"{lib.CLQ_BLOCKS_DIR}/chr{ci}_{tag}_blocks_{tag}.tsv", sep="\t"
                        ).sort_values("start_pos").reset_index(drop=True)
        for idx, r in g.iterrows():
            rows.append((f"Chr{ci}_{idx}", f"Chr{ci}", int(r.start_pos), int(r.end_pos)))
    return pd.DataFrame(rows, columns=["block", "chrom", "start", "end"]).set_index("block")


def bh(p):
    p = np.asarray(p, float); n = len(p)
    o = np.argsort(p); q = np.empty(n)
    q[o] = (p[o] * n) / (np.arange(n) + 1)
    q[o] = np.minimum.accumulate(q[o][::-1])[::-1]
    return np.clip(q, 0, 1)


def load_wza(model, cls, axis):
    f = f"{WD}/wza_{model}_{cls}_gen9_{axis}_isotonic.csv"
    if not os.path.exists(f):
        return None
    w = pd.read_csv(f).rename(columns={"index": "block"})
    w["block"] = w["block"].astype(str)
    w = w[w["Z_pVal"].notna() & (w["Z_pVal"] > 0)].copy()
    w["q"] = bh(w["Z_pVal"].to_numpy())
    return w


def main():
    spans = block_spans(); genes_df = lib.load_genes()
    gene_cache = {}

    def genes_on(block):
        if block in gene_cache:
            return gene_cache[block]
        if block not in spans.index:
            gene_cache[block] = []; return []
        s = spans.loc[block]
        gc = genes_df[(genes_df.chrom == s.chrom) & (genes_df.end >= s.start) & (genes_df.start <= s.end)]
        gene_cache[block] = list(gc.gene)
        return gene_cache[block]

    bonf_rows, fdr_rows = [], []
    missing = 0
    for model in MODELS:
        for axis in AXES:
            for cls in CLASSES:
                w = load_wza(model, cls, axis)
                if w is None:
                    missing += 1; continue
                n = len(w)
                for tier, mask, rows in [("bonf", w["Z_pVal"] < 0.05 / n, bonf_rows),
                                          ("fdr", w["q"] < 0.05, fdr_rows)]:
                    sig = w[mask]
                    for _, r in sig.iterrows():
                        b = r["block"]
                        sp = spans.loc[b] if b in spans.index else None
                        gs = genes_on(b)
                        rows.append(dict(model=model, axis=axis, cls=cls, block=b,
                                         chrom=sp.chrom if sp is not None else "",
                                         start=int(sp.start) if sp is not None else -1,
                                         end=int(sp.end) if sp is not None else -1,
                                         Z_pVal=float(r["Z_pVal"]), q=float(r["q"]),
                                         n_genes=len(gs), genes=";".join(gs)))
    print(f"missing (model,axis,cls) files: {missing} / {len(MODELS)*len(AXES)*len(CLASSES)}")

    B = pd.DataFrame(bonf_rows).sort_values(["model", "axis", "cls", "Z_pVal"])
    F = pd.DataFrame(fdr_rows).sort_values(["model", "axis", "cls", "Z_pVal"])
    B.to_csv(f"{MA}/significant_blocks_bonferroni.csv", index=False)
    F.to_csv(f"{MA}/significant_blocks_fdr.csv", index=False)
    print(f"Bonferroni: {len(B)} rows | {B.block.nunique()} unique blocks | "
          f"{B[B.n_genes>0].block.nunique()} with >=1 gene -> {MA}/significant_blocks_bonferroni.csv")
    print(f"FDR       : {len(F)} rows | {F.block.nunique()} unique blocks | "
          f"{F[F.n_genes>0].block.nunique()} with >=1 gene -> {MA}/significant_blocks_fdr.csv")

    print("\ncounts per (model, cls), summed over axes:")
    for name, df in [("Bonferroni", B), ("FDR", F)]:
        piv = df.groupby(["model", "cls"]).size().unstack("cls").reindex(index=MODELS, columns=CLASSES)
        print(f"\n{name}:\n{piv.fillna(0).astype(int).to_string()}")


if __name__ == "__main__":
    main()
