#!/usr/bin/env python
"""Gene-level summary: one row per significant gene, with function + how many models hit it.

Collapses the per-(block,gene) annotated table to ONE row per gene, unioning evidence
across blocks (2 genes span >1 significant block). For each gene reports:
  - function: symbol, entrez, name, summary, GO, functional categories, climate flag
  - DETECTION: how many of the 3 MODELS (kendall/lfmm/binomial) and 3 variant CLASSES
    (snp/smallindel/sv) and 9 model x class COMBOS flagged it BH-significant, plus the
    lists and the best (smallest) BH q-value.

Reads:  significant_genes_annotated_gen{g}_{clim}_{regime}.csv
Writes: gene_summary_gen{g}_{clim}_{regime}.csv  (sorted: most models, then most combos)

  PY=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python
  $PY analysis/grenenet_selection/r2_gea_nonsnp/phase1_replication/build_gene_summary.py
"""
from __future__ import annotations
import argparse, os, sys
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import lib

MODELS = ["kendall", "lfmm", "binomial"]
CLASSES = ["snp", "smallindel", "sv"]
COMBOS = [(m, c) for m in MODELS for c in CLASSES]
COMBO_COLS = [f"{m}_{c}" for m, c in COMBOS]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--regime", default="deg7cap2000")
    ap.add_argument("--climate", default="bio1")
    ap.add_argument("--gen", type=int, default=9)
    args = ap.parse_args()
    base = f"{lib.GEA}/r2_gea_nonsnp/phase1_replication/results"
    src = f"{base}/significant_genes_annotated_gen{args.gen}_{args.climate}_{args.regime}.csv"
    d = pd.read_csv(src)
    d = d[d.gene.notna() & (d.gene.astype(str) != "")].copy()
    for c in COMBO_COLS:                       # numeric q (NaN where not significant)
        d[c] = pd.to_numeric(d[c], errors="coerce")

    rows = []
    for gene, g in d.groupby("gene"):
        flagged = [(m, c) for (m, c) in COMBOS if g[f"{m}_{c}"].notna().any()]
        models = sorted({m for m, _ in flagged}, key=MODELS.index)
        classes = sorted({c for _, c in flagged}, key=CLASSES.index)
        qmin = np.nanmin(g[COMBO_COLS].to_numpy()) if flagged else np.nan
        f0 = g.iloc[0]
        # exact genomic region(s) per block (browser-pastable, TAIR10), unique per gene
        reg = (g.drop_duplicates("block")
                .apply(lambda r: f"{r.chrom}:{int(r.start)}-{int(r.end)}", axis=1).tolist())
        rows.append(dict(
            gene=gene, symbol=f0.get("symbol", ""), entrez=f0.get("entrez", ""),
            fn_name=f0.get("fn_name", ""),
            categories=f0.get("categories", ""),
            climate_stress_flowering=bool(f0.get("climate_stress_flowering", False)),
            n_models=len(models), models=";".join(models),
            n_classes=len(classes), classes=";".join(classes),
            n_combos=len(flagged), combos=";".join(f"{m}:{c}" for m, c in flagged),
            best_q=float(qmin) if np.isfinite(qmin) else None,   # already floored upstream
            q_floored=bool(g["q_floored"].any()) if "q_floored" in g else False,
            n_blocks=g.block.nunique(), blocks=";".join(sorted(g.block.unique())),
            region=";".join(reg),
            go_bp=f0.get("go_bp", ""), summary=str(f0.get("summary", ""))[:240]))
    out = pd.DataFrame(rows).sort_values(
        ["n_models", "n_combos", "best_q"], ascending=[False, False, True]).reset_index(drop=True)
    path = f"{base}/gene_summary_gen{args.gen}_{args.climate}_{args.regime}.csv"
    out.to_csv(path, index=False)

    print(f"{len(out)} genes | columns: {list(out.columns)}\n")
    print("detection breakdown (how many models per gene):")
    print(out.n_models.value_counts().sort_index(ascending=False).to_string())
    print(f"\nclimate/stress/flowering genes: {int(out.climate_stress_flowering.sum())}")
    print("\n== genes hit by ALL 3 models (n_models==3) ==")
    top = out[out.n_models == 3]
    print(top[["symbol", "gene", "categories", "n_models", "n_classes",
               "n_combos", "best_q", "fn_name"]].to_string(index=False))
    print(f"\n-> {path}")


if __name__ == "__main__":
    main()
