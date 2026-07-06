#!/usr/bin/env python
"""Genes overlapping the top blocks of the 31-site varlen Manhattan, for BOTH panels.

Top panel    = raw per-block test   (block_gea.csv)  -> inflated, ranked by block_p then |beta1|
Bottom panel = calibrated block-WZA (block_wza.csv)  -> ranked by block_p (most extreme |WZA|)

For each top block we list every TAIR10 protein-coding gene overlapping its
[unit_start, unit_end] span (AT-IDs from the local GFF). Run in kmate env; writes
top_block_genes_raw.csv / _wza.csv into PB_DIR and prints a readable table. AT->symbol
resolution is done in a second pass (resolve_symbols.py via Ensembl Plants REST).
"""
import os, sys
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lib

H = "results/grenenet_gea/hapfreq"
B = os.environ.get("PB_DIR", f"{H}/pipelineB_varlen")
TOPN = int(os.environ.get("TOPN", 15))


def genes_in_span(genes, chrom, lo, hi):
    gc = genes[genes.chrom == chrom]
    ov = gc[(gc.start <= hi) & (gc.end >= lo)]
    return ov.gene.tolist()


def annotate_blocks(df, genes, score_col, extra_cols):
    rows = []
    for _, r in df.iterrows():
        gl = genes_in_span(genes, r.chrom, r.unit_start, r.unit_end)
        rows.append({
            "unit": r.unit,
            "chrom": r.chrom,
            "start": int(r.unit_start),
            "end": int(r.unit_end),
            "span_kb": round((r.unit_end - r.unit_start) / 1e3, 2),
            "block_p": r.block_p,
            "q": r.q,
            **{c: r[c] for c in extra_cols},
            "n_genes": len(gl),
            "genes": ";".join(gl),
        })
    return pd.DataFrame(rows)


def main():
    genes = lib.load_genes()

    raw = pd.read_csv(f"{B}/block_gea.csv")
    raw["abs_b"] = raw.best_beta1.abs()
    raw_top = raw.sort_values(["block_p", "abs_b"], ascending=[True, False]).head(TOPN)
    raw_ann = annotate_blocks(raw_top, genes, "block_p",
                              ["n_hap_tested", "best_cluster", "best_beta1", "best_s_mean"])
    raw_ann.to_csv(f"{B}/top_block_genes_raw.csv", index=False)

    wza = pd.read_csv(f"{B}/block_wza.csv")
    wza["abs_w"] = wza.wza.abs()
    wza_top = wza.sort_values(["block_p", "abs_w"], ascending=[True, False]).head(TOPN)
    wza_ann = annotate_blocks(wza_top, genes, "block_p",
                              ["n_hap", "wza", "lead_cluster", "lead_beta1"])
    wza_ann.to_csv(f"{B}/top_block_genes_wza.csv", index=False)

    pd.set_option("display.width", 200, "display.max_colwidth", 80)
    print("=" * 100)
    print(f"TOP {TOPN} — RAW per-block test (top Manhattan panel; lambda=8.58, inflated)")
    print("=" * 100)
    print(raw_ann[["unit", "span_kb", "block_p", "q", "best_beta1", "n_genes", "genes"]]
          .to_string(index=False))
    print("\n" + "=" * 100)
    print(f"TOP {TOPN} — calibrated block-WZA (bottom Manhattan panel; lambda=1.03; none clear FDR)")
    print("=" * 100)
    print(wza_ann[["unit", "span_kb", "block_p", "q", "wza", "n_genes", "genes"]]
          .to_string(index=False))
    # unique AT-IDs for the symbol-resolution pass
    allg = sorted(set(g for s in pd.concat([raw_ann.genes, wza_ann.genes])
                      for g in s.split(";") if g))
    with open(f"{B}/top_block_atids.txt", "w") as fh:
        fh.write("\n".join(allg))
    print(f"\n[done] {len(allg)} unique AT-IDs -> {B}/top_block_atids.txt "
          f"(resolve symbols next)")


if __name__ == "__main__":
    main()
