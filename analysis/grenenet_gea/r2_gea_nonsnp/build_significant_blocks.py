#!/usr/bin/env python
"""List significant WZA blocks (A=dp, B=scoef) + every gene each block overlaps.

Takes the WZA block results (build_wza.py), computes BH-FDR + Bonferroni on the
SNP-number-corrected block p (Z_pVal), reconstructs each block's genomic SPAN from
its member SVs, and annotates with ALL TAIR10 genes overlapping that span.

Outputs (analysis/grenenet_gea/gea/wza/):
  significant_blocks.csv   every block with q<FDR_THRESH OR Bonferroni, both stats,
        cols: stat, block, chrom, start, end, span_bp, n_sv, SNPs, Z, Z_pVal, fdr,
              bonferroni, dir, n_up, n_dn, n_genes, genes, gene_names, flower_loci
  wza_{stat}_bio1.fdr.csv  full per-block table with fdr column (all blocks)

  PY=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python
  $PY analysis/grenenet_gea/r2_gea_nonsnp/build_significant_blocks.py --fdr 0.10
"""
from __future__ import annotations
import argparse, glob, os, sys
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import lib

GEA = lib.GEA; STORE = lib.AF_STORE; WZ = f"{GEA}/gea/wza"
FLOWER = {"AT4G00650": "FRI", "AT5G10140": "FLC", "AT1G65480": "FT",
          "AT2G45660": "SOC1", "AT5G61850": "LFY"}


def bh_fdr(p):
    p = np.asarray(p, float); n = np.isfinite(p).sum()
    q = np.full(len(p), np.nan)
    ok = np.where(np.isfinite(p))[0]
    o = ok[np.argsort(p[ok])]
    ranked = p[o] * n / (np.arange(1, len(o) + 1))
    q[o] = np.minimum.accumulate(ranked[::-1])[::-1].clip(max=1)
    return q


def block_spans():
    """Genomic span (chrom,start,end,n_sv) of every LD block, from its SVs."""
    idx = np.load(f"{STORE}/index_nonsnp.npz")
    size = np.abs(idx["alt_len"].astype(np.int64) - idx["ref_len"].astype(np.int64))
    nc = np.asarray(np.load(sorted(glob.glob(f"{STORE}/nc_nonsnp/*.npy"))[0]))
    mask = (size > 50) & (nc >= 150)
    block = np.load(f"{GEA}/gea/twostage_blocks.npz", allow_pickle=True)["block"].astype(str)
    z = np.load(f"{GEA}/gea/twostage_dp_bio1.npz", allow_pickle=True)
    d = pd.DataFrame(dict(block=block, chrom=z["chrom"].astype(str), pos=z["pos"]))
    d = d[d.block != ""]
    g = d.groupby("block").agg(chrom=("chrom", "first"), start=("pos", "min"),
                               end=("pos", "max"), n_sv=("pos", "size")).reset_index()
    g["span_bp"] = g.end - g.start
    return g


def genes_in_span(chrom, start, end, genes):
    gc = genes[(genes.chrom == chrom) & (genes.end >= start) & (genes.start <= end)]
    return gc.gene.tolist(), gc.name.fillna("").tolist()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fdr", type=float, default=0.10)
    ap.add_argument("--climate", default="bio1")
    args = ap.parse_args()
    genes = lib.load_genes()
    spans = block_spans()

    sig_rows = []
    for stat in ("dp", "scoef"):
        w = pd.read_csv(f"{WZ}/wza_{stat}_{args.climate}.csv")
        w["fdr"] = bh_fdr(w.Z_pVal.values)
        nb = w.Z_pVal.notna().sum()
        w["bonferroni"] = w.Z_pVal < 0.05 / nb
        w.to_csv(f"{WZ}/wza_{stat}_{args.climate}.fdr.csv", index=False)
        sigw = w[(w.fdr < args.fdr) | w.bonferroni].drop(
            columns=[c for c in ("chrom", "mid_pos") if c in w.columns])
        sig = sigw.merge(spans, on="block", how="left")
        for _, r in sig.iterrows():
            gids, gnames = genes_in_span(r.chrom, r.start, r.end, genes)
            fl = [FLOWER[g] for g in gids if g in FLOWER]
            sig_rows.append(dict(
                stat=stat, block=r.block, chrom=r.chrom, start=int(r.start),
                end=int(r.end), span_bp=int(r.span_bp), n_sv=int(r.n_sv),
                SNPs=int(r.SNPs), Z=round(r.Z, 3), Z_pVal=r.Z_pVal, fdr=round(r.fdr, 4),
                bonferroni=bool(r.bonferroni), dir=r.dir, n_up=int(r.n_up), n_dn=int(r.n_dn),
                n_genes=len(gids), genes=";".join(gids),
                gene_names=";".join(n for n in gnames if n), flower_loci=";".join(fl)))
        print(f"[{stat}] {nb:,} blocks | FDR<{args.fdr}: {(w.fdr<args.fdr).sum()} | "
              f"Bonferroni: {int(w.bonferroni.sum())}")
    out = pd.DataFrame(sig_rows).sort_values(["stat", "Z_pVal"])
    path = f"{WZ}/significant_blocks.csv"
    out.to_csv(path, index=False)
    print(f"\n-> {len(out)} significant block-records (both stats) written to {path}")
    print(f"   (unique blocks: {out.block.nunique()}; flowering-locus hits: "
          f"{(out.flower_loci!='').sum()})")


if __name__ == "__main__":
    main()
