#!/usr/bin/env python
"""Augment the SV/indel functional screen with a SNP co-signal flag.

The screen already includes ALL significant SV/indel/non-SNP block-leads (it never
applied the newpeak SNP-missed filter), so SV/indel candidates that share a locus
with SNPs are already present. This adds, per candidate, whether a genome-wide
per-class-Bonferroni-significant SNP sits within +/-2 kb (same LFMM axis basis),
so one can separate SNP-shadowed candidates from SV/indel-unique ones without
dropping either.

Reads sig_block_functional_screen.csv, adds:
  snp_cosig_2kb   bool  -- a Bonferroni SNP within 2 kb (any axis)
  snp_cosig_dist  int   -- distance (bp) to the nearest such SNP (NaN if none)
  snp_cosig_axes  str   -- axes in which that nearby SNP is significant
Writes back in place (+ prints the SNP-shadowed vs unique split).
env: kmate.  Run on a compute node.
"""
from __future__ import annotations
import os, sys
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(HERE, "..", "..")))
import lib

WZAIN = f"{lib.GEA}/r2_gea_nonsnp/phase1_replication/results/multiaxis/wza_in_clq09_tile"
AXES = [f"bio{i}" for i in range(1, 20)] + ["pc1"]
WIN = 2000
SNP_BONF_NLP = 7.606                      # 0.05 / n_snp(MAF>0.05), genome-wide


def sig_snps():
    """All genome-wide per-class-Bonferroni SNP hits across axes -> (chrom,pos,axis)."""
    rows = []
    for ax in AXES:
        f = f"{WZAIN}/lfmm_snp_gen9_{ax}.csv"
        if not os.path.exists(f):
            continue
        d = pd.read_csv(f, usecols=["chrom", "pos", "MAF", "pval"])
        d = d[d.MAF > 0.05]
        nlp = -np.log10(d.pval.clip(lower=1e-300))
        s = d[nlp > SNP_BONF_NLP]
        for ch, pos in zip(s.chrom.to_numpy(str), s.pos.to_numpy(int)):
            rows.append((ch, pos, ax))
    return pd.DataFrame(rows, columns=["chrom", "pos", "axis"])


def main():
    C = pd.read_csv(f"{HERE}/sig_block_functional_screen.csv")
    S = sig_snps()
    print(f"genome-wide Bonferroni SNP hits (MAF>0.05, any axis): {len(S)} "
          f"({S[['chrom','pos']].drop_duplicates().shape[0]} unique loci)")

    by_chrom = {ch: d.sort_values("pos") for ch, d in S.groupby("chrom")}
    cosig, dist, axes = [], [], []
    for _, r in C.iterrows():
        d = by_chrom.get(r.chrom)
        if d is None:
            cosig.append(False); dist.append(np.nan); axes.append(""); continue
        near = d[(d.pos >= r.pos - WIN) & (d.pos <= r.pos + WIN)]
        if len(near):
            cosig.append(True)
            dist.append(int((near.pos - r.pos).abs().min()))
            axes.append(",".join(sorted(near.axis.unique())))
        else:
            cosig.append(False); dist.append(np.nan); axes.append("")
    C["snp_cosig_2kb"] = cosig
    C["snp_cosig_dist"] = dist
    C["snp_cosig_axes"] = axes
    C.to_csv(f"{HERE}/sig_block_functional_screen.csv", index=False)

    fp = C[C.tier.str.startswith(("1", "2", "3"))]
    print(f"\nfunctional candidates (CDS/UTR/promoter): {len(fp)}")
    print(f"  with a genome-wide SNP within 2 kb (SNP-shadowed): {int(fp.snp_cosig_2kb.sum())}")
    print(f"  SV/indel-UNIQUE (no Bonferroni SNP within 2 kb):   {int((~fp.snp_cosig_2kb).sum())}")
    print("\n  SV-only, by SNP co-signal:")
    sv = fp[fp.vclass == "sv"]
    print(f"    shadowed {int(sv.snp_cosig_2kb.sum())} | unique {int((~sv.snp_cosig_2kb).sum())}")


if __name__ == "__main__":
    main()
