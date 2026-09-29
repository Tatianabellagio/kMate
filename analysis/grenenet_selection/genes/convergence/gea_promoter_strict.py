#!/usr/bin/env python
"""Stricter promoter definitions, and what they change.

The 1 kb promoter class is large because Arabidopsis is gene-dense: 1 kb upstream of the 28,496
genes covers 18.9 Mb -- 15.8% of the genome and 32% of everything outside gene bodies -- and the
median gap between neighbouring genes is only 909 bp. The call itself checks out (96.6% of
promoter calls sit 0-1 kb upstream of a TSS by an independent gene table; the rest are large
deletions whose SPAN reaches the window from further away). Two stricter rules are tried here:

  contained_1kb   the variant's whole REF span must lie inside the 1 kb window, not merely
                  overlap it -- so a 5 kb deletion that happens to reach a promoter no longer
                  counts as one
  overlap_500     the current overlap rule with a 500 bp window
  contained_500   both at once

Only loci the classifier currently calls promoter can change: a stricter rule cannot create
promoters. Those loci are re-tested here and, if they fail, re-assigned by the classifier's own
order -- reference TE, else within PROXIMAL_BP of a gene (proximal intergenic), else gene desert.
Gene-body tiers are untouched by construction.

Writes results/gea_universe_regions_strict.csv.gz (one tier column per rule).
env: kmate. Compute node.
"""
import os, sys
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
GEA = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, os.path.join(GEA, "genes", "dissection")); sys.path.insert(0, GEA)
import screen_sig_blocks as S                                    # noqa: E402

OUT = f"{HERE}/results"
RULES = {"contained_1kb": (1000, True), "overlap_500": (500, False), "contained_500": (500, True)}


def main():
    R = pd.read_csv(f"{OUT}/gea_universe_regions.csv.gz")
    G, _, TE = S.load_gff()
    prom = R.tier_1kb.eq("3_promoter").to_numpy()
    print(f"{len(R):,} loci; {prom.sum():,} currently called promoter (only these can change)", flush=True)
    vs = R.pos.to_numpy(np.int64); ve = vs + R.ref_len.to_numpy(np.int64) - 1
    out = {k: R.tier_1kb.to_numpy().copy() for k in RULES}
    for ch, g in G.groupby("chrom"):
        m = prom & (R.chrom.to_numpy() == ch)
        if not m.any():
            continue
        st = g.start.to_numpy(np.int64); en = g.end.to_numpy(np.int64); sd = g.strand.to_numpy(str)
        te = TE[TE.chrom == ch]
        tst, ten = te.start.to_numpy(np.int64), te.end.to_numpy(np.int64)
        idx = np.where(m)[0]
        for name, (bp, contained) in RULES.items():
            lo = np.where(sd == "+", st - bp, en + 1)
            hi = np.where(sd == "+", st - 1, en + bp)
            keep = np.zeros(len(idx), bool); te_hit = np.zeros(len(idx), bool); near = np.zeros(len(idx), np.int64)
            for j, i in enumerate(idx):
                a, b = vs[i], ve[i]
                keep[j] = ((a >= lo) & (b <= hi)).any() if contained else ((a <= hi) & (b >= lo)).any()
                if not keep[j]:
                    te_hit[j] = bool(((tst <= b) & (ten >= a)).any())
                    near[j] = np.maximum.reduce([st - b, a - en, np.zeros(len(st), np.int64)]).min()
            fall = np.where(te_hit, "6_TE",
                            np.where(near <= S.PROXIMAL_BP, "7_proximal_intergenic", "8_gene_desert"))
            out[name][idx] = np.where(keep, "3_promoter", fall)
        print(f"  {ch}: {m.sum():,} promoter loci re-tested", flush=True)
    for k in RULES:
        R[f"tier_{k}"] = out[k]
    R.to_csv(f"{OUT}/gea_universe_regions_strict.csv.gz", index=False)
    print("\npromoter loci retained by each rule (of the 1 kb overlap set):")
    for k in RULES:
        kept = (R[f"tier_{k}"] == "3_promoter").sum()
        print(f"  {k:14s} {kept:,}  ({100 * kept / prom.sum():.1f}%)")


if __name__ == "__main__":
    main()
