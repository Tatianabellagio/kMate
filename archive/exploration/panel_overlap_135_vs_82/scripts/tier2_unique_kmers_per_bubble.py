#!/usr/bin/env python -u
"""
Tier 2 D1: per-bubble unique-k-mer counts from PanGenie's pang_135 index files.

For each chrom, parse pang_135_pangenie_index_Chr{N}_kmers.tsv.gz:
  cols: #chromosome  start  end  unique_kmers  unique_kmers_overhang
where unique_kmers and unique_kmers_overhang are comma-separated 31-mer lists
(total across all alleles in the bubble; per-allele mapping lives in the binary
.cereal index).

Per-bubble counts of unique k-mers tell us about D1 (the more-accessions →
fewer-unique-k-mers-per-bubble downside). The metric of interest is the
distribution shape — particularly how many bubbles have very few (or zero)
unique kmers, since PanGenie's signal degrades when bubbles run out of unique
discriminating k-mers under its 16/allele cap.

Also: cross-reference with Tier 1's V1 records to ask whether extras-only
bubbles have more or fewer unique k-mers than shared/cactus-only bubbles.
"""
import gzip
import sys
from pathlib import Path

import numpy as np
import pandas as pd

BASE = Path(__file__).resolve().parents[2]
PG_INDEX = BASE / "pangenie_genotyping" / "data"
WORK = BASE / "panel_overlap_135_vs_82"
RES = WORK / "results"

def parse_one(chrom):
    f = PG_INDEX / f"pang_135_pangenie_index_{chrom}_kmers.tsv.gz"
    rows = []
    with gzip.open(f, "rt") as h:
        header = next(h)
        for line in h:
            parts = line.rstrip("\n").split("\t")
            # cols: chrom, start, end, kmers, overhang
            chrom_, start, end, kmers, overhang = parts
            n_kmers = 0 if kmers in ("", "nan") else kmers.count(",") + 1
            n_overhang = 0 if overhang in ("", "nan") else overhang.count(",") + 1
            rows.append((chrom_, int(start), int(end), n_kmers, n_overhang))
    df = pd.DataFrame(rows, columns=["chrom", "start", "end", "n_unique_kmers", "n_overhang_kmers"])
    df["bubble_len"] = df["end"] - df["start"]
    return df

def stratify_by_v1(df, chrom):
    """Mark each bubble: contains_v1 (any V1 record falls inside)."""
    v1_f = RES / f"v1_records_{chrom}.tsv.gz"
    if not v1_f.exists():
        return df
    v1 = pd.read_csv(v1_f, sep="\t")
    # v1.pos is the VCF POS (1-based). PanGenie bubble [start, end) is 0-based half-open
    # (relative to the bubble's position-before-ref convention).
    # For matching purposes use a conservative interval containment: V1 belongs to bubble
    # if bubble.start <= v1.pos <= bubble.end (loose).
    df = df.sort_values(["chrom", "start", "end"]).reset_index(drop=True)
    v1 = v1.sort_values(["chrom", "pos"]).reset_index(drop=True)
    # For each V1 record, find the bubble whose start <= pos <= end via vectorized search
    df["contains_v1a"] = False
    df["contains_v1b"] = False
    df["n_v1"] = 0
    # Build a per-chrom array of (start, end) for binary search
    starts = df["start"].to_numpy()
    ends = df["end"].to_numpy()
    pos = v1["pos"].to_numpy()
    cats = v1["cat"].to_numpy()
    # for each v1 pos, find the bubble idx via searchsorted on starts
    idx = np.searchsorted(starts, pos, side="right") - 1
    # check pos in [start, end]
    valid = (idx >= 0) & (pos <= ends[np.clip(idx, 0, len(ends) - 1)])
    bidx = idx[valid]
    bcat = cats[valid]
    # accumulate counts
    n_v1 = np.bincount(bidx, minlength=len(df))
    has_v1a = np.bincount(bidx[bcat == "V1a"], minlength=len(df)) > 0
    has_v1b = np.bincount(bidx[bcat == "V1b"], minlength=len(df)) > 0
    df["n_v1"] = n_v1
    df["contains_v1a"] = has_v1a
    df["contains_v1b"] = has_v1b
    return df

def main():
    out_all = []
    for chrom in ["Chr1", "Chr2", "Chr3", "Chr4", "Chr5"]:
        f = PG_INDEX / f"pang_135_pangenie_index_{chrom}_kmers.tsv.gz"
        if not f.exists():
            print(f"[skip] {f} missing")
            continue
        print(f"[parse] {chrom}")
        df = parse_one(chrom)
        print(f"  {len(df):,} bubbles, mean n_uniq={df['n_unique_kmers'].mean():.1f}, median={df['n_unique_kmers'].median()}")
        df = stratify_by_v1(df, chrom)
        out_all.append(df)
        out_per_chrom = RES / f"tier2_bubble_kmers_{chrom}.tsv.gz"
        df.to_csv(out_per_chrom, sep="\t", index=False, compression="gzip")
        print(f"  → {out_per_chrom}")
    if not out_all:
        return
    combined = pd.concat(out_all, ignore_index=True)
    # Per-strata summary
    print("\n[summary] unique-k-mer counts stratified by bubble V1 status (V1 = contains an extras-introduced ALT)")
    grp = combined.groupby([combined["contains_v1a"], combined["contains_v1b"]]).agg(
        n_bubbles=("n_unique_kmers", "size"),
        mean_n_kmers=("n_unique_kmers", "mean"),
        median_n_kmers=("n_unique_kmers", "median"),
        p5_n_kmers=("n_unique_kmers", lambda x: x.quantile(0.05)),
        p95_n_kmers=("n_unique_kmers", lambda x: x.quantile(0.95)),
        frac_zero=("n_unique_kmers", lambda x: (x == 0).mean()),
    )
    print(grp.to_string())
    grp.to_csv(RES / "tier2_bubble_kmer_summary.tsv", sep="\t")
    print(f"\n→ {RES}/tier2_bubble_kmer_summary.tsv")

if __name__ == "__main__":
    main()
