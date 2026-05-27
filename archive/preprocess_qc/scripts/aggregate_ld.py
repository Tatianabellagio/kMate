"""
aggregate_ld.py — bin raw LD pair TSVs into distance × (class1, class2) summary.

Reads preprocess_qc/output/ld/{grenenet_snp,merged_snp,merged_sv}_chr{1-5}.tsv,
streams pairs (no in-memory concat), bins each pair by log-distance and by
(class1, class2), accumulates running median/quartile/count statistics, and
writes preprocess_qc/output/ld/ld_summary.tsv:

  source  chrom  class1  class2  bin_lo  bin_hi  bin_mid  n  median_r2  q25_r2  q75_r2  mean_r2

Bin edges: log-spaced 1..1e6 bp, 25 bins.

This compresses ~25 GB of raw pairs to ~50 KB summary.
"""
import sys, glob
import numpy as np
import pandas as pd
from collections import defaultdict
from pathlib import Path

LD_DIR = Path("preprocess_qc/output/ld")
OUT    = LD_DIR / "ld_summary.tsv"
BINS = np.logspace(0, 6, 26)
BIN_LO = BINS[:-1]; BIN_HI = BINS[1:]

# (source, chrom, class1, class2, bin_idx) -> list of r²
buckets = defaultdict(list)

sources = {
    "grenenet_snp": "GrENE-Net 231 (SNP-only catalog)",
    "merged_snp":   "Production merged 231 (SNP anchors)",
    "merged_sv":    "Production merged 231 (SV anchors)",
}

for src in sources:
    for chrom in range(1, 6):
        path = LD_DIR / f"{src}_chr{chrom}.tsv"
        if not path.exists():
            print(f"  SKIP missing: {path}", file=sys.stderr); continue
        print(f"  reading {path}", file=sys.stderr)
        # Stream in chunks to avoid OOM on 2 GB files
        for chunk in pd.read_csv(path, sep="\t", chunksize=2_000_000):
            d = chunk.distance.values
            r2 = chunk.r2.values
            c1 = chunk.class1.values
            c2 = chunk.class2.values
            bins = np.clip(np.searchsorted(BIN_HI, d, side="right"), 0, len(BIN_HI)-1)
            for i in range(len(d)):
                key = (src, chrom, c1[i], c2[i], int(bins[i]))
                buckets[key].append(r2[i])

print(f"\n  buckets: {len(buckets):,}", file=sys.stderr)
rows = []
for (src, chrom, c1, c2, bi), vals in buckets.items():
    arr = np.asarray(vals, dtype=np.float32)
    rows.append({
        "source":     sources[src],
        "source_key": src,
        "chrom":      chrom,
        "class1":     c1,
        "class2":     c2,
        "bin_lo":     BIN_LO[bi],
        "bin_hi":     BIN_HI[bi],
        "bin_mid":    np.sqrt(BIN_LO[bi] * BIN_HI[bi]),
        "n":          len(arr),
        "median_r2":  float(np.median(arr)),
        "q25_r2":     float(np.quantile(arr, 0.25)),
        "q75_r2":     float(np.quantile(arr, 0.75)),
        "mean_r2":    float(arr.mean()),
    })

summary = pd.DataFrame(rows)
summary.to_csv(OUT, sep="\t", index=False)
print(f"  wrote {len(summary):,} rows -> {OUT}", file=sys.stderr)
print(f"  sources/chroms covered:", file=sys.stderr)
for (s, c), n in summary.groupby(["source_key","chrom"]).n.sum().items():
    print(f"    {s} chr{c}: {int(n):,} pairs", file=sys.stderr)
