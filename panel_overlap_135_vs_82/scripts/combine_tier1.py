#!/usr/bin/env python -u
"""Combine per-chrom Tier 1 summary tables into a genome-wide table."""
import pandas as pd
from pathlib import Path

RES = Path("/carnegie/nobackup/scratch/tbellagio/hapfire_sv/panel_overlap_135_vs_82/results")

dfs = []
for c in ["Chr1", "Chr2", "Chr3", "Chr4", "Chr5"]:
    f = RES / f"tier1_summary_{c}.tsv"
    if not f.exists():
        print(f"missing {f}")
        continue
    df = pd.read_csv(f, sep="\t")
    # The chrom col is a single value spread, but the index col is 'cat'
    # Re-read properly
    df = pd.read_csv(f, sep="\t", index_col=0)
    df["chrom"] = c
    df = df.reset_index().set_index(["chrom", "cat"])
    dfs.append(df)

if not dfs:
    raise SystemExit("no inputs")
combined = pd.concat(dfs)
# Sum across chroms keeping the "cat" rows aligned
# Drop "TOTAL" row to recompute cleanly
combined = combined[combined.index.get_level_values("cat") != "TOTAL"]
combined = combined.drop(columns=[c for c in combined.columns if c in ("chrom",)], errors="ignore")

genome = combined.groupby("cat").sum()
genome["TOTAL_per_cat"] = genome.sum(axis=1)

# Top-line genome-wide
print("Genome-wide per-(record, ALT) row counts, by category × variant type")
print(genome.to_string())

# Fractions
total_all = genome.loc[["V1a", "V1b", "cactus_only", "empty", "shared"], "TOTAL_per_cat"].sum()
print(f"\nTotal per-ALT records across pang_135 (5 chroms): {total_all:,}")
v1_total = genome.loc[["V1a", "V1b"], "TOTAL_per_cat"].sum()
print(f"V1 (extras-introduced ALTs) total: {v1_total:,}  ({100*v1_total/total_all:.2f}%)")
print(f"  V1a (new biallelic bubble): {genome.loc['V1a', 'TOTAL_per_cat']:,}  ({100*genome.loc['V1a','TOTAL_per_cat']/total_all:.2f}%)")
print(f"  V1b (new ALT at multi-allelic): {genome.loc['V1b', 'TOTAL_per_cat']:,}  ({100*genome.loc['V1b','TOTAL_per_cat']/total_all:.2f}%)")
print(f"shared (≥1 cactus AND ≥1 extras): {genome.loc['shared', 'TOTAL_per_cat']:,}  ({100*genome.loc['shared','TOTAL_per_cat']/total_all:.2f}%)")
print(f"cactus_only (≥1 cactus, 0 extras): {genome.loc['cactus_only', 'TOTAL_per_cat']:,}  ({100*genome.loc['cactus_only','TOTAL_per_cat']/total_all:.2f}%)")
print(f"empty: {genome.loc['empty', 'TOTAL_per_cat']:,}")

# Save combined
out = RES / "tier1_summary_combined.tsv"
genome.to_csv(out, sep="\t")
print(f"\n→ {out}")

# Same numbers by variant type
print("\n--- V1 fraction by variant type (genome-wide):")
vtypes = ["SNP", "INS<50", "DEL<50", "INS>=50", "DEL>=50", "MNP/other"]
vt_table = pd.DataFrame(index=vtypes, columns=["V1a", "V1b", "cactus_only", "shared", "empty"])
for v in vtypes:
    for cat in vt_table.columns:
        if cat in genome.index and v in genome.columns:
            vt_table.loc[v, cat] = int(genome.loc[cat, v])
vt_table = vt_table.fillna(0).astype(int)
vt_table["TOTAL"] = vt_table.sum(axis=1)
vt_table["V1_frac"] = (vt_table["V1a"] + vt_table["V1b"]) / vt_table["TOTAL"]
vt_table["V1a_frac"] = vt_table["V1a"] / vt_table["TOTAL"]
vt_table["V1b_frac"] = vt_table["V1b"] / vt_table["TOTAL"]
print(vt_table.to_string())
vt_table.to_csv(RES / "tier1_summary_by_vtype.tsv", sep="\t")
print(f"\n→ {RES / 'tier1_summary_by_vtype.tsv'}")
