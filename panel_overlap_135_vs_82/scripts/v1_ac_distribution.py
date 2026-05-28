#!/usr/bin/env python -u
"""V1 ALT records — distribution of ac_53 (number of the 53 extras carrying the ALT).
Singletons (ac_53==1) are higher-risk for noise (sequencing artifact, het-mask edge case).
Multi-carrier V1 (ac_53>=2) are more likely real variation."""
import gzip
import pandas as pd
from pathlib import Path

RES = Path(__file__).resolve().parents[1] / "results"

dfs = []
for c in ["Chr1", "Chr2", "Chr3", "Chr4", "Chr5"]:
    f = RES / f"v1_records_{c}.tsv.gz"
    if not f.exists():
        continue
    df = pd.read_csv(f, sep="\t")
    dfs.append(df)
df = pd.concat(dfs, ignore_index=True)

print(f"Total V1 records (genome-wide): {len(df):,}")
print(f"  V1a: {(df['cat']=='V1a').sum():,}")
print(f"  V1b: {(df['cat']=='V1b').sum():,}")

print("\nac_53 distribution (number of the 53 extras carrying the ALT):")
ac_dist = df["ac_53"].value_counts().sort_index()
print(ac_dist.head(15).to_string())
print(f"  ...")
print(f"  >= 5 carriers: {(df['ac_53']>=5).sum():,}")
print(f"  >= 10 carriers: {(df['ac_53']>=10).sum():,}")
print(f"  == 1 (singleton): {(df['ac_53']==1).sum():,}")
print(f"  fraction singleton: {(df['ac_53']==1).mean():.3f}")

print("\nSingleton rate by V1 sub-category × variant type:")
df["is_singleton"] = df["ac_53"] == 1
grp = df.groupby(["cat", "vtype"]).agg(
    n=("pos", "size"),
    n_singleton=("is_singleton", "sum"),
)
grp["frac_singleton"] = grp["n_singleton"] / grp["n"]
print(grp.to_string())
grp.to_csv(RES / "tier1_v1_ac_distribution.tsv", sep="\t")
