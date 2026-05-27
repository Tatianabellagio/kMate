#!/usr/bin/env python -u
import pandas as pd
import numpy as np
from pathlib import Path

RES = Path("/carnegie/nobackup/scratch/tbellagio/hapfire_sv/panel_overlap_135_vs_82/results")

ac = pd.read_csv(RES / "ac_split_Chr1.tsv.gz", sep="\t")
pg = pd.read_csv(RES / "ac_pg151_Chr1.tsv.gz", sep="\t")
m = ac.merge(pg[["chrom", "pos", "alt_idx", "ac_pg", "an_pg"]],
             on=["chrom", "pos", "alt_idx"], how="left")
m["pg_present"] = ~m["ac_pg"].isna()
m["ac_pg"] = m["ac_pg"].fillna(0).astype(int)
m["pg_carrier"] = m["ac_pg"] > 0

m["cat"] = np.where(
    (m["ac_82"] > 0) & (m["ac_53"] > 0), "shared",
    np.where(
        (m["ac_82"] > 0) & (m["ac_53"] == 0), "cactus_only",
        np.where(
            (m["ac_82"] == 0) & (m["ac_53"] > 0),
            np.where(m["n_alt"] == 1, "V1a", "V1b"),
            "empty",
        )
    )
)
v1 = m[m["cat"].isin(["V1a", "V1b"])].copy()
v1["ac_53_bin"] = pd.cut(v1["ac_53"], bins=[0, 1, 4, 9, 100],
                         labels=["1 (singleton)", "2-4", "5-9", "10+"], include_lowest=False)

# vtype
def vt(r, a):
    if r == 1 and a == 1: return "SNP"
    d = a - r
    if d == 0: return "MNP/other"
    elif d > 0: return "INS>=50" if d >= 50 else "INS<50"
    else: return "DEL>=50" if -d >= 50 else "DEL<50"
v1["vtype"] = [vt(r, a) for r, a in zip(v1["ref_len"], v1["alt_len"])]

print("Chr1 V1 records — PG carrier rate by ac_53 bin")
print()
print("Overall (V1a + V1b):")
overall = v1.groupby("ac_53_bin", observed=True).agg(
    n_records=("pos", "size"),
    pg_present=("pg_present", "sum"),
    n_pg_carrier=("pg_carrier", "sum"),
    mean_ac_pg=("ac_pg", "mean"),
)
overall["frac_pg_carrier"] = overall["n_pg_carrier"] / overall["n_records"]
overall["frac_pg_carrier_given_present"] = overall["n_pg_carrier"] / overall["pg_present"].replace(0, np.nan)
print(overall.to_string())

print()
print("Split V1a / V1b:")
splt = v1.groupby(["cat", "ac_53_bin"], observed=True).agg(
    n_records=("pos", "size"),
    n_pg_carrier=("pg_carrier", "sum"),
    pg_present=("pg_present", "sum"),
)
splt["frac_pg_carrier"] = splt["n_pg_carrier"] / splt["n_records"]
splt["frac_pg_carrier_given_present"] = splt["n_pg_carrier"] / splt["pg_present"].replace(0, np.nan)
print(splt.to_string())

print()
print("By variant type, all V1 (no ac_53 stratification):")
vt_table = v1.groupby("vtype").agg(
    n_records=("pos", "size"),
    n_pg_carrier=("pg_carrier", "sum"),
    pg_present=("pg_present", "sum"),
)
vt_table["frac_pg_carrier"] = vt_table["n_pg_carrier"] / vt_table["n_records"]
vt_table["frac_pg_carrier_given_present"] = vt_table["n_pg_carrier"] / vt_table["pg_present"].replace(0, np.nan)
print(vt_table.to_string())
