#!/usr/bin/env python -u
"""
Of the 151 PG founders' total ALT-allele mass, what fraction goes to:
  - shared (ALT carried by cactus_82 AND extras_53)
  - cactus_only (ALT carried by cactus_82 only)
  - V1a (extras-introduced biallelic bubble)
  - V1b (extras-introduced ALT at multi-allelic bubble)
  - empty (ALT carried by neither in pang_135; PG sometimes still calls due to ./.→. noise)

stratified by variant type (SNP / INS<50 / DEL<50 / INS≥50 / DEL≥50 / MNP/other).

Mass = sum over (record, alt_idx) of ac_pg = total ALT-allele calls across the
151 PG founders. This is what counts toward variation "carried by the 151."
"""
import pandas as pd
import numpy as np
from pathlib import Path

RES = Path(__file__).resolve().parents[1] / "results"

def vt(r, a):
    if r == 1 and a == 1: return "SNP"
    d = a - r
    if d == 0: return "MNP/other"
    elif d > 0: return "INS>=50" if d >= 50 else "INS<50"
    else: return "DEL>=50" if -d >= 50 else "DEL<50"

frames = []
for c in ["Chr1", "Chr2", "Chr3", "Chr4", "Chr5"]:
    ac_f = RES / f"ac_split_{c}.tsv.gz"
    pg_f = RES / f"ac_pg151_{c}.tsv.gz"
    if not ac_f.exists() or not pg_f.exists():
        continue
    print(f"[load] {c}")
    ac = pd.read_csv(ac_f, sep="\t")
    pg = pd.read_csv(pg_f, sep="\t")
    m = ac.merge(pg[["chrom", "pos", "alt_idx", "ac_pg"]],
                 on=["chrom", "pos", "alt_idx"], how="left")
    m["ac_pg"] = m["ac_pg"].fillna(0).astype(int)
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
    m["vtype"] = [vt(r, a) for r, a in zip(m["ref_len"], m["alt_len"])]
    frames.append(m[["cat", "vtype", "ac_pg"]])

df = pd.concat(frames, ignore_index=True)
print(f"\nTotal per-ALT rows: {len(df):,}")
print(f"Total PG ALT mass across 151 founders: {df['ac_pg'].sum():,}")

# Sum ALT mass per (cat, vtype)
mass = df.groupby(["cat", "vtype"])["ac_pg"].sum().unstack(fill_value=0)
mass["TOTAL"] = mass.sum(axis=1)
mass.loc["TOTAL"] = mass.sum(axis=0)
print("\n=== PG ALT mass (sum of ac_pg across 151 founders), by category × vtype ===")
print(mass.to_string())

# Fractions: per-vtype, what fraction of mass goes to each category
pct = mass.div(mass.loc["TOTAL"], axis=1) * 100
print("\n=== PER-VTYPE % of PG ALT mass by category (column-wise, sums to 100) ===")
print(pct.round(2).to_string())

# Fold V1a + V1b → "V1 (extras-introduced)"
df["cat2"] = df["cat"].replace({"V1a": "V1_extras_introduced", "V1b": "V1_extras_introduced"})
mass2 = df.groupby(["cat2", "vtype"])["ac_pg"].sum().unstack(fill_value=0)
mass2["TOTAL"] = mass2.sum(axis=1)
mass2.loc["TOTAL"] = mass2.sum(axis=0)
pct2 = mass2.div(mass2.loc["TOTAL"], axis=1) * 100
print("\n=== Same, with V1a+V1b folded ('V1_extras_introduced') ===")
print(mass2.to_string())
print("\nPER-VTYPE %:")
print(pct2.round(2).to_string())

mass.to_csv(RES / "pg_mass_attribution.tsv", sep="\t")
pct.to_csv(RES / "pg_mass_attribution_pct.tsv", sep="\t")
print(f"\n→ {RES}/pg_mass_attribution.tsv")
print(f"→ {RES}/pg_mass_attribution_pct.tsv")
