#!/usr/bin/env python -u
"""
Final headline report — once PG merge has emitted ac_pg151_Chr*.tsv.gz for
all chroms, this combines V1 records with PG carrier counts and emits the
single-table summary that answers the user's question.

Question: for the 151 short-read founders, what fraction of pang_135 records
that exist only because of the 53 extras are actually USED by PanGenie (i.e.
≥1 of the 151 is called as ALT carrier)?

Output: tier1_headline_combined.tsv  +  one printed summary.
"""
import gzip
import sys
import pandas as pd
import numpy as np
from pathlib import Path

RES = Path("/carnegie/nobackup/scratch/tbellagio/hapfire_sv/panel_overlap_135_vs_82/results")

def load_chrom(chrom):
    ac_f = RES / f"ac_split_{chrom}.tsv.gz"
    pg_f = RES / f"ac_pg151_{chrom}.tsv.gz"
    if not ac_f.exists() or not pg_f.exists():
        return None
    # check the PG file is complete by checking last byte sequence isn't 0
    # actually just try and catch
    print(f"[load] {chrom}")
    try:
        ac = pd.read_csv(ac_f, sep="\t", dtype={
            "chrom": "string", "pos": "int64",
            "ref_len": "int64", "n_alt": "int32", "alt_idx": "int32",
            "alt_len": "int64",
            "ac_82": "int32", "an_82": "int32",
            "ac_53": "int32", "an_53": "int32"})
        pg = pd.read_csv(pg_f, sep="\t", dtype={
            "chrom": "string", "pos": "int64",
            "ref_len": "int64", "n_alt": "int32", "alt_idx": "int32",
            "alt_len": "int64",
            "ac_pg": "int32", "an_pg": "int32"})
    except (EOFError, OSError) as e:
        print(f"[skip {chrom}] {type(e).__name__}: {e} — file likely still being written")
        return None
    m = ac.merge(pg[["chrom", "pos", "alt_idx", "ac_pg", "an_pg"]],
                 on=["chrom", "pos", "alt_idx"], how="left")
    m["pg_present"] = ~m["ac_pg"].isna()
    m["ac_pg"] = m["ac_pg"].fillna(0).astype(int)
    m["an_pg"] = m["an_pg"].fillna(0).astype(int)
    m["pg_carrier"] = m["ac_pg"] > 0
    # add cat
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
    # vtype
    def vt(r, a):
        if r == 1 and a == 1: return "SNP"
        d = a - r
        if d == 0: return "MNP/other"
        elif d > 0: return "INS>=50" if d >= 50 else "INS<50"
        else: return "DEL>=50" if -d >= 50 else "DEL<50"
    m["vtype"] = [vt(r, a) for r, a in zip(m["ref_len"], m["alt_len"])]
    # add ac_53 bin for V1 records
    m["ac_53_bin"] = pd.cut(m["ac_53"], bins=[-1, 0, 1, 4, 9, 100], labels=["0", "1", "2-4", "5-9", "10+"])
    return m

def main():
    chroms = []
    for c in ["Chr1", "Chr2", "Chr3", "Chr4", "Chr5"]:
        m = load_chrom(c)
        if m is not None:
            chroms.append(m)
        else:
            print(f"[skip] {c} not ready")
    if not chroms:
        sys.exit("no inputs")
    df = pd.concat(chroms, ignore_index=True)
    print(f"\nGenome-wide rows: {len(df):,}")

    # ========== Headline 1: of V1 records, what fraction has ≥1 PG carrier? ==========
    head = df.groupby(["cat", "vtype"]).agg(
        n_records=("pos", "size"),
        pg_present=("pg_present", "sum"),
        n_pg_carrier=("pg_carrier", "sum"),
        mean_ac_pg=("ac_pg", "mean"),
    )
    head["frac_pg_carrier"] = head["n_pg_carrier"] / head["n_records"]
    head["frac_pg_present"] = head["pg_present"] / head["n_records"]
    print("\n=== Headline 1: PG-carrier fraction by category × vtype ===")
    print(head.to_string())
    head.to_csv(RES / "tier1_headline_combined.tsv", sep="\t")

    # ========== Headline 2: V1 records only, stratified by ac_53 bin ==========
    v1 = df[df["cat"].isin(["V1a", "V1b"])]
    head2 = v1.groupby(["cat", "ac_53_bin"], observed=True).agg(
        n_records=("pos", "size"),
        n_pg_carrier=("pg_carrier", "sum"),
        mean_ac_pg=("ac_pg", "mean"),
    )
    head2["frac_pg_carrier"] = head2["n_pg_carrier"] / head2["n_records"]
    print("\n=== Headline 2: V1 ALTs by ac_53 bin (more extras carrying it → more likely PG also carries it) ===")
    print(head2.to_string())
    head2.to_csv(RES / "tier1_headline_v1_by_ac53.tsv", sep="\t")

    # ========== Headline 3: V1 ALTs by vtype, multi-carrier only (ac_53 ≥ 2) ==========
    v1m = v1[v1["ac_53"] >= 2]
    head3 = v1m.groupby("vtype").agg(
        n_records=("pos", "size"),
        n_pg_carrier=("pg_carrier", "sum"),
        mean_ac_pg=("ac_pg", "mean"),
    )
    head3["frac_pg_carrier"] = head3["n_pg_carrier"] / head3["n_records"]
    print("\n=== Headline 3: V1 ALTs with ac_53 >= 2, by vtype (the 'realistic' value-add) ===")
    print(head3.to_string())

    # ========== Top-line single number ==========
    n_v1 = len(v1)
    n_v1_pg = (v1["pg_carrier"]).sum()
    print("\n=== TOP-LINE ANSWER ===")
    print(f"Total V1 records (pang_135 ALTs not in 82-cactus): {n_v1:,}  (= {100*n_v1/len(df):.2f}% of all per-ALT records)")
    print(f"  Of these, {n_v1_pg:,} ({100*n_v1_pg/n_v1:.2f}%) have ≥1 of the 151 PG founders carrying the ALT.")
    print(f"  → That's the 'realized' value-add of the 53 extras for PanGenie genotyping the 151.")
    v1_strict = v1[v1["ac_53"] >= 2]
    print(f"\nIf we restrict to V1 with ≥2 extras carriers (less noise):")
    print(f"  {len(v1_strict):,} V1 records, {(v1_strict['pg_carrier']).sum():,} ({100*v1_strict['pg_carrier'].mean():.2f}%) PG-carried.")

if __name__ == "__main__":
    main()
