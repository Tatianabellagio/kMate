#!/usr/bin/env python
"""
Tier 1 analysis — classify per-ALT records into:
  shared    : ac_82 > 0  AND  ac_53 > 0
  cactus_only: ac_82 > 0  AND  ac_53 == 0
  V1        : ac_82 == 0 AND  ac_53 >  0  ("extras-introduced ALT")
  empty     : ac_82 == 0 AND  ac_53 == 0  (after subsetting; shouldn't have many)

Within V1, split:
  V1a : biallelic record (n_alt == 1) AND ALT is V1 → "new bubble"
  V1b : multi-allelic record (n_alt >  1) AND THIS ALT is V1 → "new allele at existing bubble"

Then if PG AC table is available, compute fraction of V1/V1a/V1b records where ≥1 of
the 151 PG founders carry the ALT (ac_pg > 0). This is the headline number.

Output:
  results/tier1_summary_{chrom}.tsv
  results/tier1_summary_combined.tsv  (after all chroms)
"""
import gzip
import sys
import pandas as pd
import numpy as np
from pathlib import Path

WORK = Path("/carnegie/nobackup/scratch/tbellagio/hapfire_sv/panel_overlap_135_vs_82")
RES = WORK / "results"

def categorize(df):
    """add columns: cat, vtype.
    cat ∈ {shared, cactus_only, V1a, V1b, empty}
    vtype ∈ {SNP, INS<50, DEL<50, INS>=50, DEL>=50, MNP/other}
    """
    df = df.copy()
    df["delta"] = df["alt_len"] - df["ref_len"]
    cat = np.where(
        (df["ac_82"] > 0) & (df["ac_53"] > 0), "shared",
        np.where(
            (df["ac_82"] > 0) & (df["ac_53"] == 0), "cactus_only",
            np.where(
                (df["ac_82"] == 0) & (df["ac_53"] > 0),
                np.where(df["n_alt"] == 1, "V1a", "V1b"),
                "empty",
            )
        )
    )
    df["cat"] = cat
    # variant type per ALT
    def vt(ref_len, alt_len):
        if ref_len == 1 and alt_len == 1:
            return "SNP"
        d = alt_len - ref_len
        if d == 0:
            return "MNP/other"
        elif d > 0:
            return "INS>=50" if d >= 50 else "INS<50"
        else:
            return "DEL>=50" if -d >= 50 else "DEL<50"
    df["vtype"] = [vt(r, a) for r, a in zip(df["ref_len"], df["alt_len"])]
    return df

def summarize(df, chrom):
    df = categorize(df)
    # per-cat × vtype counts
    cnt = df.groupby(["cat", "vtype"]).size().unstack(fill_value=0)
    cnt["TOTAL"] = cnt.sum(axis=1)
    cnt.loc["TOTAL"] = cnt.sum(axis=0)
    cnt["chrom"] = chrom
    return cnt, df

def main(chrom):
    f = RES / f"ac_split_{chrom}.tsv.gz"
    if not f.exists():
        print(f"missing {f}", file=sys.stderr)
        return None
    print(f"[load] {f}")
    df = pd.read_csv(f, sep="\t", dtype={
        "chrom": "string", "pos": "int64",
        "ref_len": "int64", "n_alt": "int32", "alt_idx": "int32",
        "alt_len": "int64",
        "ac_82": "int32", "an_82": "int32",
        "ac_53": "int32", "an_53": "int32"})
    print(f"[load] {len(df):,} per-ALT rows")
    cnt, df2 = summarize(df, chrom)
    out_cnt = RES / f"tier1_summary_{chrom}.tsv"
    cnt.to_csv(out_cnt, sep="\t")
    print(f"[out] {out_cnt}")
    print(cnt.to_string())

    # save the V1 set for join with PG (smaller file)
    v1 = df2[df2["cat"].isin(["V1a", "V1b"])][
        ["chrom", "pos", "n_alt", "alt_idx", "ref_len", "alt_len", "ac_53", "an_53", "cat", "vtype"]
    ]
    out_v1 = RES / f"v1_records_{chrom}.tsv.gz"
    v1.to_csv(out_v1, sep="\t", index=False, compression="gzip")
    print(f"[out] {len(v1):,} V1 records → {out_v1}")

    # Try PG join if available
    pg_file = RES / f"ac_pg151_{chrom}.tsv.gz"
    if pg_file.exists():
        print(f"[pg] joining {pg_file}")
        pg = pd.read_csv(pg_file, sep="\t", dtype={
            "chrom": "string", "pos": "int64",
            "ref_len": "int64", "n_alt": "int32", "alt_idx": "int32",
            "alt_len": "int64",
            "ac_pg": "int32", "an_pg": "int32"})
        m = df2.merge(pg[["chrom", "pos", "alt_idx", "ac_pg", "an_pg"]],
                     on=["chrom", "pos", "alt_idx"], how="left")
        # Was the (chrom, pos, alt_idx) found in PG?
        m["pg_present"] = ~m["ac_pg"].isna()
        m["ac_pg"] = m["ac_pg"].fillna(0).astype(int)
        m["an_pg"] = m["an_pg"].fillna(0).astype(int)
        m["pg_carrier"] = m["ac_pg"] > 0

        # Headline by cat × vtype
        head = m.groupby(["cat", "vtype"]).agg(
            n_records=("pos", "size"),
            pg_present=("pg_present", "sum"),
            n_pg_carrier=("pg_carrier", "sum"),
            mean_ac_pg=("ac_pg", "mean"),
        )
        head["frac_pg_carrier"] = head["n_pg_carrier"] / head["n_records"]
        out_head = RES / f"tier1_headline_{chrom}.tsv"
        head.to_csv(out_head, sep="\t")
        print(f"[out] {out_head}")
        print(head.to_string())
    else:
        print(f"[pg] PG AC table not yet available: {pg_file}")

if __name__ == "__main__":
    for c in (sys.argv[1:] if len(sys.argv) > 1 else ["Chr1", "Chr2", "Chr3", "Chr4", "Chr5"]):
        try:
            main(c)
            print()
        except Exception as e:
            print(f"FAIL {c}: {e}", file=sys.stderr)
