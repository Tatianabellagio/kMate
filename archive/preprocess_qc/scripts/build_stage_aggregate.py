"""
build_stage_aggregate.py — assemble per-ecotype stage-by-stage retention
across both panels by joining: ENA manifest, Trimmomatic, Clumpify, QC counts.

The pipeline is:
  ENA reads (manifest)
    -> raw fastqs on disk (size sanity)
    -> Trimmomatic input  (Input Read Pairs / Input Reads)
    -> Trimmomatic output (Both Surviving for PE, Surviving for SE)
    -> Clumpify input     (Reads In) — for PE this is 2 * Both Surviving
    -> Clumpify output    (Reads Out)
    -> QC counted on disk (n_r1 + n_r2 from zcat)

This script joins everything and emits cross-checks. Anomaly columns:
  ena_to_trim_in_pct   should be ~100 (small drops are paired-only / SE conversion)
  trim_surv_pct        Trimmomatic survival
  dedup_dropped_pct    Clumpify dup pct
  clump_vs_qc_diff     |reads_out - qc_counted| / reads_out  (should be ~0)

Layout-aware handling: PE Trimmomatic output is "Both Surviving" pairs, but
Clumpify Reads In counts pairs as 2 reads.
"""
import argparse, sys
import pandas as pd
import numpy as np

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ena-main",  required=True)
    ap.add_argument("--ena-loo",   required=True)
    ap.add_argument("--trim-main", required=True)
    ap.add_argument("--trim-loo",  required=True)
    ap.add_argument("--dup-main",  required=True)
    ap.add_argument("--dup-loo",   required=True)
    ap.add_argument("--qc",        required=True)
    ap.add_argument("--out",       required=True)
    args = ap.parse_args()

    parts = []
    for panel, ena_path, trim_path, dup_path in [
        ("main", args.ena_main, args.trim_main, args.dup_main),
        ("loo",  args.ena_loo,  args.trim_loo,  args.dup_loo),
    ]:
        ena = pd.read_csv(ena_path, sep="\t").rename(columns={"ecotype_id":"ecotype"})
        ena["ecotype"] = ena["ecotype"].astype(str)
        ena_keep = ena[["ecotype","read_count","base_count","library_layout","instrument_model"]].copy()

        trim = pd.read_csv(trim_path, sep="\t")
        trim["ecotype"] = trim["ecotype"].astype(str)

        dup = pd.read_csv(dup_path, sep="\t")[["ecotype","reads_in","reads_out","dup_count","dup_pct"]]
        dup.columns = ["ecotype","clump_in","clump_out","dup_count","dup_pct"]
        dup["ecotype"] = dup["ecotype"].astype(str)

        df = (trim.merge(ena_keep, on="ecotype", how="left")
                  .merge(dup, on="ecotype", how="left"))
        df["panel"] = panel
        parts.append(df)

    df = pd.concat(parts, ignore_index=True)

    # QC merged in
    qc = pd.read_csv(args.qc, sep="\t")[
        ["panel","ecotype","layout","gz_r1","gz_r2","n_r1","n_r2","pair_ok",
         "bases_total","cov_est","lmean_r1","lmean_r2"]].copy()
    qc["ecotype"] = qc["ecotype"].astype(str)
    qc["qc_reads_total"] = qc.n_r1 + qc.n_r2.fillna(0)
    df = df.merge(qc, on=["panel","ecotype","layout"], how="left")

    # Stage-by-stage interpretable counts. For PE, "reads" usually means read PAIRS
    # at Trimmomatic input/output but PAIRS doubled at Clumpify. Normalize to pairs.
    is_pe = df.layout == "PE"
    df["trim_out"] = np.where(is_pe, df.trim_both_surv, df.trim_surv)
    # Clumpify counts each end of a pair separately. For PE comparison, convert
    # clump_in to "pairs equivalent" = clump_in/2.
    df["clump_in_pairs"]  = np.where(is_pe, df.clump_in/2,  df.clump_in)
    df["clump_out_pairs"] = np.where(is_pe, df.clump_out/2, df.clump_out)
    df["qc_reads_pairs"]  = np.where(is_pe, df.qc_reads_total/2, df.qc_reads_total)

    # Retention percentages (each is post / pre at that stage)
    df["pct_ena_to_trim_in"] = df.trim_in / df.read_count * 100
    df["pct_trim_surv"]      = df.trim_out / df.trim_in * 100
    df["pct_clump_link"]     = df.clump_in_pairs / df.trim_out * 100   # should be ~100; pre/post Clumpify input check
    df["pct_clump_keep"]     = df.clump_out_pairs / df.clump_in_pairs * 100
    df["pct_qc_link"]        = df.qc_reads_pairs / df.clump_out_pairs * 100  # should be ~100; file integrity vs Clumpify Reads Out

    # Bool flags for anomalies
    df["flag_ena_link"]   = (df.pct_ena_to_trim_in.between(95, 105) == False) & df.read_count.notna()
    df["flag_low_trim"]   = df.pct_trim_surv < 50
    df["flag_high_dup"]   = df.dup_pct > 50
    df["flag_clump_link"] = (df.pct_clump_link.between(95, 105) == False) & df.clump_in.notna()
    df["flag_qc_link"]    = (df.pct_qc_link.between(95, 105) == False) & df.qc_reads_pairs.notna()
    df["flag_low_cov"]    = df.cov_est < 5

    df.to_csv(args.out, sep="\t", index=False)
    print(f"wrote {len(df)} rows -> {args.out}")

    # Quick console report
    flags = ["flag_ena_link","flag_low_trim","flag_high_dup","flag_clump_link","flag_qc_link","flag_low_cov"]
    print("\n=== anomaly counts ===")
    print(df[flags].sum().to_string())
    for f in flags:
        bad = df[df[f].fillna(False)]
        if len(bad):
            print(f"\n--- {f} ---")
            cols = ["panel","ecotype","layout","read_count","trim_in","trim_out",
                    "clump_in_pairs","clump_out_pairs","qc_reads_pairs","cov_est",
                    "pct_ena_to_trim_in","pct_trim_surv","pct_clump_link","pct_clump_keep","pct_qc_link","dup_pct"]
            print(bad[cols].to_string(index=False))

if __name__ == "__main__":
    main()
