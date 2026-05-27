"""For each ecotype in the manifest, check whether its preprocess outputs are
on disk. Emit (panel, ecotype, layout, status, files) per row.

Layout=PE if either RUN's _1/_2 raw files exist or RAW format implies PE
(xwu_BAM lane-concat'd). Layout=SE if a single fastq is the only raw input.
"""
import sys, os, argparse
from pathlib import Path

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--panel", required=True)
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--prep-dir", required=True)
    ap.add_argument("--raw-dir", required=True)
    ap.add_argument("--prep-style", choices=["main","loo"], required=True,
                    help="filename convention: main = `${eco}_1.dedup` / `${eco}.dedup`; loo = `${eco}_1P_dedup` / `${eco}_dedup`")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    prep = Path(args.prep_dir)
    raw  = Path(args.raw_dir)

    rows = []
    with open(args.manifest) as f:
        header = f.readline()
        cols = header.rstrip("\n").split("\t")
        for line in f:
            fields = line.rstrip("\n").split("\t")
            row = dict(zip(cols, fields))
            # Different manifest schemas: main has source_type/ecotype/run, LOO is ecotype/run/...
            if "ecotype_id" in row:
                eco = row["ecotype_id"]
            elif "ecotype" in row:
                eco = row["ecotype"]
            else:
                eco = fields[0] if args.prep_style=="loo" else fields[1]
            run = row.get("run") or row.get("run_accession") or (fields[1] if args.prep_style=="loo" else fields[2])
            src = row.get("source_type") or "ENA"

            # Determine layout from raw files on disk
            eco_raw = raw / eco
            if (eco_raw / f"{run}_1.fastq.gz").exists() and (eco_raw / f"{run}_2.fastq.gz").exists():
                layout = "PE"
            elif (eco_raw / f"{run}.fastq.gz").exists():
                layout = "SE"
            elif src == "xwu_BAM" and (eco_raw / f"{eco}_1.fastq.gz").exists():
                layout = "PE"
            else:
                layout = "MISSING_RAW"

            # Expected output filenames per style+layout
            if args.prep_style == "main":
                pe_r1 = prep / f"{eco}_1.dedup.fq.gz"
                pe_r2 = prep / f"{eco}_2.dedup.fq.gz"
                se    = prep / f"{eco}.dedup.fq.gz"
            else:  # loo
                pe_r1 = prep / f"{eco}_1P_dedup.fq.gz"
                pe_r2 = prep / f"{eco}_2P_dedup.fq.gz"
                se    = prep / f"{eco}_dedup.fq.gz"

            if layout == "PE":
                ok = pe_r1.exists() and pe_r2.exists()
                files = f"{pe_r1.name};{pe_r2.name}"
            elif layout == "SE":
                ok = se.exists()
                files = se.name
            else:
                ok = False
                files = "NA"

            status = "OK" if ok else "MISSING_OUTPUT" if layout != "MISSING_RAW" else "MISSING_RAW"
            rows.append((args.panel, eco, run, src, layout, status, files))

    with open(args.out, "w") as f:
        f.write("panel\tecotype\trun\tsource\tlayout\tstatus\tfiles\n")
        for r in rows:
            f.write("\t".join(map(str, r)) + "\n")

    from collections import Counter
    c = Counter(r[5] for r in rows)
    print(f"[{args.panel}] {len(rows)} ecotypes; status counts: {dict(c)}")
    if c.get("MISSING_OUTPUT",0) or c.get("MISSING_RAW",0):
        print(f"  MISSING:")
        for r in rows:
            if r[5] != "OK":
                print("   ", r)

if __name__ == "__main__":
    main()
