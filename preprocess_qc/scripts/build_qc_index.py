"""
build_qc_index.py — emit a per-ecotype index of preprocessed outputs for the
SLURM array QC job. Reads presence_<panel>.tsv (status=OK rows only) and
writes qc_index_<panel>.tsv with one row per ecotype:

    panel  ecotype  layout  file_r1  file_r2

For SE rows, file_r2 is empty.
"""
import argparse
from pathlib import Path

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--presence", required=True)
    ap.add_argument("--prep-dir", required=True)
    ap.add_argument("--prep-style", choices=["main","loo"], required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    prep = Path(args.prep_dir).resolve()
    rows = []
    with open(args.presence) as f:
        header = f.readline().rstrip("\n").split("\t")
        for line in f:
            r = dict(zip(header, line.rstrip("\n").split("\t")))
            if r["status"] != "OK":
                continue
            eco = r["ecotype"]; layout = r["layout"]
            if args.prep_style == "main":
                f1 = prep / (f"{eco}_1.dedup.fq.gz" if layout == "PE" else f"{eco}.dedup.fq.gz")
                f2 = prep / (f"{eco}_2.dedup.fq.gz" if layout == "PE" else "")
            else:
                f1 = prep / (f"{eco}_1P_dedup.fq.gz" if layout == "PE" else f"{eco}_dedup.fq.gz")
                f2 = prep / (f"{eco}_2P_dedup.fq.gz" if layout == "PE" else "")
            rows.append((r["panel"], eco, layout, str(f1), str(f2) if layout == "PE" else ""))

    with open(args.out, "w") as f:
        f.write("panel\tecotype\tlayout\tfile_r1\tfile_r2\n")
        for r in rows:
            f.write("\t".join(r) + "\n")
    print(f"wrote {len(rows)} rows -> {args.out}")

if __name__ == "__main__":
    main()
