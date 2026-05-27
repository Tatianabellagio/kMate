"""
parse_dup_rates.py — scrape Clumpify "Duplicates Found: N (P%)" + Reads
In/Out from preprocess SLURM logs and emit per-(panel, ecotype) TSV.

Each preprocess SLURM array task produces both .out (the script's own echos
including the ecotype banner) and .err (Trimmomatic + Clumpify stderr, which
is where Reads In/Out/Duplicates Found land). We merge both halves per
(jobid, array_idx) before emitting.
"""
import re, argparse, glob, os
from pathlib import Path
from collections import defaultdict

DUP_RE = re.compile(r"Duplicates Found:\s+(\d+)\s+([\d\.]+)%")
IN_RE  = re.compile(r"Reads In:\s+(\d+)")
OUT_RE = re.compile(r"Reads Out:\s+(\d+)")
ECO_MAIN_RE = re.compile(r"\] (\d+): source=")
ECO_LOO_RE  = re.compile(r"\] LOO (\d+): R1=")
FNAME_RE = re.compile(r"^(?:loo_)?prep_(\d+)_(\d+)\.(out|err)$")

def parse_log(path):
    """Yield (key, value) pairs found in one log file."""
    found = {}
    try:
        with open(path, errors="replace") as f:
            for line in f:
                if "eco" not in found:
                    m = ECO_MAIN_RE.search(line) or ECO_LOO_RE.search(line)
                    if m: found["eco"] = m.group(1)
                m = IN_RE.search(line)
                if m: found["reads_in"] = int(m.group(1))
                m = OUT_RE.search(line)
                if m: found["reads_out"] = int(m.group(1))
                m = DUP_RE.search(line)
                if m:
                    found["dup_count"] = int(m.group(1))
                    found["dup_pct"]   = float(m.group(2))
    except FileNotFoundError:
        pass
    return found

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--logs-dir", required=True)
    ap.add_argument("--prefix",   required=True, help="e.g. 'prep_' or 'loo_prep_'")
    ap.add_argument("--panel",    required=True)
    ap.add_argument("--out",      required=True)
    args = ap.parse_args()

    # Group .out + .err under the same (jobid, array_idx) key
    grouped = defaultdict(dict)
    for path in glob.glob(f"{args.logs_dir}/{args.prefix}*.[oe]??"):
        name = os.path.basename(path)
        m = FNAME_RE.match(name)
        if not m: continue
        jid, idx, ext = m.group(1), m.group(2), m.group(3)
        key = (int(jid), int(idx))
        grouped[key][ext] = path

    # For each task, parse both halves; later jobs overwrite earlier (so a re-fire
    # of a previously failed array task wins).
    by_eco = {}
    for (jid, idx), parts in sorted(grouped.items()):
        rec = {}
        for ext, p in parts.items():
            rec.update(parse_log(p))
        if "eco" not in rec or "dup_pct" not in rec:
            continue
        log_files = ",".join(f"{args.prefix}{jid}_{idx}.{e}" for e in sorted(parts.keys()))
        by_eco[rec["eco"]] = (rec.get("reads_in"), rec.get("reads_out"),
                              rec.get("dup_count"), rec["dup_pct"], log_files)

    with open(args.out, "w") as f:
        f.write("panel\tecotype\treads_in\treads_out\tdup_count\tdup_pct\tlog\n")
        for eco, (rin, rout, dn, dp, log) in sorted(by_eco.items(), key=lambda kv: int(kv[0])):
            f.write(f"{args.panel}\t{eco}\t{rin}\t{rout}\t{dn}\t{dp}\t{log}\n")
    print(f"wrote {len(by_eco)} rows -> {args.out}")

if __name__ == "__main__":
    main()
