"""
parse_trim_stats.py — scrape Trimmomatic Input + Both Surviving / Surviving
from preprocess SLURM .err logs. Merge per (jobid, array_idx) so .out
ecotype banner ties to .err Trimmomatic numbers, matching parse_dup_rates.py.

Trimmomatic emits one of:
  PE: "Input Read Pairs: 10270404 Both Surviving: 69 (0.00%) Forward Only Surviving: 4448 (0.04%) Reverse Only Surviving: 145 (0.00%) Dropped: 10265742 (99.95%)"
  SE: "Input Reads: 26359926 Surviving: 563133 (2.14%) Dropped: 25796793 (97.86%)"

Output rows: panel, ecotype, layout, trim_in, trim_both_surv, trim_fwd_only,
             trim_rev_only, trim_dropped, trim_drop_pct
"""
import re, argparse, glob, os
from collections import defaultdict

PE_RE = re.compile(
    r"Input Read Pairs:\s+(\d+)\s+Both Surviving:\s+(\d+)\s+\(([\d\.]+)%\)\s+"
    r"Forward Only Surviving:\s+(\d+)\s+\([\d\.]+%\)\s+"
    r"Reverse Only Surviving:\s+(\d+)\s+\([\d\.]+%\)\s+"
    r"Dropped:\s+(\d+)\s+\(([\d\.]+)%\)")
SE_RE = re.compile(
    r"Input Reads:\s+(\d+)\s+Surviving:\s+(\d+)\s+\(([\d\.]+)%\)\s+"
    r"Dropped:\s+(\d+)\s+\(([\d\.]+)%\)")
ECO_MAIN_RE = re.compile(r"\] (\d+): source=")
ECO_LOO_RE  = re.compile(r"\] LOO (\d+): R1=")
FNAME_RE = re.compile(r"^(?:loo_)?prep_(\d+)_(\d+)\.(out|err)$")

def parse(path):
    found = {}
    try:
        with open(path, errors="replace") as f:
            for line in f:
                if "eco" not in found:
                    m = ECO_MAIN_RE.search(line) or ECO_LOO_RE.search(line)
                    if m: found["eco"] = m.group(1)
                m = PE_RE.search(line)
                if m:
                    found.update(layout="PE", trim_in=int(m.group(1)),
                                 trim_both_surv=int(m.group(2)), trim_both_pct=float(m.group(3)),
                                 trim_fwd_only=int(m.group(4)), trim_rev_only=int(m.group(5)),
                                 trim_dropped=int(m.group(6)), trim_drop_pct=float(m.group(7)))
                m = SE_RE.search(line)
                if m and "layout" not in found:
                    found.update(layout="SE", trim_in=int(m.group(1)),
                                 trim_surv=int(m.group(2)), trim_surv_pct=float(m.group(3)),
                                 trim_dropped=int(m.group(4)), trim_drop_pct=float(m.group(5)))
    except FileNotFoundError:
        pass
    return found

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--logs-dir", required=True)
    ap.add_argument("--prefix",   required=True)
    ap.add_argument("--panel",    required=True)
    ap.add_argument("--out",      required=True)
    args = ap.parse_args()

    grouped = defaultdict(dict)
    for path in glob.glob(f"{args.logs_dir}/{args.prefix}*.[oe]??"):
        name = os.path.basename(path)
        m = FNAME_RE.match(name)
        if not m: continue
        jid, idx, ext = int(m.group(1)), int(m.group(2)), m.group(3)
        grouped[(jid, idx)][ext] = path

    by_eco = {}
    for (jid, idx), parts in sorted(grouped.items()):
        rec = {}
        for ext, p in parts.items():
            rec.update(parse(p))
        if "eco" not in rec or "trim_in" not in rec:
            continue
        # later (jid, idx) wins (re-fired tasks override)
        by_eco[rec["eco"]] = rec

    cols = ["panel","ecotype","layout","trim_in",
            "trim_both_surv","trim_fwd_only","trim_rev_only",
            "trim_surv",
            "trim_dropped","trim_drop_pct"]
    with open(args.out, "w") as f:
        f.write("\t".join(cols) + "\n")
        for eco, rec in sorted(by_eco.items(), key=lambda kv: int(kv[0])):
            row = [args.panel, eco, rec.get("layout",""),
                   rec.get("trim_in",""),
                   rec.get("trim_both_surv",""),
                   rec.get("trim_fwd_only",""),
                   rec.get("trim_rev_only",""),
                   rec.get("trim_surv",""),
                   rec.get("trim_dropped",""),
                   rec.get("trim_drop_pct","")]
            f.write("\t".join(str(x) for x in row) + "\n")
    print(f"wrote {len(by_eco)} rows -> {args.out}")

if __name__ == "__main__":
    main()
