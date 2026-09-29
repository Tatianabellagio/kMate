#!/usr/bin/env python3
"""Recover kMate speed/compute for the poolsize_depth sweep WITHOUT rerunning
anything: every 07j (p80)/07h (p231) job's stdout already logs a per-stage
wall-time breakdown (k-mer count / EM / total), and SLURM's own `sacct`
accounting (retained since 2025-06) has Elapsed/MaxRSS/TotalCPU per job ID,
which is embedded in the log filename (`<script>_psd_<jobid>.out`). This joins
the two so kMate's speed/compute can be compared against hapFIRE's without any
new instrumentation.

Run with the `basic` env:
    /global/home/users/tbellg/miniforge3/envs/basic/bin/python \
        benchmarks/poolsize_depth/scripts/recover_kmate_speed.py
"""
import os, re, glob, subprocess
import pandas as pd

ROOT = "/global/scratch/users/tbellg/kmate"
OUT = f"{ROOT}/benchmarks/poolsize_depth/results"
os.makedirs(OUT, exist_ok=True)

LOG_GLOBS = {
    "p80":  [f"{ROOT}/benchmarks/p80/logs/07j_psd_*.out",
             f"{ROOT}/benchmarks/p80/scripts/logs/07j_psd_*.out"],
    "p231": [f"{ROOT}/benchmarks/p231/logs/07h_psd_*.out",
             f"{ROOT}/benchmarks/p231/scripts/logs/07h_psd_*.out"],
}

HEADER_RE = re.compile(r"p(?:80|231) --unit chrom\s+N=(\d+)\s+cov=(\d+)\s+seed=(\d+)")
COUNT_RE = re.compile(r"\[Chr1\]\s+(\d+)s count:")
EM_RE = re.compile(r"EM solved in \d+ iters \[(\d+)s\]")
TOTAL_RE = re.compile(r"TOTAL:\s*(\d+)s")
JOBID_RE = re.compile(r"_(\d+)\.out$")

rows = []
for panel, patterns in LOG_GLOBS.items():
    for pattern in patterns:
        for path in glob.glob(pattern):
            jobid_m = JOBID_RE.search(os.path.basename(path))
            if not jobid_m:
                continue
            jobid = jobid_m.group(1)
            text = open(path).read()
            hm = HEADER_RE.search(text)
            if not hm:
                continue
            n, cov, seed = map(int, hm.groups())
            cm, em, tm = COUNT_RE.search(text), EM_RE.search(text), TOTAL_RE.search(text)
            rows.append(dict(
                panel=panel, N=n, depth=cov, seed=seed, jobid=jobid,
                count_s=int(cm.group(1)) if cm else None,
                em_s=int(em.group(1)) if em else None,
                driver_total_s=int(tm.group(1)) if tm else None,
            ))

df = pd.DataFrame(rows).drop_duplicates(subset=["panel", "N", "depth", "seed"], keep="last")
print(f"parsed {len(df)} kMate poolsize_depth log records")

# batch-query sacct for every job ID at once
jobids = df["jobid"].tolist()
fmt = "JobID,Elapsed,TotalCPU,MaxRSS,State"
out = subprocess.run(
    ["sacct", "-j", ",".join(jobids), "--format", fmt, "--parsable2", "--noheader"],
    capture_output=True, text=True, check=True,
).stdout
sacct_rows = {}
for line in out.splitlines():
    parts = line.split("|")
    if len(parts) < 5:
        continue
    raw_jid, elapsed, cpu, maxrss, state = parts[:5]
    is_batch = raw_jid.endswith(".batch")
    jid = raw_jid.split(".")[0]
    prev = sacct_rows.get(jid)
    # MaxRSS lives on the .batch sub-step, not the parent job row -- prefer it.
    if prev is None or (is_batch and not prev.get("_is_batch")):
        sacct_rows[jid] = dict(elapsed=elapsed, total_cpu=cpu, max_rss=maxrss,
                               state=state, _is_batch=is_batch)


def hms_to_s(hms):
    if not hms or hms in ("", "Unknown"):
        return None
    parts = hms.replace("-", ":").split(":")
    parts = [float(p) for p in parts]
    while len(parts) < 3:
        parts.insert(0, 0)
    d = 0
    if "-" in hms:
        d, rest = hms.split("-")
        d = int(d)
        parts = [float(p) for p in rest.split(":")]
        while len(parts) < 3:
            parts.insert(0, 0)
    h, m, s = parts[-3:]
    return d * 86400 + h * 3600 + m * 60 + s


def rss_to_mb(rss):
    if not rss:
        return None
    rss = rss.strip()
    if rss.endswith("K"):
        return float(rss[:-1]) / 1024
    if rss.endswith("M"):
        return float(rss[:-1])
    if rss.endswith("G"):
        return float(rss[:-1]) * 1024
    try:
        return float(rss) / 1024
    except ValueError:
        return None


df["elapsed_s"] = df["jobid"].map(lambda j: hms_to_s(sacct_rows.get(j, {}).get("elapsed")))
df["cpu_s"] = df["jobid"].map(lambda j: hms_to_s(sacct_rows.get(j, {}).get("total_cpu")))
df["max_rss_mb"] = df["jobid"].map(lambda j: rss_to_mb(sacct_rows.get(j, {}).get("max_rss")))
df["state"] = df["jobid"].map(lambda j: sacct_rows.get(j, {}).get("state"))

table_path = f"{OUT}/kmate_speed_table.tsv"
df.sort_values(["panel", "N", "depth", "seed"]).to_csv(table_path, sep="\t", index=False)
print(f"wrote {table_path} ({len(df)} rows)")
print(df.groupby(["panel", "N", "depth"])[["elapsed_s", "max_rss_mb"]].mean().round(2))
