#!/usr/bin/env python
"""Compute/runtime stats for a SLURM array (e.g. the kMate pilot).

Pulls sacct accounting for every array task, aggregates across the main/.batch/
.extern steps, and reports wall-time, CPU-hours, CPU efficiency, peak RSS, real
job span, achieved concurrency/throughput, and a projection to the full cohort.

Usage:
  python scripts/pilot_compute_stats.py <jobid> [--cohort 2168] [--ncpus 8]
"""
import argparse, subprocess, sys
from statistics import median, mean

def dur_to_sec(s):
    """Parse sacct duration 'D-HH:MM:SS[.mmm]' / 'HH:MM:SS' / 'MM:SS.mmm' -> sec."""
    if not s or s in ("", "INVALID"): return None
    days=0
    if "-" in s: d,s=s.split("-",1); days=int(d)
    parts=s.split(":")
    parts=[float(p) for p in parts]
    while len(parts)<3: parts.insert(0,0.0)
    h,m,sec=parts
    return days*86400+h*3600+m*60+sec

def rss_to_mb(s):
    if not s: return None
    u=s[-1]
    try:
        if u=="K": return float(s[:-1])/1024
        if u=="M": return float(s[:-1])
        if u=="G": return float(s[:-1])*1024
        return float(s)/1024/1024  # bytes
    except ValueError:
        return None

def parse_dt(s):
    if not s or s in ("Unknown","None"): return None
    # sacct: 2026-06-01T13:16:49
    import datetime
    try: return datetime.datetime.strptime(s,"%Y-%m-%dT%H:%M:%S")
    except ValueError: return None

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("jobid")
    ap.add_argument("--cohort", type=int, default=2168, help="full set size to project")
    ap.add_argument("--ncpus", type=int, default=8, help="cores/task (for projection)")
    args=ap.parse_args()

    fmt="JobID,State,ElapsedRaw,TotalCPU,NCPUS,MaxRSS,Start,End"
    out=subprocess.run(["sacct","-j",args.jobid,"-P","-n","-o",fmt],
                       capture_output=True,text=True).stdout
    tasks={}
    for line in out.strip().split("\n"):
        if not line: continue
        jid,state,eraw,tcpu,ncpus,maxrss,start,end=line.split("|")
        base=jid.split(".")[0]
        t=tasks.setdefault(base,{"state":None,"eraw":None,"ncpus":None,
                                 "tcpu":None,"maxrss":None,"start":None,"end":None})
        if "." not in jid:   # main task line
            t["state"]=state; t["eraw"]=int(eraw) if eraw.isdigit() else None
            t["ncpus"]=int(ncpus) if ncpus.isdigit() else None
            t["start"]=parse_dt(start); t["end"]=parse_dt(end)
        # TotalCPU / MaxRSS usually live on .batch
        c=dur_to_sec(tcpu);
        if c is not None: t["tcpu"]=max(t["tcpu"] or 0, c)
        r=rss_to_mb(maxrss)
        if r is not None: t["maxrss"]=max(t["maxrss"] or 0, r)

    comp=[t for t in tasks.values() if t["state"]=="COMPLETED"]
    failed=[t for t in tasks.values() if t["state"] and t["state"]!="COMPLETED"]
    walls=[t["eraw"] for t in comp if t["eraw"]]
    cpus =[t["tcpu"] for t in comp if t["tcpu"]]
    rss  =[t["maxrss"] for t in comp if t["maxrss"]]
    ncpus=[t["ncpus"] for t in comp if t["ncpus"]] or [args.ncpus]
    nc=ncpus[0]

    print(f"=== COMPUTE STATS  job {args.jobid} ===")
    print(f"tasks: {len(tasks)} total | COMPLETED {len(comp)} | other/failed {len(failed)}")
    if failed:
        from collections import Counter
        print("  non-complete states:", dict(Counter(t['state'] for t in failed)))
    if not walls:
        print("no completed tasks with timing yet."); return

    tot_cpu_h=sum(cpus)/3600 if cpus else 0
    tot_corealloc_h=sum(w*nc for w in walls)/3600
    eff=100*sum(cpus)/sum(w*nc for w in walls) if cpus else float("nan")

    print(f"\n-- per-sample WALL (s): min {min(walls)}  median {int(median(walls))}  "
          f"mean {int(mean(walls))}  max {max(walls)}")
    print(f"-- per-sample WALL (min): median {median(walls)/60:.1f}  max {max(walls)/60:.1f}")
    print(f"-- cores/task: {nc}")
    print(f"-- CPU time actually used: {tot_cpu_h:.1f} core-h  (sum TotalCPU)")
    print(f"-- core-hours ALLOCATED:   {tot_corealloc_h:.1f} core-h  (Elapsed x NCPUS)")
    print(f"-- CPU efficiency: {eff:.0f}%  (used / allocated)")
    if rss:
        print(f"-- peak RSS/task (GB): median {median(rss)/1024:.1f}  max {max(rss)/1024:.1f}")

    # real job span + achieved concurrency
    starts=[t["start"] for t in tasks.values() if t["start"]]
    ends=[t["end"] for t in comp if t["end"]]
    if starts and ends:
        span=(max(ends)-min(starts)).total_seconds()
        conc=tot_corealloc_h/(span/3600)/nc  # avg tasks running concurrently
        print(f"\n-- REAL job span: {span/60:.1f} min  ({len(comp)} samples)")
        print(f"-- achieved avg concurrency: ~{conc:.0f} samples at once")
        print(f"-- throughput: {len(comp)/(span/3600):.0f} samples/hour")

    # cohort projection
    per_core_h=tot_corealloc_h/len(comp)
    per_cpu_h = tot_cpu_h/len(comp) if cpus else 0
    print(f"\n=== PROJECTION to full cohort N={args.cohort} (global, {nc} cores/task) ===")
    print(f"-- core-hours allocated: ~{per_core_h*args.cohort:.0f} core-h "
          f"(~{per_cpu_h*args.cohort:.0f} core-h actually used)")
    med=median(walls)
    for C in (50,100,200):
        wall_h=args.cohort*med/3600/C
        print(f"   at {C:>3} concurrent samples: ~{wall_h:.1f} h wall")

if __name__=="__main__":
    main()
