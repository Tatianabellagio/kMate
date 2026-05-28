#!/usr/bin/env python3
"""
Score any set of h_sweep.npz outputs in scratch/g0_sweep_h_test/ against the g0
per-rep truth (taken from g0_sweep_per_founder.tsv, which records founder->truth
per sim; truth is method-independent).

Reports per (sim, method): MAE, RMSE, max single-founder error, off-target leakage,
and cactus-class estimated vs true mass. RMSE + maxErr are the metrics that punish
"a few ecotypes badly wrong" — the failure mode that hurts downstream AF.

Usage: score_g0_sweep.py [tag1 tag2 ...]   (default: filt2 subsamp protect1 protect2)
"""
import csv, math, sys, glob, os
from collections import defaultdict
from pathlib import Path
import numpy as np

ROOT = str(Path(__file__).resolve().parents[2])
OUT_DIR = f"{ROOT}/scratch/g0_sweep_h_test"
PF = f"{ROOT}/scratch/g0_sweep_per_founder.tsv"

tags = sys.argv[1:] or ["filt2", "subsamp", "protect1", "protect2"]

# truth + side per (sim, founder)
truth = defaultdict(dict); side = defaultdict(dict)
with open(PF) as f:
    for d in csv.DictReader(f, delimiter="\t"):
        truth[d["sim"]][d["founder"]] = float(d["truth"])
        side[d["sim"]][d["founder"]] = d["side"]
sims = sorted(truth)


def score(sim, tag):
    p = f"{OUT_DIR}/{tag}_{sim}.h_sweep.npz"
    if not os.path.exists(p):
        return None
    d = np.load(p, allow_pickle=True)
    founders = np.asarray(d["founders"]).astype(str)
    h = np.asarray(d["h_per_alpha"])[0].astype(float)
    t = np.array([truth[sim].get(f, np.nan) for f in founders])
    ok = ~np.isnan(t)
    h, t, founders = h[ok], t[ok], founders[ok]
    ae = np.abs(h - t)
    mae = ae.mean(); rmse = math.sqrt((ae**2).mean()); mx = ae.max()
    leak = float(h[t == 0].sum())
    cmask = np.array([side[sim].get(f) == "cactus" for f in founders])
    return mae, rmse, mx, leak, float(h[cmask].sum()), float(t[cmask].sum())


print(f"{'sim':<20}{'method':<10}{'MAE':>9}{'RMSE':>9}{'maxErr':>9}{'leak':>8}{'C_est':>7}{'C_tru':>7}")
print("-" * 79)
for sim in sims:
    base = None
    for tag in tags:
        r = score(sim, tag)
        if r is None:
            continue
        mae, rmse, mx, leak, ce, ct = r
        flag = ""
        if base is not None:
            flag = f"  RMSE {(rmse/base-1)*100:+.0f}% vs {tags[0]}"
        else:
            base = rmse
        print(f"{sim:<20}{tag:<10}{mae:>9.5f}{rmse:>9.5f}{mx:>9.4f}{leak:>8.3f}{ce:>7.3f}{ct:>7.3f}{flag}")
    print()
