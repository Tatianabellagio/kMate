#!/usr/bin/env python3
"""Ash-style vg-SV evaluation: score vg on ITS OWN native variant set (the 82-acc
deconstruct snarls), not kMate's foreign 135-asm atomic panel. Truth = pool AF
projected from the deconstruct's per-founder GTs; match vg<->truth by snarl traversal
(exact, same node space). Emits an npz of (truth, vg) for plotting."""
import re, gzip, subprocess, sys
import numpy as np, pandas as pd

GDIR = "/global/scratch/users/tbellg/pang/pang_1001gplus/pang/output"
DECON = f"{GDIR}/pang_1001gplus_82acc.vcf.gz"
VGVCF = sys.argv[1]   # vg call -a output
OUT = sys.argv[2]     # npz
POOLW = sys.argv[3] if len(sys.argv) > 3 else \
    "benchmarks/p80/sims/cov10_n80_g0_s42_hotspots_p80_chr1/pool_weights.tsv"
NODES = re.compile(r"\d+")
def tkey(s): return tuple(int(x) for x in NODES.findall(s))

ren = pd.read_csv("data/sv_panel_to_accession_id.tsv", sep="\t", dtype=str)
asm2acc = dict(zip(ren.Assembly_ID, ren.Accession_ID))
pw = pd.read_csv(POOLW, sep="\t", dtype={"founder": str})
w_by_acc = dict(zip(pw.founder, pw.weight.astype(float)))
samples = subprocess.run(["bcftools", "query", "-l", DECON], capture_output=True, text=True).stdout.split()
wvec = np.array([w_by_acc.get(asm2acc.get(s, ""), 0.0) for s in samples])
Sw = wvec.sum()
print(f"founders weighted: {(wvec>0).sum()}/{len(samples)} Sw={Sw:.3f}", flush=True)

truth = {}
with gzip.open(DECON, "rt") as fh:
    for ln in fh:
        if ln[0] == "#":
            continue
        f = ln.rstrip("\n").split("\t")
        ref, alts = f[3], f[4].split(",")
        m = re.search(r"AT=([^;\t]*)", f[7])
        if not m:
            continue
        travs = m.group(1).split(",")
        gts = None
        for i, a in enumerate(alts):
            if abs(len(a) - len(ref)) >= 50 and i + 1 < len(travs):
                if gts is None:
                    gts = [s.partition(":")[0] for s in f[9:]]
                tgt = str(i + 1)
                carr = np.fromiter((1.0 if g == tgt else 0.0 for g in gts), float, len(gts))
                truth[tkey(travs[i + 1])] = float((wvec * carr).sum() / Sw)
print(f"deconstruct SV alt-traversals: {len(truth):,}", flush=True)

vgaf = {}
with open(VGVCF) as fh:
    for ln in fh:
        if ln[0] == "#":
            continue
        f = ln.rstrip("\n").split("\t")
        m = re.search(r"AT=([^;\t]*)", f[7])
        if not m:
            continue
        travs = m.group(1).split(",")
        fd = dict(zip(f[8].split(":"), f[9].split(":")))
        try:
            adv = [int(x) for x in fd.get("AD", "").split(",")]
        except ValueError:
            continue
        tot = sum(adv)
        if tot <= 0 or len(adv) != len(travs):
            continue
        for i in range(1, len(travs)):
            vgaf[tkey(travs[i])] = adv[i] / tot

keys = [k for k in truth if k in vgaf]
t = np.array([truth[k] for k in keys]); v = np.array([vgaf[k] for k in keys])
r2 = 1 - ((v - t) ** 2).sum() / ((t - t.mean()) ** 2).sum()
print(f"MATCHED {len(keys):,}/{len(truth):,} = {len(keys)/len(truth):.1%} | "
      f"vg R2={r2:.4f} MAE={np.abs(v-t).mean():.4f} r={np.corrcoef(v,t)[0,1]:.4f}", flush=True)
np.savez(OUT, truth=t, vg=v)
