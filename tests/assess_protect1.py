#!/usr/bin/env python3
"""
Assessment of the protect1 cn_full strategy for h accuracy.
Produces (1) per-ecotype k-mer count redistribution plots (filt2 vs protect1 vs
subsampMedian), (2) h-accuracy plots on g0 sims and real SEEDMIX (8-rep average
toward ~uniform), and prints a SEEDMIX score table.

Outputs PNGs to notebooks/plots/protect1_assessment/.
"""
import csv, json, math, glob, os
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = str(Path(__file__).resolve().parents[2])
PLOTS = f"{ROOT}/notebooks/plots/protect1_assessment"; os.makedirs(PLOTS, exist_ok=True)
split = json.load(open(f"{ROOT}/data/founder_split_cactus_pg.json"))
CAC = set(map(str, split["cactus"]))

def load_tags(v):
    """Return per-founder counts in AC strata: ac2, ac3-4, ac5-10, ac>=11, total."""
    d = {}
    with open(f"{ROOT}/scratch/diag/per_founder_tags_{v}.tsv") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            tot=int(r["n_total"]); a2=int(r["n_ac2"]); a4=int(r["n_ac_le4"]); a10=int(r["n_ac_le10"])
            d[r["founder"]] = {"ac2":a2, "ac3_4":a4-a2, "ac5_10":a10-a4, "ac11+":tot-a10, "total":tot}
    return d

def load_h(p):
    d = np.load(p, allow_pickle=True)
    f = np.asarray(d["founders"]).astype(str)
    h = np.asarray(d["h_per_alpha"])[0].astype(float)
    return dict(zip(f, h))

METHODS = [("filt2","filt2"), ("subsampMedian_refilt2","subsampMedian"),
           ("subsampProtect1_refilt2","protect1"), ("classmatchPG_refilt2","classmatchPG")]

# ---------- Plot 1: per-ecotype k-mer count by AC STRATUM, method x class ----------
# (total count is already class-balanced; the imbalance + redistribution is in the
#  AC composition — cactus rich in low-AC discriminating tags.)
tags = {lab: load_tags(v) for v, lab in METHODS}
founders = sorted(tags["filt2"])
STRATA = ["ac2","ac3_4","ac5_10","ac11+"]
fig, axes = plt.subplots(2, 4, figsize=(20, 9), sharex="col")
colors = {"filt2":"#888", "subsampMedian":"#d62728", "protect1":"#2ca02c", "classmatchPG":"#9467bd"}
for ci, cls in enumerate(["cactus", "PG"]):
    fs = [f for f in founders if (f in CAC) == (cls=="cactus")]
    for si, st in enumerate(STRATA):
        ax = axes[ci, si]
        for _, lab in METHODS:
            vals = np.array([tags[lab][f][st] for f in fs])
            ax.hist(vals, bins=30, histtype="step", lw=2, color=colors[lab],
                    label=f"{lab} (med {np.median(vals):.0f})")
        ax.set_title(f"{cls} | {st}"); ax.legend(fontsize=8)
        if si==0: ax.set_ylabel(f"{cls}\n# founders")
        if ci==1: ax.set_xlabel(f"k-mers/founder ({st})")
fig.suptitle("Per-ecotype k-mer count by AC stratum — filt2 vs subsampMedian vs protect1\n"
             "(protect1 keeps ac2 = discriminating tags; subsamples ac>=3)")
fig.tight_layout(); fig.savefig(f"{PLOTS}/1_per_ecotype_kmer_by_AC.png", dpi=110); plt.close(fig)

# bar: cactus vs PG median per stratum per method
fig, axes = plt.subplots(1, len(STRATA), figsize=(20,5))
x = np.arange(len(METHODS)); w=0.35
for si, st in enumerate(STRATA):
    ax=axes[si]
    for i, cls in enumerate(["cactus","PG"]):
        fs=[f for f in founders if (f in CAC)==(cls=="cactus")]
        meds=[np.median([tags[lab][f][st] for f in fs]) for _,lab in METHODS]
        ax.bar(x+(i-0.5)*w, meds, w, label=cls)
    ax.set_xticks(x); ax.set_xticklabels([l for _,l in METHODS], rotation=20)
    ax.set_title(st);
    if si==0: ax.set_ylabel("median k-mers/founder"); ax.legend()
fig.suptitle("Median per-ecotype k-mers by AC stratum, cactus vs PG")
fig.tight_layout(); fig.savefig(f"{PLOTS}/2_median_kmer_by_AC_class.png", dpi=110); plt.close(fig)

# ---------- SEEDMIX: 8-rep average h vs uniform & vs hapFIRE ----------
# hapFIRE 8-rep avg
hf_reps=[]
for s in range(1,9):
    p=f"/global/scratch/projects/fc_moilab/projects/grenenet-phase1/frequency/hapFIRE_frequencies/seed_mix/s{s}_ecotype_frequency.txt"
    if os.path.exists(p):
        hf={}
        for line in open(p):
            q=line.split()
            if len(q)>=2:
                try: hf[q[0]]=float(q[1])
                except: pass
        hf_reps.append(hf)
hf_founders=sorted(hf_reps[0]) if hf_reps else []
hf_avg={f:np.mean([r.get(f,0) for r in hf_reps]) for f in hf_founders} if hf_reps else {}

sm_tags=["filt2","subsamp","f2subsamp","protect1","classmatch"]
sm_avg={}
for tag in sm_tags:
    reps=[load_h(p) for p in sorted(glob.glob(f"{ROOT}/scratch/seedmix_h_test/{tag}_S*.h_sweep.npz"))]
    if not reps: continue
    fs=sorted(reps[0])
    sm_avg[tag]={f:np.mean([r.get(f,np.nan) for r in reps]) for f in fs}
    # per-rep stability
    stab=np.mean([np.std([r[f] for r in reps]) for f in fs])
    n=len(fs); unif=1.0/n
    def rmse(a,b,keys): v=[a[k]-b[k] for k in keys if k in a and k in b]; return math.sqrt(np.mean(np.square(v)))
    def mae(a,b,keys): v=[abs(a[k]-b[k]) for k in keys if k in a and k in b]; return np.mean(v)
    uni={f:unif for f in fs}
    cmass=sum(sm_avg[tag][f] for f in fs if f in CAC)
    line=f"  {tag:<9} nreps={len(reps)} cactus_mass={cmass:.3f} per-rep_std={stab:.5f}"
    line+=f"  RMSE_vs_unif={rmse(sm_avg[tag],uni,fs):.5f}"
    if hf_avg: line+=f"  RMSE_vs_hapFIRE={rmse(sm_avg[tag],hf_avg,fs):.5f}"
    print(line)

# SEEDMIX scatter: protect1 avg-h vs uniform, colored by class
if "protect1" in sm_avg:
    fig, ax = plt.subplots(figsize=(7,7))
    fs=sorted(sm_avg["protect1"]); unif=1.0/len(fs)
    h=np.array([sm_avg["protect1"][f] for f in fs]); cl=np.array([f in CAC for f in fs])
    ax.scatter(np.full(cl.sum(),unif), h[cl], s=14, alpha=.6, label="cactus", color="#1f77b4")
    ax.scatter(np.full((~cl).sum(),unif), h[~cl], s=14, alpha=.6, label="PG", color="#ff7f0e")
    ax.axhline(unif, ls="--", c="k", lw=1, label="uniform 1/231")
    ax.set_xlabel("uniform 1/231"); ax.set_ylabel("protect1 avg-h (8 SEEDMIX reps)")
    ax.set_title("SEEDMIX 8-rep avg h (protect1) vs uniform"); ax.legend()
    fig.tight_layout(); fig.savefig(f"{PLOTS}/3_seedmix_protect1_vs_uniform.png", dpi=110); plt.close(fig)

print(f"\nPNGs -> {PLOTS}")
