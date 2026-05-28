#!/usr/bin/env python3
"""Isolate the M-step normalizer. Same filt2 matrix; factorial over
weight in {uniform omega=1, per-bubble omega=1/m_b} x normalizer in
{global total, per-founder denom_f=sum_k omega_k cn_fk}. 4 conditions x regimes.
Reuses existing filt2 counts (no jellyfish)."""
import numpy as np, csv, json, math, glob, os
from pathlib import Path
from scipy.sparse import load_npz

ROOT=str(Path(__file__).resolve().parents[3])
CAC=set(map(str,json.load(open(f"{ROOT}/data/founder_split_cactus_pg.json"))["cactus"]))
truth={}
for d in csv.DictReader(open(f"{ROOT}/scratch/g0_sweep_per_founder.tsv"),delimiter="\t"):
    truth.setdefault(d["sim"],{})[d["founder"]]=float(d["truth"])
hf_reps=[]
for s in range(1,9):
    p=f"/global/scratch/projects/fc_moilab/projects/grenenet-phase1/frequency/hapFIRE_frequencies/seed_mix/s{s}_ecotype_frequency.txt"
    if os.path.exists(p):
        hd={}
        for ln in open(p):
            q=ln.split()
            if len(q)>=2:
                try: hd[q[0]]=float(q[1])
                except: pass
        hf_reps.append(hd)

import argparse
ap=argparse.ArgumentParser(); ap.add_argument("--base",default="filt2",choices=["filt2","subsamp"])
args=ap.parse_args()
if args.base=="subsamp":
    CNDIR="cn_full_231_v3qc_v3_subsampMedian_refilt2"; CTAG="subsamp"; OUTSUB="normfact_subsamp"
else:
    CNDIR="cn_full_231_v3qc_v3_filt2"; CTAG="filt2"; OUTSUB="normfact"

print(f"loading {args.base} cn ...",flush=True)
cn=load_npz(f"{ROOT}/poolfreq/data/{CNDIR}/cn_Chr1.cn.npz").tocsr()
meta=np.load(f"{ROOT}/poolfreq/data/{CNDIR}/cn_Chr1.meta.npz",allow_pickle=True)
fo_all=np.asarray(meta["founders"]).astype(str)
bid=np.asarray(meta["bubble_id"]).astype(np.int64)
m_b_full=np.bincount(bid)[bid].astype(np.float64)
F,K=cn.shape
print(f"  F={F} K={K:,}",flush=True)

def em(cn_nz, counts_nz, omega, normalizer, max_iter=300, tol=1e-7):
    F=cn_nz.shape[0]; h=np.full(F,1.0/F,np.float32)
    wc=(omega*counts_nz).astype(np.float32); total=float(wc.sum())
    denom_f=np.maximum(cn_nz@omega,1e-12).astype(np.float32) if normalizer=="perfounder" else None
    for it in range(max_iter):
        mu=np.maximum(h@cn_nz,1e-7); cw=wc/mu; num=h*(cn_nz@cw)
        hn = num/denom_f if normalizer=="perfounder" else num/total
        s=hn.sum()
        if s<=0 or not np.isfinite(s): break
        hn/=s
        if np.linalg.norm(hn-h)<tol: h=hn; break
        h=hn
    return h.astype(np.float64)

def metrics(fo,h,t):
    cac=np.array([x in CAC for x in fo]); ae=np.abs(h-t)
    return float(h[cac].sum()), math.sqrt(np.mean(ae**2)), float(h[t==0].sum())

OUT=f"{ROOT}/scratch/h_fixes/{OUTSUB}"; os.makedirs(OUT,exist_ok=True)
TAG={("uniform","global"):"u_global",("uniform","perfounder"):"u_perf",
     ("1/m_b","global"):"mb_global",("1/m_b","perfounder"):"mb_perf"}
SIMS=["g0_n231_rep0_rand","g0_n50_rep0_cact","g0_n50_rep2_pg"]   # simulated reads only
def counts_path(sim,rep=None):
    if sim=="SEEDMIX": return f"{ROOT}/scratch/seedmix_h_test/{CTAG}_S{rep}.counts.npy"
    return f"{ROOT}/scratch/g0_sweep_h_test/{CTAG}_{sim}.counts.npy"

print(f"\n{'sim':<18}{'weight':<10}{'normalizer':<12}{'cactus':>8}{'RMSE':>9}{'spurious':>9}")
print("-"*68)
for sim in SIMS:
    # densify nz once per sim using the union of nonzero across reps (SEEDMIX) or the single counts
    if sim=="SEEDMIX":
        cset=[np.load(counts_path(sim,r)) for r in range(1,9)]
        nz=np.zeros(K,bool)
        for c in cset: nz|= c>0
    else:
        c0=np.load(counts_path(sim)); nz=c0>0
    cn_nz=np.asarray(cn[:,nz].todense(),dtype=np.float32)
    mb_nz=m_b_full[nz]
    weights={"uniform":np.ones(nz.sum(),np.float32),"1/m_b":(1.0/mb_nz).astype(np.float32)}
    # truth/ref aligned to fo_all
    if sim=="SEEDMIX":
        ref=np.array([np.mean([r.get(x,0) for r in hf_reps]) for x in fo_all])
        t=ref
    else:
        t=np.array([truth[sim].get(x,np.nan) for x in fo_all]); ok=~np.isnan(t)
        # all founders should be present; assume ok all true
    for wname,om in weights.items():
        for norm in ["global","perfounder"]:
            tag=TAG[(wname,norm)]
            if sim=="SEEDMIX":
                hs=[]
                for r in range(1,9):
                    cc=np.load(counts_path(sim,r))[nz].astype(np.float32)
                    hr=em(cn_nz,cc,om,norm); hs.append(hr)
                    np.savez(f"{OUT}/{tag}_S{r}.npz",founders=fo_all,h=hr)
                h=np.mean(hs,axis=0)
            else:
                cc=c0[nz].astype(np.float32); h=em(cn_nz,cc,om,norm)
                np.savez(f"{OUT}/{tag}_{sim}.npz",founders=fo_all,h=h)
            cm,rmse,spur=metrics(fo_all,h,t)
            print(f"{sim:<18}{wname:<10}{norm:<12}{cm:>8.3f}{rmse:>9.5f}{spur:>9.3f}",flush=True)
    del cn_nz
print("\n(baseline filt2 = uniform/global. cactus truth: n231=0.346, n50_cact=0.80, n50_pg=0.10, SEEDMIX~0.33)")
