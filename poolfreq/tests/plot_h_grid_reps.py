#!/usr/bin/env python3
"""Replicate-AVERAGED h grid (sorted-by-truth style). 4 methods x 3 regimes.
For each (method,regime): average per-founder h and per-founder truth over the 5 seeds,
then plot sorted by avg-truth; gray=truth, blue=cactus est, orange=PG est."""
import numpy as np, json, glob, os, math
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
ROOT="/carnegie/nobackup/scratch/tbellagio/hapfire_sv"; HDIR=f"{ROOT}/scratch/h_fixes/replicates/h"
REPS=f"{ROOT}/sims/visor_freqk/g0_reps"
CAC=set(map(str,json.load(open(f"{ROOT}/data/founder_split_cactus_pg.json"))["cactus"]))
SEEDS=[101,102,103,104,105]
# (label, h-file tag)
M=[("filt2","filt2"),("subsampMedian","subsampMedian"),
   ("fact w=1/m_b global","fact_w=1m_b_global"),("subsamp+1/m_b global","subsamp1m_b_global")]
# (regime label, sim-name template)
R=[("n=231 (all)","g0_n231_rep{s}_rand"),("n=50 cact-heavy","g0_n50_rep{s}_cact"),
   ("n=50 pg-heavy","g0_n50_rep{s}_pg")]

def load_truth(sim):
    th={}
    p=f"{REPS}/{sim}/h_truth.tsv"
    if not os.path.exists(p): return None
    for ln in open(p):
        q=ln.split()
        if q and q[0]!="founder": th[q[0]]=float(q[1])
    return th

def avg_over_seeds(tag,tmpl):
    hs=[]; ths=[]; fo_ref=None
    for s in SEEDS:
        sim=tmpl.format(s=s); f=f"{HDIR}/{tag}__{sim}.npz"
        if not os.path.exists(f): continue
        d=np.load(f,allow_pickle=True); fo=np.asarray(d["founders"]).astype(str)
        fo_ref=fo; hs.append(np.asarray(d["h"]).astype(float))
        th=load_truth(sim); ths.append(np.array([th.get(x,0.0) for x in fo]))
    if not hs: return None
    return fo_ref, np.mean(hs,axis=0), np.mean(ths,axis=0), len(hs)

fig,ax=plt.subplots(len(M),len(R),figsize=(5.2*len(R),3.0*len(M)))
for ri,(mlab,tag) in enumerate(M):
    for ci,(rlab,tmpl) in enumerate(R):
        a=ax[ri,ci]; res=avg_over_seeds(tag,tmpl)
        if res is None: a.text(.5,.5,"n/a",ha="center",transform=a.transAxes); continue
        fo,h,t,nrep=res; cac=np.array([x in CAC for x in fo])
        order=np.argsort(t); h,t,cac=h[order],t[order],cac[order]; x=np.arange(len(fo))
        a.scatter(x,t,s=6,c="0.6",zorder=1,label="truth")
        a.scatter(x[cac],h[cac],s=12,marker="x",c="#1f77b4",zorder=3,label="cactus")
        a.scatter(x[~cac],h[~cac],s=12,marker="x",c="#ff7f0e",zorder=2,label="PG")
        lvl=t[t>0].max() if (t>0).any() else None
        if lvl: a.axhline(lvl,ls="--",lw=.8,c="r")
        rmse=math.sqrt(np.mean((h-t)**2)); spur=float(h[t==0].sum()); ct=t[cac].sum(); ce=h[cac].sum()
        a.set_title(f"{rlab} — {mlab}  (avg of {nrep} reps)\nRMSE={rmse:.5f} spurious={spur:.3f} cact {ct:.2f}->{ce:.2f}",fontsize=8)
        a.set_ylabel("avg h",fontsize=8)
        if ri==0 and ci==0: a.legend(fontsize=7,markerscale=1.5,loc="upper left")
for ci in range(len(R)): ax[-1,ci].set_xlabel("founder (sorted by avg truth h)")
fig.tight_layout(); out=f"{ROOT}/notebook/plots/H_VECTOR_grid_reps.png"
fig.savefig(out,dpi=110,bbox_inches="tight"); print("wrote",out)
