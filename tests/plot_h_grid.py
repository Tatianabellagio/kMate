#!/usr/bin/env python3
"""H-vector grid, 'sorted-by-truth' style. cols=[g0_n231, g0_n50_cact, SEEDMIX], rows=methods.
x = founder rank sorted by truth h ascending (hapFIRE for SEEDMIX); gray=truth, blue=cactus est,
orange=PG est, dashed=present-founder truth level. Title: regime/method + MAE, spurious_h(absent),
cact_sum truth->est."""
import numpy as np, csv, json, math, glob, os
from pathlib import Path
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

ROOT=str(Path(__file__).resolve().parents[2]); HF=f"{ROOT}/scratch/h_fixes"
G0S="g0_sweep_h_test"; SMD="seedmix_h_test"
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

def ex(d,kind):
    if kind=="h": return np.asarray(d["h"]).astype(float)
    if kind=="alpha0": return np.asarray(d["h_per_alpha"])[0].astype(float)
    if kind.startswith("g:"): return np.asarray(d["h_per_gamma"])[int(kind[2:])].astype(float)
    if kind.startswith("r:"): return np.asarray(d["h_per_rho"])[int(kind[2:])].astype(float)
    if kind.startswith("k:"): return np.asarray(d[kind[2:]]).astype(float)

M=[
 ("raw (no filter, no correction)", lambda s:f"{ROOT}/scratch/{G0S}/raw_{s}.h_sweep.npz", f"{ROOT}/scratch/{SMD}/raw_S*.h_sweep.npz","alpha0"),
 ("filt2 (baseline)", lambda s:f"{ROOT}/scratch/{G0S}/filt2_{s}.h_sweep.npz", f"{ROOT}/scratch/{SMD}/filt2_S*.h_sweep.npz","alpha0"),
 ("subsampMedian", lambda s:f"{ROOT}/scratch/{G0S}/subsamp_{s}.h_sweep.npz", f"{ROOT}/scratch/{SMD}/subsamp_S*.h_sweep.npz","alpha0"),
 ("effective-n rho=0.9", lambda s:f"{HF}/correlation/{s}.corr.npz", f"{HF}/correlation/seedmix_S*.corr.npz","r:2"),
 # --- normalizer factorial: same filt2 matrix, weight x normalizer ---
 ("fact: w=1, global (=filt2)", lambda s:f"{HF}/normfact/u_global_{s}.npz", f"{HF}/normfact/u_global_S*.npz","h"),
 ("fact: w=1, per-founder", lambda s:f"{HF}/normfact/u_perf_{s}.npz", f"{HF}/normfact/u_perf_S*.npz","h"),
 ("fact: w=1/m_b, global", lambda s:f"{HF}/normfact/mb_global_{s}.npz", f"{HF}/normfact/mb_global_S*.npz","h"),
 ("fact: w=1/m_b, per-founder", lambda s:f"{HF}/normfact/mb_perf_{s}.npz", f"{HF}/normfact/mb_perf_S*.npz","h"),
 # --- subsampMedian + 1/m_b on top (the requested final combination) ---
 ("subsamp + 1/m_b, global", lambda s:f"{HF}/normfact_subsamp/mb_global_{s}.npz", f"{HF}/normfact_subsamp/mb_global_S*.npz","h"),
 ("subsamp + 1/m_b, per-founder (=oldC2)", lambda s:f"{HF}/normfact_subsamp/mb_perf_{s}.npz", f"{HF}/normfact_subsamp/mb_perf_S*.npz","h"),
 # --- de-focused (commented per request) ---
 #("protect1", lambda s:f"{ROOT}/scratch/{G0S}/protect1_{s}.h_sweep.npz", f"{ROOT}/scratch/{SMD}/protect1_S*.h_sweep.npz","alpha0"),
 #("per-bubble (1/m_b)", lambda s:f"{HF}/perbubble/g0_{s}.npz", f"{HF}/perbubble/seedmix_S*.npz","h"),
 #("inverse-AC g=0.5", lambda s:f"{HF}/invac/g0_{s}.invac.npz", f"{HF}/invac/seedmix_S*.invac.npz","k:h_gamma0.5"),
 #("down-wt-rare g=-1", lambda s:f"{HF}/combined_neg/{s}.h_sweep.npz", f"{HF}/combined_neg/S*.h_sweep.npz","g:0"),
 #("balanced-SNP", lambda s:f"{HF}/balanced/balanced_{s}.npz", f"{HF}/balanced/balanced_S*.npz","h"),
 #("WINNER perbub-on-subsamp", lambda s:f"{HF}/perbubble_on_subsamp/{s}.h_sweep.npz", f"{HF}/perbubble_on_subsamp/S*.h_sweep.npz","g:0"),
]
COLS=[("g0_n231_rep0_rand","n=231 (all)"),("g0_n50_rep0_cact","n=50, cact-heavy"),
      ("g0_n50_rep2_pg","n=50, pg-heavy")]   # simulated reads only (SEEDMIX dropped)

def get_g0(fn,kind,sim):
    if not os.path.exists(fn): return None
    d=np.load(fn,allow_pickle=True); fo=np.asarray(d["founders"]).astype(str); h=ex(d,kind)
    t=np.array([truth[sim].get(x,np.nan) for x in fo]); ok=~np.isnan(t)
    return fo[ok],h[ok],t[ok]
def get_sm(glob_pat,kind):
    fs=sorted(glob.glob(glob_pat))
    if len(fs)<8 or not hf_reps: return None
    reps=[]
    for p in fs:
        d=np.load(p,allow_pickle=True); reps.append(dict(zip(np.asarray(d["founders"]).astype(str),ex(d,kind))))
    fo=sorted(reps[0]); h=np.array([np.mean([r.get(x,np.nan) for r in reps]) for x in fo])
    ref=np.array([np.mean([r.get(x,0) for r in hf_reps]) for x in fo])  # hapFIRE as truth proxy
    return np.array(fo),h,ref

def panel(ax,fo,h,t,title,present_level):
    cac=np.array([x in CAC for x in fo])
    order=np.argsort(t); h,t,cac=h[order],t[order],cac[order]
    x=np.arange(len(fo))
    ax.scatter(x,t,s=6,c="0.6",label="truth",zorder=1)
    ax.scatter(x[cac],h[cac],s=12,marker="x",c="#1f77b4",label="cactus",zorder=3)
    ax.scatter(x[~cac],h[~cac],s=12,marker="x",c="#ff7f0e",label="PG",zorder=2)
    if present_level: ax.axhline(present_level,ls="--",lw=.8,c="r")
    rmse=math.sqrt(np.mean((h-t)**2)); spur=float(h[t==0].sum()); ct=t[cac].sum(); ce=h[cac].sum()
    ax.set_title(f"{title}\nRMSE={rmse:.5f}  spurious(absent)={spur:.3f}  cact_sum {ct:.2f}->{ce:.2f}",fontsize=8)
    return rmse

fig,ax=plt.subplots(len(M),len(COLS),figsize=(5.2*len(COLS),3.0*len(M)))
for ri,(name,g0f,smg,kind) in enumerate(M):
    for ci,(sim,clab) in enumerate(COLS):
        a=ax[ri,ci]
        if sim=="__SEEDMIX__":
            r=get_sm(smg,kind)
            if r is None: a.text(.5,.5,"n/a",ha="center",transform=a.transAxes); continue
            fo,h,ref=r; panel(a,fo,h,ref,f"{clab} — {name}",1.0/len(fo))
        else:
            r=get_g0(g0f(sim),kind,sim)
            if r is None: a.text(.5,.5,"n/a",ha="center",transform=a.transAxes); continue
            fo,h,t=r; lvl=t[t>0].max() if (t>0).any() else None
            panel(a,fo,h,t,f"{clab} — {name}",lvl)
        if ri==0 and ci==0: a.legend(fontsize=7,markerscale=1.5,loc="upper left")
        a.set_ylabel("h",fontsize=8)
for ci in range(len(COLS)): ax[-1,ci].set_xlabel("founder (sorted by truth h asc)")
fig.tight_layout()
out=f"{ROOT}/notebook/plots/H_VECTOR_grid_methods.png"
fig.savefig(out,dpi=105,bbox_inches="tight"); print("wrote",out)
