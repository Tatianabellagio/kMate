#!/usr/bin/env python3
"""Unified scorer across all h-fix methods + baselines. Handles each method's npz format."""
import csv, json, math, glob, os
import numpy as np

ROOT="/carnegie/nobackup/scratch/tbellagio/hapfire_sv"
HF=f"{ROOT}/scratch/h_fixes"
split=json.load(open(f"{ROOT}/data/founder_split_cactus_pg.json")); CAC=set(map(str,split["cactus"]))
G0=["g0_n231_rep0_rand","g0_n50_rep0_cact","g0_n50_rep1_bal","g0_n50_rep2_pg","g0_n200_rep0_rand"]

# g0 truth
truth={}
for d in csv.DictReader(open(f"{ROOT}/scratch/g0_sweep_per_founder.tsv"),delimiter="\t"):
    truth.setdefault(d["sim"],{})[d["founder"]]=float(d["truth"])
# hapFIRE 8-rep avg + uniform
hf_reps=[]
for s in range(1,9):
    p=f"/carnegie/nobackup/scratch/xwu/GrENE_net/hapFIRE_frequencies/seed_mix/s{s}_ecotype_frequency.txt"
    if os.path.exists(p):
        hd={}
        for ln in open(p):
            q=ln.split()
            if len(q)>=2:
                try: hd[q[0]]=float(q[1])
                except: pass
        hf_reps.append(hd)

def hvec(d):
    f=np.asarray(d["founders"]).astype(str)
    if "h" in d: return f, np.asarray(d["h"]).astype(float)
    return f, None  # sweep handled separately

def g0_metrics(founders,h,sim):
    t=np.array([truth[sim].get(x,np.nan) for x in founders]); ok=~np.isnan(t)
    h=h[ok]; t=t[ok]; fo=founders[ok]
    cm=float(h[np.array([x in CAC for x in fo])].sum())
    ae=np.abs(h-t); leak=float(h[t==0].sum())
    return cm, ae.mean(), math.sqrt((ae**2).mean()), ae.max(), leak

def sm_avg(files, gi=None, key=None):
    reps=[]
    for p in files:
        d=np.load(p,allow_pickle=True); f=np.asarray(d["founders"]).astype(str)
        h = np.asarray(d[key])[gi].astype(float) if key else np.asarray(d["h"]).astype(float)
        reps.append(dict(zip(f,h)))
    fs=sorted(reps[0]); avg={x:np.mean([r.get(x,np.nan) for r in reps]) for x in fs}
    unif=1.0/len(fs)
    cm=sum(avg[x] for x in fs if x in CAC)
    ru=math.sqrt(np.mean([(avg[x]-unif)**2 for x in fs]))
    rh=float('nan')
    if hf_reps:
        hfa={x:np.mean([r.get(x,0) for r in hf_reps]) for x in fs}
        rh=math.sqrt(np.mean([(avg[x]-hfa[x])**2 for x in fs]))
    return cm,ru,rh

# ---- method resolvers: (g0 file, seedmix files, h-extractor) ----
def files_g0(method, sim):
    return {"perbubble":f"{HF}/perbubble/g0_{sim}.npz",
            "combined":f"{HF}/combined/{sim}.h_sweep.npz",
            "combined_neg":f"{HF}/combined_neg/{sim}.h_sweep.npz",
            "perbub_on_subsamp":f"{HF}/perbubble_on_subsamp/{sim}.h_sweep.npz",
            "neg_subsamp":f"{HF}/combined_neg_subsamp/{sim}.h_sweep.npz",
            "correlation":f"{HF}/correlation/{sim}.corr.npz",
            "balanced":f"{HF}/balanced/balanced_{sim}.npz",
            "whitening":f"{HF}/whitening/{sim}.whiten.npz",
            "invac":f"{HF}/invac/{sim}.invac.npz"}[method]
def files_sm(method):
    pat={"perbubble":f"{HF}/perbubble/seedmix_S*.npz",
         "combined":f"{HF}/combined/S*.h_sweep.npz",
         "combined_neg":f"{HF}/combined_neg/S*.h_sweep.npz",
         "perbub_on_subsamp":f"{HF}/perbubble_on_subsamp/S*.h_sweep.npz",
         "neg_subsamp":f"{HF}/combined_neg_subsamp/S*.h_sweep.npz",
         "correlation":f"{HF}/correlation/seedmix_S*.corr.npz",
         "balanced":f"{HF}/balanced/balanced_S*.npz",
         "whitening":f"{HF}/whitening/*S*.whiten.npz",
         "invac":f"{HF}/invac/*S*.invac.npz"}[method]
    return sorted(glob.glob(pat))

SWEEP={"combined":("h_per_gamma","gammas"),"combined_neg":("h_per_gamma","gammas"),
       "perbub_on_subsamp":("h_per_gamma","gammas"),"neg_subsamp":("h_per_gamma","gammas"),
       "correlation":("h_per_rho","rhos")}

def run_method(method):
    f0=files_g0(method,"g0_n231_rep0_rand")
    if not os.path.exists(f0): return None
    if method in SWEEP:
        key,gname=SWEEP[method]
        d=np.load(f0,allow_pickle=True); grid=np.asarray(d[gname])
        print(f"\n### {method}: sweep {gname}={list(grid)} ###")
        for gi,g in enumerate(grid):
            fo=np.asarray(d["founders"]).astype(str); h=np.asarray(d[key])[gi]
            cm,mae,rmse,mx,leak=g0_metrics(fo,h,"g0_n231_rep0_rand")
            smf=files_sm(method); scm,sru,srh=(sm_avg(smf,gi,key) if len(smf)>=8 else (float('nan'),)*3)
            print(f"  {gname}={g:<5} g0n231: cactus={cm:.3f} RMSE={rmse:.5f} | SEEDMIX: cactus={scm:.3f} RMSE_unif={sru:.5f}")
    else:
        d=np.load(f0,allow_pickle=True); fo,h=hvec(d)
        cm,mae,rmse,mx,leak=g0_metrics(fo,h,"g0_n231_rep0_rand")
        smf=files_sm(method); scm,sru,srh=(sm_avg(smf) if len(smf)>=8 else (float('nan'),)*3)
        print(f"\n### {method}: g0n231 cactus={cm:.3f} RMSE={rmse:.5f} leak={leak:.3f} | SEEDMIX(n={len(smf)}) cactus={scm:.3f} RMSE_unif={sru:.5f} RMSE_hf={srh:.5f}")

print("BASELINES: filt2 g0n231 cactus=0.503 RMSE=0.00283 | SEEDMIX cactus=0.518 RMSE_unif=0.00318")
print("           subsampMedian g0n231 cactus=0.353 RMSE=0.00225 | SEEDMIX cactus=0.370 RMSE_unif=0.00255")
print("           protect1 g0n231 cactus=0.415 RMSE=0.00220")
print("           TARGET: g0n231 cactus=0.346, SEEDMIX cactus~0.35")
for m in ["balanced","perbubble","correlation","combined","combined_neg","perbub_on_subsamp","neg_subsamp","invac","whitening"]:
    try: run_method(m)
    except Exception as e: print(f"\n### {m}: ERROR {e}")
