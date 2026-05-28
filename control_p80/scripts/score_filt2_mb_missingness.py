#!/usr/bin/env python3
"""filt2+1/m_b vs filt2-uniform AF accuracy under 3 missingness conditions.
n_called (of 80 founders) comes from the est TSV (h-independent panel QC, so the
mask is identical across weights -> fair A/B). missing_frac = (80-n_called)/80.
Conditions: full | miss<=0.90 | miss<=0.50. Writes summary TSV with miss_cond."""
import os, math
from pathlib import Path
import numpy as np, pandas as pd
CTRL=str(Path(__file__).resolve().parents[1])
KEYS=["chrom","pos","ref_len","alt_len"]
REGIMES=["n50_g0","n231_g0","n50_g1","n231_g1","n50_g3","n50_g3_dom500"]
F=80
CONDS=[("full",1.01),("miss_le_90",0.90),("miss_le_50",0.50)]
def vclass(r,a): return np.where(np.maximum(r,a)>=50,"SV",np.where((r==1)&(a==1),"SNP","indel"))
def truth_path(reg):
    n=reg.split("_g")[0][1:]; g=reg.split("_g")[1].replace("_dom500","")
    sub="hotspots_dom500_p80_chr1" if "dom500" in reg else "hotspots_p80_chr1"
    return f"{CTRL}/sims/cov10_n{n}_g{g}_s42_{sub}/recomb_truth.tsv.gz"
def est_path(reg,w):
    od="cactus_em_global_filt2_mb" if w=="inv_mb" else "cactus_em_global_filt2_uniform"
    wt="filt2mb" if w=="inv_mb" else "filt2u"
    return f"{CTRL}/results/{od}/{reg}/p80_{wt}_{reg}_cov10_s42.tsv"
def met(est,tru):
    d=est-tru; ss=np.sum(d**2); st=np.sum((tru-tru.mean())**2)
    return dict(n=len(d),MAE=np.mean(np.abs(d)),RMSE=math.sqrt(np.mean(d**2)),
                R2=(1-ss/st) if st>0 else float("nan"),bias=np.mean(d),
                outlier=np.mean(np.abs(d)>0.10))
rows=[]
for reg in REGIMES:
    tp=truth_path(reg)
    if not os.path.exists(tp): print(f"[skip] {reg}"); continue
    tr=pd.read_csv(tp,sep="\t").dropna(subset=["truth_af"]).drop_duplicates(KEYS)
    tr["vclass"]=vclass(tr.ref_len.values,tr.alt_len.values)
    for w in ["uniform","inv_mb"]:
        ep=est_path(reg,w)
        if not os.path.exists(ep): print(f"[miss est] {reg} {w}"); continue
        es=pd.read_csv(ep,sep="\t").drop_duplicates(KEYS)
        es["miss_frac"]=(F-es["n_called"])/F
        m=tr.merge(es[KEYS+["alt_freq","miss_frac"]],on=KEYS,how="inner")
        for cond,thr in CONDS:
            mm=m[m.miss_frac<=thr]
            for cls in ["ALL","SNP","indel","SV"]:
                sub=mm if cls=="ALL" else mm[mm.vclass==cls]
                if len(sub)==0: continue
                mt=met(sub.alt_freq.values.astype(float),sub.truth_af.values.astype(float))
                rows.append(dict(regime=reg,weight=w,miss_cond=cond,cls=cls,**mt))
df=pd.DataFrame(rows)
out=f"{CTRL}/results/filt2_mb_missingness_summary.tsv"
df.to_csv(out,sep="\t",index=False); print("wrote",out,df.shape)
piv=df[(df.cls=='ALL')].pivot_table(index="regime",columns=["miss_cond","weight"],values="MAE")
pd.set_option("display.width",240,"display.max_columns",40)
print("\n=== ALL-class MAE: regime x (miss_cond, weight) ===");
print(piv.reindex(columns=["full","miss_le_90","miss_le_50"],level=0).round(5))
