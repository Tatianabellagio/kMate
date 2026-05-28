#!/usr/bin/env python3
"""Score filt2+1/m_b (inv_mb) vs filt2 (uniform) global EM on p80 recomb sims.
JOIN BY RECORD INDEX: est, truth, cn_var are all in identical cn_var record order,
so est[i]<->truth[i]. The (chrom,pos,ref_len,alt_len) 4-key is NOT unique (p80 cn_var
has 4.6% dup keys; atomized 231 has 38%) -> a key-merge mis-joins multiallelic ALTs.
We assert identical (pos,ref_len,alt_len) row order then compare by position.
Class from ref_len/alt_len. Writes filt2_mb_vs_uniform_summary.tsv."""
import os, math
import numpy as np, pandas as pd
CTRL="/carnegie/nobackup/scratch/tbellagio/hapfire_sv/control_p80"
REGIMES=["n50_g0","n80_g0","n50_g1","n80_g1","n50_g3","n50_g3_dom500"]
def subdir(reg):
    if reg=="n50_g3_dom500": return "cov10_n50_g3_s42_hotspots_dom500_p80_chr1"
    n,g=reg.split("_g"); return f"cov10_{n}_g{g}_s42_hotspots_p80_chr1"
def est_path(reg,w):
    od="cactus_em_global_filt2_mb" if w=="inv_mb" else "cactus_em_global_filt2_uniform"
    wt="filt2mb" if w=="inv_mb" else "filt2u"
    return f"{CTRL}/results/{od}/{reg}/p80_{wt}_{reg}_cov10_s42.tsv"
def met(e,t):
    d=e-t; ss=np.sum(d**2); st=np.sum((t-t.mean())**2)
    return dict(n=len(d),MAE=np.mean(np.abs(d)),RMSE=math.sqrt(np.mean(d**2)),
                R2=(1-ss/st) if st>0 else float("nan"),bias=np.mean(d),outlier=np.mean(np.abs(d)>0.10))
rows=[]
for reg in REGIMES:
    tp=f"{CTRL}/sims/{subdir(reg)}/recomb_truth.tsv.gz"
    if not os.path.exists(tp): print(f"[skip] {reg}: no truth"); continue
    tr=pd.read_csv(tp,sep="\t")
    rl=tr["ref_len"].values.astype(int); al=tr["alt_len"].values.astype(int)
    vcls=np.where(np.maximum(rl,al)>=50,"SV",np.where((rl==1)&(al==1),"SNP","indel"))
    for w in ["uniform","inv_mb"]:
        ep=est_path(reg,w)
        if not os.path.exists(ep): print(f"[wait] {reg} {w}: est not ready"); continue
        es=pd.read_csv(ep,sep="\t")
        if len(es)!=len(tr): print(f"[ERR] {reg} {w}: len mismatch -- SKIP"); continue
        if not (np.array_equal(es["pos"].values,tr["pos"].values)
                and np.array_equal(es["ref_len"].values.astype(int),rl)
                and np.array_equal(es["alt_len"].values.astype(int),al)):
            print(f"[ERR] {reg} {w}: row order mismatch -- SKIP"); continue
        e=es["alt_freq"].values; t=tr["truth_af"].values; valid=np.isfinite(e)&np.isfinite(t)
        for cls in ["ALL","SNP","indel","SV"]:
            msk=valid if cls=="ALL" else (valid&(vcls==cls))
            if msk.sum()==0: continue
            rows.append(dict(regime=reg,weight=w,cls=cls,**met(e[msk],t[msk])))
df=pd.DataFrame(rows)
out=f"{CTRL}/results/filt2_mb_vs_uniform_summary.tsv"; df.to_csv(out,sep="\t",index=False)
print("wrote",out,df.shape,"(record-index join)")
piv=df[df.cls=="ALL"].pivot(index="regime",columns="weight",values=["MAE","RMSE","outlier"])
pd.set_option("display.width",200,"display.max_columns",30)
print("\n=== ALL-class, uniform vs inv_mb (record-index join) ==="); print(piv.round(5))
