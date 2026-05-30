#!/usr/bin/env python3
"""Score kMate p231 AF estimates vs per-record truth, both var_pa arms.
JOIN BY RECORD INDEX (not the 4-key): est, truth and the var_pa meta are all emitted
in identical var_pa record order, so est[i]<->truth[i]<->var_pa record i. The
(chrom,pos,ref_len,alt_len) 4-key is NOT unique on the multiallelic-heavy atomized
panel (38% dup keys) -> a key-merge mis-joins/collapses multiallelic ALTs. We assert
identical (pos,ref_len,alt_len) row order, then compare by position.
Atomized arm: all SNP-class. Raw arm: SNP/indel/SV. Writes p231_kmate_summary.tsv."""
import os, math
from pathlib import Path
import numpy as np, pandas as pd
CTRL=str(Path(__file__).resolve().parents[1])
REGIMES=["n50_g0","n231_g0","n50_g1","n231_g1","n50_g3","n50_g3_dom500"]
def subdir(reg):
    if reg=="n50_g3_dom500": return "cov10_n50_g3_s42_hotspots_dom500_p231_chr1"
    n,g=reg.split("_g"); return f"cov10_{n}_g{g}_s42_hotspots_p231_chr1"
def met(e,t):
    d=e-t; ss=np.sum(d**2); st=np.sum((t-t.mean())**2)
    return dict(n=len(d),MAE=np.mean(np.abs(d)),RMSE=math.sqrt(np.mean(d**2)),
                R2=(1-ss/st) if st>0 else float("nan"),bias=np.mean(d),outlier=np.mean(np.abs(d)>0.10))
rows=[]
for reg in REGIMES:
    for cnvar in ["atomized","raw"]:
        tp=f"{CTRL}/sims/{subdir(reg)}/recomb_truth_{cnvar}.tsv.gz"
        if not os.path.exists(tp): print(f"[skip] {reg} {cnvar}: no truth"); continue
        tr=pd.read_csv(tp,sep="\t")
        rl=tr["ref_len"].values.astype(int); al=tr["alt_len"].values.astype(int)
        vcls=np.where(np.maximum(rl,al)>=50,"SV",np.where((rl==1)&(al==1),"SNP","indel"))
        for weight,wtag in [("inv_mb","filt2mb"),("uniform","filt2u")]:
            ep=f"{CTRL}/results/kmate_global_{wtag}_{cnvar}/{reg}/p231_{wtag}_{cnvar}_{reg}_cov10_s42.tsv"
            if not os.path.exists(ep): continue
            es=pd.read_csv(ep,sep="\t")
            if len(es)!=len(tr):
                print(f"[ERR] {reg} {cnvar} {weight}: len est={len(es)} != truth={len(tr)} -- SKIP"); continue
            # ORDER GATE: both must be in identical var_pa record order
            if not (np.array_equal(es["pos"].values,tr["pos"].values)
                    and np.array_equal(es["ref_len"].values.astype(int),rl)
                    and np.array_equal(es["alt_len"].values.astype(int),al)):
                print(f"[ERR] {reg} {cnvar} {weight}: row order mismatch -- SKIP"); continue
            e=es["alt_freq"].values; t=tr["truth_af"].values
            valid=np.isfinite(e)&np.isfinite(t)
            classes=["ALL","SNP","indel","SV"] if cnvar=="raw" else ["ALL","SNP"]
            for cls in classes:
                msk=valid if cls=="ALL" else (valid&(vcls==cls))
                if msk.sum()==0: continue
                mt=met(e[msk],t[msk])
                rows.append(dict(regime=reg,cnvar=cnvar,weight=weight,cls=cls,join="index",**mt))
df=pd.DataFrame(rows)
out=f"{CTRL}/results/p231_kmate_summary.tsv"; df.to_csv(out,sep="\t",index=False)
print("wrote",out,df.shape)
if len(df):
    pd.set_option("display.width",240,"display.max_columns",40,"display.max_rows",200)
    print("\n=== ALL-class ==="); print(df[df.cls=="ALL"][["regime","cnvar","weight","n","MAE","RMSE","R2","outlier"]].to_string(index=False))
    print("\n=== raw-arm by class ==="); print(df[(df.cnvar=="raw")&(df.cls!="ALL")][["regime","weight","cls","n","MAE","RMSE","outlier"]].to_string(index=False))
