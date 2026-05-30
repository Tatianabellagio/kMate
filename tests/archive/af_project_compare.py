#!/usr/bin/env python3
"""Project h (estimated + true) through the latest arch3 var_pa -> per-record AF,
compare est vs truth. For the 2 finalist methods on the g0 replicate sims.
AF_r = (h . var_pa)_r / (h . var_called)_r   (MAR projection)."""
import numpy as np, glob, os, math
from pathlib import Path
from scipy.sparse import load_npz
ROOT=str(Path(__file__).resolve().parents[2])
CV=f"{ROOT}/panel/arch3/chr1/var_pa_231_arch3_chr1"
print("loading var_pa ...",flush=True)
var_pa=load_npz(f"{CV}.var_pa.npz").tocsc()
cn_called=load_npz(f"{CV}.var_called.npz").tocsc()
meta=np.load(f"{CV}.meta.npz",allow_pickle=True)
fo=np.asarray(meta["founders"]).astype(str)
rl=np.asarray(meta["ref_len"]); al=np.asarray(meta["alt_len"])
is_snp=(rl==1)&(al==1); is_sv=(rl>=50)|(al>=50)
fidx={f:i for i,f in enumerate(fo)}
R=var_pa.shape[1]; print(f"  founders={len(fo)} records={R:,} SNP={is_snp.sum():,} SV={is_sv.sum():,}",flush=True)

def proj(h):
    num=var_pa.T.dot(h); den=cn_called.T.dot(h)
    return num,den

def hvec_from_truth(sim):
    h=np.zeros(len(fo))
    for ln in open(f"{ROOT}/sims/visor_freqk/g0_reps/{sim}/h_truth.tsv"):
        p=ln.split()
        if p and p[0]!="founder" and p[0] in fidx: h[fidx[p[0]]]=float(p[1])
    return h

METHODS={"fact w=1/m_b global":"fact_w=1m_b_global","subsamp+1/m_b global":"subsamp1m_b_global"}
def regime(s): return "n231" if "n231" in s else ("n50_cact" if "cact" in s else "n50_pg")
sims=sorted(set(os.path.basename(p).split("__")[1][:-4]
            for p in glob.glob(f"{ROOT}/scratch/h_fixes/replicates/h/filt2__*.npz")))

out=f"{ROOT}/scratch/h_fixes/replicates/af_results.tsv"
with open(out,"w") as fh:
    fh.write("method\tsim\tregime\tAF_MAE\tAF_RMSE\tAF_MAE_SNP\tAF_MAE_SV\tn_valid\n")
    for sim in sims:
        ht=hvec_from_truth(sim); numt,dent=proj(ht)
        vt=dent>1e-9; truth=np.where(vt,numt/np.maximum(dent,1e-12),np.nan)
        for mlab,tag in METHODS.items():
            hp=f"{ROOT}/scratch/h_fixes/replicates/h/{tag}__{sim}.npz"
            if not os.path.exists(hp): continue
            he=np.asarray(np.load(hp,allow_pickle=True)["h"]).astype(float)
            nume,dene=proj(he); ve=dene>1e-9; est=np.where(ve,nume/np.maximum(dene,1e-12),np.nan)
            valid=vt&ve; d=np.abs(est-truth)
            mae=np.nanmean(d[valid]); rmse=math.sqrt(np.nanmean((d[valid])**2))
            msnp=np.nanmean(d[valid&is_snp]); msv=np.nanmean(d[valid&is_sv])
            fh.write(f"{mlab}\t{sim}\t{regime(sim)}\t{mae:.6f}\t{rmse:.6f}\t{msnp:.6f}\t{msv:.6f}\t{valid.sum()}\n")
            print(f"  {sim:<22}{mlab:<22} AF_MAE={mae:.5f} RMSE={rmse:.5f} SNP={msnp:.5f} SV={msv:.5f}",flush=True)
print(f"wrote {out}",flush=True)
