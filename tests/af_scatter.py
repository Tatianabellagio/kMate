import numpy as np, os, math
from pathlib import Path
from scipy.sparse import load_npz
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
ROOT=str(Path(__file__).resolve().parents[2]); CV=f"{ROOT}/arch3/chr1/cn_var_231_arch3_chr1"
print("loading cn_var...",flush=True)
cn_var=load_npz(f"{CV}.cn_var.npz").tocsc(); cn_called=load_npz(f"{CV}.cn_var_called.npz").tocsc()
fo=np.asarray(np.load(f"{CV}.meta.npz",allow_pickle=True)["founders"]).astype(str)
fidx={f:i for i,f in enumerate(fo)}
def proj(h): return cn_var.T.dot(h), cn_called.T.dot(h)
def htrue(sim):
    h=np.zeros(len(fo))
    for ln in open(f"{ROOT}/sims/visor_freqk/g0_reps/{sim}/h_truth.tsv"):
        p=ln.split()
        if p and p[0]!="founder" and p[0] in fidx: h[fidx[p[0]]]=float(p[1])
    return h
METH=[("fact w=1/m_b global","fact_w=1m_b_global"),("subsamp+1/m_b global","subsamp1m_b_global")]
COLS=[("n=231","g0_n231_rep101_rand"),("n=50 cact","g0_n50_rep101_cact"),("n=50 pg","g0_n50_rep101_pg")]
fig,ax=plt.subplots(len(METH),len(COLS),figsize=(5*len(COLS),4.5*len(METH)))
for ci,(rlab,sim) in enumerate(COLS):
    nt,dt=proj(htrue(sim)); vt=dt>1e-9; truth=np.where(vt,nt/np.maximum(dt,1e-12),np.nan)
    for ri,(mlab,tag) in enumerate(METH):
        a=ax[ri,ci]; hp=f"{ROOT}/scratch/h_fixes/replicates/h/{tag}__{sim}.npz"
        he=np.asarray(np.load(hp,allow_pickle=True)["h"]).astype(float)
        ne,de=proj(he); ve=de>1e-9; est=np.where(ve,ne/np.maximum(de,1e-12),np.nan)
        m=vt&ve&~np.isnan(truth)&~np.isnan(est); x=truth[m]; y=est[m]
        mae=np.abs(x-y).mean(); r2=np.corrcoef(x,y)[0,1]**2
        a.hexbin(x,y,gridsize=80,bins="log",cmap="viridis",mincnt=1)
        a.plot([0,1],[0,1],"r--",lw=1)
        a.set_xlim(0,1); a.set_ylim(0,1); a.set_xlabel("truth AF"); a.set_ylabel("estimated AF")
        a.set_title(f"{rlab} — {mlab}\nAF MAE={mae:.5f}  R2={r2:.4f}  n={m.sum():,}",fontsize=9)
fig.tight_layout(); out=f"{ROOT}/notebook/plots/af_true_vs_est.png"
fig.savefig(out,dpi=110,bbox_inches="tight"); print("wrote",out)
