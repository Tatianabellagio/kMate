#!/usr/bin/env python3
"""Run the 4 focus methods on replicate g0 sims, GLOBAL normalizer.
For one cn base (filt2 or subsamp): load cn once, then for each sim count reads
against the base kmer_index and run EM with w in {1, 1/m_b}. Append scored rows.
Methods: filt2 / fact-w=1/m_b-global (base=filt2); subsampMedian / subsamp+1/m_b-global (base=subsamp)."""
import sys, os, glob, csv, json, math, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))+"/../src")
import numpy as np
from scipy.sparse import load_npz
from kmer_count import count_kmers_in_fasta

ROOT="/carnegie/nobackup/scratch/tbellagio/hapfire_sv"
base=sys.argv[1]                       # filt2 | subsamp
sims=sorted(glob.glob(sys.argv[2]))    # sim dirs
out_tsv=sys.argv[3]
CAC=set(map(str,json.load(open(f"{ROOT}/data/founder_split_cactus_pg.json"))["cactus"]))
CNDIR={"filt2":"cn_full_231_v3qc_v3_filt2","subsamp":"cn_full_231_v3qc_v3_subsampMedian_refilt2"}[base]
NAME={"filt2":{"u":"filt2","mb":"fact w=1/m_b global"},
      "subsamp":{"u":"subsampMedian","mb":"subsamp+1/m_b global"}}[base]

print(f"[{base}] loading cn ...",flush=True)
cn=load_npz(f"{ROOT}/poolfreq/data/{CNDIR}/cn_Chr1.cn.npz").tocsr()
meta=np.load(f"{ROOT}/poolfreq/data/{CNDIR}/cn_Chr1.meta.npz",allow_pickle=True)
founders=np.asarray(meta["founders"]).astype(str)
kmer_index=list(np.asarray(meta["kmer_index"]).astype(str))
bid=np.asarray(meta["bubble_id"]).astype(np.int64); m_b=np.bincount(bid)[bid].astype(np.float64)
F,K=cn.shape; print(f"  F={F} K={K:,}",flush=True)

def em_global(cn_nz,counts_nz,omega,max_iter=300,tol=1e-7):
    h=np.full(cn_nz.shape[0],1.0/cn_nz.shape[0],np.float32)
    wc=(omega*counts_nz).astype(np.float32); tot=float(wc.sum())
    for _ in range(max_iter):
        mu=np.maximum(h@cn_nz,1e-7); num=h*(cn_nz@(wc/mu)); hn=num/tot; s=hn.sum()
        if s<=0 or not np.isfinite(s): break
        hn/=s
        if np.linalg.norm(hn-h)<tol: h=hn; break
        h=hn
    return h.astype(np.float64)

def score(h,truthd):
    t=np.array([truthd.get(x,np.nan) for x in founders]); ok=~np.isnan(t)
    h2,t2,fo=h[ok],t[ok],founders[ok]; cac=np.array([x in CAC for x in fo])
    return float(h2[cac].sum()), math.sqrt(np.mean((h2-t2)**2)), float(h2[t2==0].sum())

HDIR=f"{ROOT}/scratch/h_fixes/replicates/h"; os.makedirs(HDIR,exist_ok=True)
CDIR=f"{ROOT}/scratch/h_fixes/replicates/counts"; os.makedirs(CDIR,exist_ok=True)
rows=[]
for sd in sims:
    name=os.path.basename(sd); r1,r2=f"{sd}/reads/r1.fq",f"{sd}/reads/r2.fq"
    if not os.path.exists(r1): print(f"  skip {name} (no reads)"); continue
    th={}
    for ln in open(f"{sd}/h_truth.tsv"):
        p=ln.split()
        if p and p[0]!="founder": th[p[0]]=float(p[1])
    cpath=f"{CDIR}/{base}_{name}.counts.npy"
    if os.path.exists(cpath):
        counts=np.load(cpath)
    else:
        t=time.time(); cd=count_kmers_in_fasta([r1,r2],kmer_index,k=31,threads=8,hash_size="3G")
        counts=np.array([cd[km] for km in kmer_index],dtype=np.int64); np.save(cpath,counts)
    nz=counts>0; cn_nz=np.asarray(cn[:,nz].todense(),np.float32); cnz=counts[nz].astype(np.float32); mbz=m_b[nz]
    print(f"  {name}: counted nz={nz.sum():,} [{time.time()-t:.0f}s]",flush=True)
    for wkey,om in [("u",np.ones(nz.sum(),np.float32)),("mb",(1.0/mbz).astype(np.float32))]:
        h=em_global(cn_nz,cnz,om); cm,rmse,spur=score(h,th)
        tag=NAME[wkey].replace(" ","_").replace("/","").replace("+","")
        np.savez(f"{HDIR}/{tag}__{name}.npz",founders=founders,h=h)
        rows.append((NAME[wkey],name,round(cm,3),round(rmse,5),round(spur,3)))
        print(f"    {NAME[wkey]:<24} cactus={cm:.3f} RMSE={rmse:.5f} spur={spur:.3f}",flush=True)
    del cn_nz
hdr=not os.path.exists(out_tsv)
with open(out_tsv,"a") as f:
    if hdr: f.write("method\tsim\tcactus\tRMSE\tspurious\n")
    for r in rows: f.write("\t".join(map(str,r))+"\n")
print(f"[{base}] wrote {len(rows)} rows -> {out_tsv}",flush=True)
