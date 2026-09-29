import numpy as np, scipy.sparse as sp, sys, time
sys.path.insert(0,"src")
import pandas as pd
from scipy.optimize import nnls
from kmate.em_solver import solve_em
from kmate.block_em import BlockSpec, assign_kmers_to_blocks
from kmate.kmer_count import count_kmers_in_fasta
ROOT="/global/scratch/users/tbellg/kmate"
pfx=f"{ROOT}/data/kmer_pa_231_arch3_filt2inv/kmer_pa_Chr1"
Ksp=sp.load_npz(pfx+".kmer_pa.npz").tocsc(); F=Ksp.shape[0]
m=np.load(pfx+".meta.npz",allow_pickle=True); ki=m["kmer_index"]
bid=np.asarray(m["bubble_id"]).astype(np.int64); bchrom=np.asarray(m["bubble_chrom"]).astype(str)
bstart=np.asarray(m["bubble_start"]).astype(np.int64); bend=np.asarray(m["bubble_end"]).astype(np.int64)
bt=pd.read_csv(f"{ROOT}/benchmarks/p231/data/ld_blocks_r2_0.10_Chr1.tsv",sep="\t")
blocks=[BlockSpec("Chr1",int(s),int(e)) for s,e in zip(bt.start_pos,bt.end_pos)]
kmer_block=assign_kmers_to_blocks(bid,bchrom,bstart,bend,blocks)
# var_pa (raw) for AF projection + truth
vp=sp.load_npz(f"{ROOT}/panel/arch3/chr1/var_pa_231_arch3_chr1.var_pa.npz").astype(np.float32).tocsr()
vc=sp.load_npz(f"{ROOT}/panel/arch3/chr1/var_pa_231_arch3_chr1.var_called.npz").astype(np.float32).tocsr()
effn=lambda h:1.0/np.sum(h**2)
def af_mae(h,truth):
    num=h@vp; den=np.maximum(h@vc,1e-9); af=np.asarray(num/den).ravel()
    ok=np.isfinite(af)&np.isfinite(truth); return float(np.mean(np.abs(af[ok]-truth[ok])))

def fit_all(c_nz,Knz,kfw,gh):
    R={}
    R["EM current"]=solve_em(c_nz,Knz,1.0,max_iter=200,tol=1e-7,normalize="per_founder",kfw=kfw)[0]
    R["anchor->unif w0.3"]=solve_em(c_nz,Knz,1.0,max_iter=400,tol=1e-7,normalize="per_founder",kfw=kfw,prior_h=np.full(F,1/F,np.float32),prior_weight=0.3)[0]
    R["anchor->global w0.3"]=solve_em(c_nz,Knz,1.0,max_iter=400,tol=1e-7,normalize="per_founder",kfw=kfw,prior_h=gh.astype(np.float32),prior_weight=0.3)[0]
    R["Dirichlet 1.5"]=solve_em(c_nz,Knz,1.0,max_iter=400,tol=1e-7,normalize="per_founder",kfw=kfw,dirichlet_alpha=1.5)[0]
    # HARP-style NNLS: min ||K^T h - (c/lambda)||, h>=0, then renorm to simplex
    lam=c_nz.sum()/max(Knz.sum(),1e-9); y=(c_nz/max(lam,1e-9)).astype(np.float64)
    hq,_=nnls(Knz.T.astype(np.float64), y, maxiter=500); hq=hq/max(hq.sum(),1e-12)
    R["HARP NNLS"]=hq.astype(np.float32)
    return R

for pool,sub in [("n231_g0","cov10_n231_g0_s42_hotspots_p231_chr1"),
                 ("n50_g0", "cov10_n50_g0_s42_hotspots_p231_chr1")]:
    cd=count_kmers_in_fasta([f"{ROOT}/benchmarks/p231/sims/{sub}/reads/r1.fq",f"{ROOT}/benchmarks/p231/sims/{sub}/reads/r2.fq"],list(ki),k=31,threads=8,hash_size="3G")
    counts=np.array([cd[k] for k in ki],dtype=np.float64)
    truth=pd.read_csv(f"{ROOT}/benchmarks/p231/sims/{sub}/recomb_truth_raw.tsv.gz",sep="\t",usecols=["truth_af"])["truth_af"].values.astype(float)
    pw=pd.read_csv(f"{ROOT}/benchmarks/p231/sims/{sub}/pool_weights.tsv",sep="\t")
    present=set(pw[pw["weight"]>0]["founder"].astype(str)); fo=np.asarray(m["founders"]).astype(str)
    absent_mask=np.array([f not in present for f in fo])
    # global h (for anchor->global): EM on all k-mers
    allidx=np.flatnonzero(kmer_block>=0); nzall=allidx[counts[allidx]>0]
    Kall=np.ascontiguousarray(np.asarray(Ksp[:,nzall].todense(),np.float32)); kfwall=np.asarray(Ksp[:,allidx].sum(1)).ravel().astype(np.float32)
    gh=solve_em(counts[nzall].astype(np.float32),Kall,1.0,max_iter=200,tol=1e-7,normalize="per_founder",kfw=kfwall)[0]
    print(f"\n########## POOL {pool}  (true present={len(present)}, true eff_n~{1/np.sum(pw.weight.values**2 / pw.weight.values.sum()**2):.0f}) ##########",flush=True)
    for unit,idx in [("CHROM",allidx),("blk7(cent)",np.flatnonzero(kmer_block==7))]:
        nz=idx[counts[idx]>0]; Knz=np.ascontiguousarray(np.asarray(Ksp[:,nz].todense(),np.float32)); c=counts[nz].astype(np.float32)
        kfw=np.asarray(Ksp[:,idx].sum(1)).ravel().astype(np.float32)
        print(f"  --- {unit} (n_nz={len(nz):,}) ---",flush=True)
        for name,h in fit_all(c,Knz,kfw,gh).items():
            absmass=float(h[absent_mask].sum())
            print(f"    {name:22s} eff_n={effn(h):6.1f}  AF_MAE={af_mae(h,truth):.4f}  mass_on_ABSENT={absmass:.3f}",flush=True)
print("DONE-REG")
