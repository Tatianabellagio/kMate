import numpy as np, scipy.sparse as sp, sys
sys.path.insert(0,"src")
import pandas as pd
from kmate.em_solver import solve_em
from kmate.block_em import BlockSpec, assign_kmers_to_blocks, solve_em_per_block
from kmate.kmer_count import count_kmers_in_fasta

ROOT="/global/scratch/users/tbellg/kmate"
pfx=f"{ROOT}/data/kmer_pa_231_arch3_filt2inv/kmer_pa_Chr1"
Ksp=sp.load_npz(pfx+".kmer_pa.npz").tocsc(); F=Ksp.shape[0]
m=np.load(pfx+".meta.npz",allow_pickle=True); ki=m["kmer_index"]
bid=np.asarray(m["bubble_id"]).astype(np.int64); bchrom=np.asarray(m["bubble_chrom"]).astype(str)
bstart=np.asarray(m["bubble_start"]).astype(np.int64); bend=np.asarray(m["bubble_end"]).astype(np.int64)
sub="cov10_n231_g0_s42_hotspots_p231_chr1"
cd=count_kmers_in_fasta([f"{ROOT}/benchmarks/p231/sims/{sub}/reads/r1.fq",f"{ROOT}/benchmarks/p231/sims/{sub}/reads/r2.fq"],list(ki),k=31,threads=8,hash_size="3G")
counts=np.array([cd[k] for k in ki],dtype=np.float64)
bt=pd.read_csv(f"{ROOT}/benchmarks/p231/data/ld_blocks_r2_0.10_Chr1.tsv",sep="\t")
blocks=[BlockSpec("Chr1",int(s),int(e)) for s,e in zip(bt.start_pos,bt.end_pos)]
kmer_block=assign_kmers_to_blocks(bid,bchrom,bstart,bend,blocks)
kfw_fullchrom=np.asarray(Ksp.sum(1)).ravel().astype(np.float32)
effn=lambda h:1.0/np.sum(h**2)

for b in (7,0):
    idxs=np.flatnonzero(kmer_block==b); msk=counts[idxs]>0; nz=idxs[msk]
    Kd=np.asarray(Ksp[:,idxs].todense(),dtype=np.float32); Knz=np.ascontiguousarray(Kd[:,msk]); c=counts[nz].astype(np.float32)
    kfw_block=(Kd@np.ones(len(idxs),np.float32)).astype(np.float32)
    print(f"\n===== blk{b}  observed_nz={len(nz):,} =====",flush=True)
    for lab,kw in [("per_founder + BLOCK kfw (prod)",dict(normalize="per_founder",kfw=kfw_block)),
                   ("per_founder + FULL-CHROM kfw",   dict(normalize="per_founder",kfw=kfw_fullchrom)),
                   ("multinomial (OLD)",              dict(normalize="global")),
                   ("per_founder + Dirichlet a=1.05", dict(normalize="per_founder",kfw=kfw_block,dirichlet_alpha=1.05)),
                   ("per_founder + Dirichlet a=1.5",  dict(normalize="per_founder",kfw=kfw_block,dirichlet_alpha=1.5)),
                   ("per_founder + anchor->uniform w=0.3", dict(normalize="per_founder",kfw=kfw_block,prior_h=np.full(F,1/F,np.float32),prior_weight=0.3))]:
        h,info=solve_em(c,Knz,1.0,max_iter=600,tol=1e-8,omega=None,**kw)
        print(f"  {lab:40s} eff_n={effn(h):6.1f} iters={info['iterations']} conv={info['converged']}",flush=True)

# equivalence: solve_em_per_block (1 block = blk7) vs standalone per_founder+block-kfw
print("\n===== equivalence: block-path vs standalone (blk7) =====",flush=True)
Kdense=np.asarray(Ksp.todense(),dtype=np.float32)
kb=np.full(Ksp.shape[1],-1,np.int32); idx7=np.flatnonzero(kmer_block==7); kb[idx7]=0
hbk,st,gh=solve_em_per_block(counts.astype(np.float32),Kdense,kb,1,1.0,em_max_iter=600,tol=1e-8,
                             min_kmers_per_block=1,local_only=True,normalize="per_founder")
msk=counts[idx7]>0
Kd7=np.ascontiguousarray(Kdense[:,idx7][:,msk]); kfw7=(Kdense[:,idx7]@np.ones(len(idx7),np.float32)).astype(np.float32)
hstd,_=solve_em(counts[idx7[msk]].astype(np.float32),Kd7,1.0,max_iter=600,tol=1e-8,normalize="per_founder",kfw=kfw7)
print(f"  block-path eff_n={effn(hbk[0]):.1f}  standalone eff_n={effn(hstd):.1f}  max|Δh|={np.abs(hbk[0]-hstd).max():.2e}",flush=True)
print("DONE-BATTERY")
