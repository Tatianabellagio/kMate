import numpy as np, scipy.sparse as sp, sys
sys.path.insert(0,"src")
import pandas as pd
from kmate.em_solver import solve_em
from kmate.block_em import BlockSpec, assign_kmers_to_blocks
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
effn=lambda h:1.0/np.sum(h**2)

def fit(idxs, mi):
    msk=counts[idxs]>0; nz=idxs[msk]
    Kd=np.ascontiguousarray(np.asarray(Ksp[:,nz].todense(),dtype=np.float32))
    kfw=np.asarray(Ksp[:,idxs].sum(1)).ravel().astype(np.float32)
    h,info=solve_em(counts[nz].astype(np.float32),Kd,1.0,max_iter=mi,tol=1e-7,normalize="per_founder",kfw=kfw)
    return len(nz),info["iterations"],info["converged"],effn(h)

print("PRODUCTION settings (tol=1e-7, max_iter=200)")
print("unit          n_nz        iters  conv    eff_n")
allidx=np.flatnonzero(kmer_block>=0)
n,it,cv,en=fit(allidx,200); print(f"CHROM(all)   {n:>9,}  {it:>4d}  {str(cv):>5s}  {en:6.1f}")
for b in range(len(blocks)):
    n,it,cv,en=fit(np.flatnonzero(kmer_block==b),200)
    print(f"blk{b:<2d}         {n:>9,}  {it:>4d}  {str(cv):>5s}  {en:6.1f}")
print("\nASYMPTOTE (max_iter=3000): does chrom converge if allowed?")
for name,idx in [("CHROM",allidx),("blk0",np.flatnonzero(kmer_block==0)),("blk7",np.flatnonzero(kmer_block==7))]:
    n,it,cv,en=fit(idx,3000); print(f"  {name:6s} iters={it} conv={cv} eff_n={en:.1f}")
print("DONE-CONVMAP")
