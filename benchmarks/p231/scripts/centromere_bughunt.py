import numpy as np, scipy.sparse as sp, sys, time
sys.path.insert(0, "src")
import pandas as pd
from kmate.em_solver import solve_em, haploblock_collapse_indices
from kmate.block_em import BlockSpec, assign_kmers_to_blocks
from kmate.kmer_count import count_kmers_in_fasta

ROOT="/global/scratch/users/tbellg/kmate"
pfx=f"{ROOT}/data/kmer_pa_231_arch3_filt2inv/kmer_pa_Chr1"
K=sp.load_npz(pfx+".kmer_pa.npz").tocsc()
m=np.load(pfx+".meta.npz",allow_pickle=True)
ki=m["kmer_index"]
bid=np.asarray(m["bubble_id"]).astype(np.int64); bchrom=np.asarray(m["bubble_chrom"]).astype(str)
bstart=np.asarray(m["bubble_start"]).astype(np.int64); bend=np.asarray(m["bubble_end"]).astype(np.int64)
F=K.shape[0]

# count reads
sub="cov10_n231_g0_s42_hotspots_p231_chr1"
r1=f"{ROOT}/benchmarks/p231/sims/{sub}/reads/r1.fq"; r2=f"{ROOT}/benchmarks/p231/sims/{sub}/reads/r2.fq"
t=time.time(); cd=count_kmers_in_fasta([r1,r2],list(ki),k=31,threads=8,hash_size="3G")
counts=np.array([cd[k] for k in ki],dtype=np.float64); print(f"counted in {time.time()-t:.0f}s; total nz {int((counts>0).sum()):,}",flush=True)

bt=pd.read_csv(f"{ROOT}/benchmarks/p231/data/ld_blocks_r2_0.10_Chr1.tsv",sep="\t")
blocks=[BlockSpec("Chr1",int(s),int(e)) for s,e in zip(bt.start_pos,bt.end_pos)]
kmer_block=assign_kmers_to_blocks(bid,bchrom,bstart,bend,blocks)

def effn(h): return 1.0/np.sum(h**2)
print("\n=== observed (nonzero) k-mer count per block ===")
for b in range(len(blocks)):
    idxs=np.flatnonzero(kmer_block==b); nz=idxs[counts[idxs]>0]
    print(f"  blk{b:2d} panel={len(idxs):>9,}  observed_nz={len(nz):>9,}  ({100*len(nz)/max(1,len(idxs)):.1f}%)")

print("\n=== CENTROMERE blk7: fit the SAME k-mers several ways ===")
b=7; idxs=np.flatnonzero(kmer_block==b); nz=idxs[counts[idxs]>0]
Kb_dense=np.asarray(K[:,idxs].todense(),dtype=np.float32)
Knz=np.ascontiguousarray(Kb_dense[:,counts[idxs]>0]); c_nz=counts[nz].astype(np.float32)
kfw=(Kb_dense@np.ones(len(idxs),np.float32)).astype(np.float32)
print(f"  blk7 panel={len(idxs):,} observed_nz={len(nz):,}")
for lab,kw in [("NEW per_founder+kfw", dict(normalize="per_founder",kfw=kfw)),
               ("NEW per_founder (no kfw arg)", dict(normalize="per_founder")),
               ("OLD multinomial(global)", dict(normalize="global")),
               ("NEW per_founder 800it", dict(normalize="per_founder",kfw=kfw,max_iter=800))]:
    mi=kw.pop("max_iter",300)
    h,info=solve_em(c_nz,Knz,1.0,max_iter=mi,tol=1e-7,omega=None,**kw)
    print(f"    {lab:32s} eff_n={effn(h):6.1f}  iters={info['iterations']}  conv={info['converged']}")
print("DONE-BUGHUNT")
