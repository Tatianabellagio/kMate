"""Measure hapFIRE's OWN blocks on the kMate k-mer panel (chr1).
hapFIRE defines fine bigLD blocks on SNPs; here we ask: over those exact block
boundaries, how many distinct *k-mer* haplotypes does the 231-founder filt2inv
panel resolve? (i.e. what would kMate get if handed hapFIRE's blocks.)
Compares hapFIRE's SNP-based distinct-hap count (block_n_uniq) against the
k-mer-based distinct-hap count on the same intervals.
"""
import numpy as np, scipy.sparse as sp, time
ROOT="/global/scratch/users/tbellg/kmate"
PAN=f"{ROOT}/data/kmer_pa_231_arch3_filt2inv/kmer_pa_Chr1"
t0=time.time()
meta=np.load(f"{PAN}.meta.npz",allow_pickle=True)
F=len(np.asarray(meta["founders"]))
bs=np.asarray(meta["bubble_start"]).astype(np.int64); be=np.asarray(meta["bubble_end"]).astype(np.int64)
pos=((bs+be)//2)
K=sp.load_npz(f"{PAN}.kmer_pa.npz").astype(bool).tocsc()

d=np.load(f"{ROOT}/data/hapfire_block_index.npz", allow_pickle=True)
ch=np.asarray(d["block_chrom"]).astype(str)
m=(ch=="1")|(ch=="Chr1")
a=d["block_pos_start"][m].astype(np.int64); b=d["block_pos_end"][m].astype(np.int64)
nuniq_snp=d["block_n_uniq"][m].astype(int)
print(f"hapFIRE chr1 blocks: {m.sum()}  panel k-mers: {K.shape[1]:,}  [{time.time()-t0:.0f}s]",flush=True)

def distinct_haps(cols):
    if len(cols)==0: return 0
    sub=K[:,cols].tocsr(); ip=sub.indptr; idx=sub.indices
    return len({hash(idx[ip[i]:ip[i+1]].tobytes()) for i in range(F)})

nk=np.empty(m.sum(),int); nhk=np.empty(m.sum(),int)
order=np.argsort(pos)
spos=pos[order]
for j in range(m.sum()):
    lo=np.searchsorted(spos,a[j],"left"); hi=np.searchsorted(spos,b[j],"right")
    cols=order[lo:hi]
    nk[j]=len(cols); nhk[j]=distinct_haps(cols)
def q(x,p): return np.percentile(x,p)
print(f"\nOver hapFIRE's OWN {m.sum()} chr1 blocks (median width ~963bp):")
print(f"  k-mers/block:            med={np.median(nk):.0f} p90={q(nk,90):.0f}  (blocks with 0 k-mers: {(nk==0).mean()*100:.0f}%)")
print(f"  distinct K-MER haps:     med={np.median(nhk):.0f} p90={q(nhk,90):.0f} max={nhk.max()}")
print(f"  distinct SNP haps (hapFIRE block_n_uniq): med={np.median(nuniq_snp):.0f} p90={q(nuniq_snp,90):.0f} max={nuniq_snp.max()}")
print(f"  => k-mer panel resolves {np.median(nhk):.0f} vs hapFIRE's {np.median(nuniq_snp):.0f} distinct haps at the SAME (~1kb) blocks")
print(f"[{time.time()-t0:.0f}s]")
