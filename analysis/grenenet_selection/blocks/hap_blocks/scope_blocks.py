"""Panel-only scoping: block-size tradeoff for a haplotype-level kMate.
Per block, on the 231-founder Chr1 filt2inv panel, measure
  (a) # DISTINCT k-mer haplotypes among the 231 founders (what kMate can resolve;
      the 231 -> N identifiability collapse), and
  (b) # panel k-mers in the block (the per-block identifiability budget).
Compare candidate block schemes (fixed sizes + the hapFIRE/grenenet LD partitions)
to find the compromise: distinct-haplotypes meaningfully < 231 while k-mers/block
stay adequate. NO cohort run, NO kMate — just the panel matrix.
"""
import numpy as np, scipy.sparse as sp, os, json, time

ROOT="/global/scratch/users/tbellg/kmate"
PAN=f"{ROOT}/data/kmer_pa_231_arch3_filt2inv/kmer_pa_Chr1"
OUT=f"{ROOT}/analysis/grenenet_selection/blocks/hap_blocks"; os.makedirs(OUT,exist_ok=True)

t0=time.time()
meta=np.load(f"{PAN}.meta.npz",allow_pickle=True)
fo=np.asarray(meta["founders"]).astype(str); F=len(fo)
bs=np.asarray(meta["bubble_start"]).astype(np.int64); be=np.asarray(meta["bubble_end"]).astype(np.int64)
pos=((bs+be)//2)                                  # per-k-mer genomic position
K=sp.load_npz(f"{PAN}.kmer_pa.npz").astype(bool).tocsc()
print(f"loaded Chr1 panel: F={F} K={K.shape[1]:,} pos {pos.min()}-{pos.max()} [{time.time()-t0:.0f}s]",flush=True)

def distinct_haps(cols):
    """# distinct founder k-mer-presence patterns over `cols` (block k-mers)."""
    if len(cols)==0: return 0,0
    sub=K[:,cols].tocsr()
    sigs=set()
    ip=sub.indptr; idx=sub.indices
    for i in range(F):
        sigs.add(hash(idx[ip[i]:ip[i+1]].tobytes()))
    return len(sigs), len(cols)

def blocks_fixed(win):
    lo,hi=pos.min(),pos.max()
    edges=np.arange(lo, hi+win, win)
    bidx=np.searchsorted(edges, pos, side="right")-1
    return [np.flatnonzero(bidx==b) for b in range(bidx.max()+1)]

def blocks_from_partition(path, chrom="1"):
    out=[]
    for ln in open(path):
        p=ln.split()
        if p[0] not in (chrom,f"Chr{chrom}"): continue
        a,b=int(p[1]),int(p[2])
        out.append(np.flatnonzero((pos>=a)&(pos<=b)))
    return out

def summarize(name, blocks):
    nd=[]; nk=[]
    for cols in blocks:
        if len(cols)==0: continue
        d,k=distinct_haps(cols); nd.append(d); nk.append(k)
    nd=np.array(nd); nk=np.array(nk)
    r=dict(scheme=name, n_blocks=int(len(nd)),
           dhap_median=float(np.median(nd)), dhap_mean=float(nd.mean()),
           dhap_p90=float(np.percentile(nd,90)), dhap_max=int(nd.max()),
           frac_blocks_lt200=float((nd<200).mean()), frac_lt150=float((nd<150).mean()),
           kmers_median=float(np.median(nk)), kmers_p10=float(np.percentile(nk,10)))
    print(f"[{name:16s}] blocks={r['n_blocks']:5d}  distinct-hap/block med={r['dhap_median']:.0f} "
          f"mean={r['dhap_mean']:.0f} p90={r['dhap_p90']:.0f} max={r['dhap_max']}  "
          f"(<200: {100*r['frac_blocks_lt200']:.0f}%, <150: {100*r['frac_lt150']:.0f}%)  "
          f"kmers/block med={r['kmers_median']:.0f} p10={r['kmers_p10']:.0f}  [{time.time()-t0:.0f}s]",flush=True)
    return r

res=[]
for win,name in [(50_000,"fixed_50kb"),(100_000,"fixed_100kb"),(250_000,"fixed_250kb"),
                 (500_000,"fixed_500kb"),(1_000_000,"fixed_1Mb"),(2_000_000,"fixed_2Mb")]:
    res.append(summarize(name, blocks_fixed(win)))
GRENE="/global/scratch/projects/fc_moilab/projects/grenenet-phase1/drive_zenodo/data-intermediate/1001g_grenet_genomewide_partition.txt"
if os.path.exists(GRENE):
    res.append(summarize("grenenet_bigLD", blocks_from_partition(GRENE)))
json.dump(res, open(f"{OUT}/scope_blocks_chr1.json","w"), indent=1)
print(f"\nwrote {OUT}/scope_blocks_chr1.json  ({time.time()-t0:.0f}s)")
