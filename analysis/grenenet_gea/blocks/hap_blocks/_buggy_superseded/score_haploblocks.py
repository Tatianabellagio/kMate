"""CORRECT reframe scoring: per block, collapse 231 founders -> distinct
HAPLOBLOCKS (founders with identical k-mer presence over the block's k-mers),
then score the HAPLOTYPE-CLASS frequencies -- NOT 231 founders.

Per block b:
  - K_b = # distinct founder k-mer signatures over block b's k-mers (the actual
    number of identifiable things -- 'not 231').
  - true class freq  = (#founders in class)/231   (equimolar g0 truth)
  - est  class freq  = sum of the block's 231-founder h over the class members
    (the EM class-sum is the identifiable haplotype frequency; the within-class
     split is arbitrary, the sum is not).
Reports, per block size: distribution of K_b (resolvable haplotypes) and how well
the haploblock frequencies are recovered (class-freq MAE), vs the equimolar truth.
"""
import numpy as np, scipy.sparse as sp, json
ROOT="/global/scratch/users/tbellg/kmate"
OUT=f"{ROOT}/analysis/grenenet_gea/blocks/hap_blocks/bench_g0_231"
KMPRE=f"{ROOT}/benchmarks/p231/data/kmer_pa_p231_filt2inv/kmer_pa_Chr1"
F=231; TRUE_PER=1.0/F

km=np.load(f"{KMPRE}.meta.npz",allow_pickle=True)
bpos=((km["bubble_start"].astype(np.int64)+km["bubble_end"].astype(np.int64))//2)
K=sp.load_npz(f"{KMPRE}.kmer_pa.npz").astype(bool).tocsc()
order=np.argsort(bpos); spos=bpos[order]

def block_classes(a,b):
    """founder class labels + sizes over k-mers in [a,b]."""
    lo=np.searchsorted(spos,a,"left"); hi=np.searchsorted(spos,b,"right")
    cols=order[lo:hi]
    if len(cols)==0: return None,0
    sub=K[:,cols].tocsr(); ip=sub.indptr; idx=sub.indices
    sig=[hash(idx[ip[i]:ip[i+1]].tobytes()) for i in range(F)]
    _,lab=np.unique(sig,return_inverse=True)
    return lab,len(cols)

print(f"{'config':8s} {'blocks':>6s} {'Kb_med':>6s} {'Kb_p90':>6s} {'Kb_max':>6s} "
      f"{'TVD_med':>8s} {'TVD_p90':>8s} {'worstcls':>8s} {'relerr':>7s} {'nkmer':>6s}")
for r2 in ["0.10","0.20","0.30","0.40"]:
    z=np.load(f"{OUT}/alltype_r2_{r2}.h_blocks_per_chrom.npz",allow_pickle=True)
    hb=z["Chr1_h_blocks"]; bs=z["Chr1_block_start"]; be=z["Chr1_block_end"]
    Kbs=[]; tvd=[]; worst=[]; rel=[]; nk=[]
    for j in range(len(hb)):
        h=hb[j]
        if not np.isfinite(h).all(): continue
        h=h/h.sum()
        lab,ncol=block_classes(int(bs[j]),int(be[j]))
        if lab is None: continue
        Kb=lab.max()+1
        est=np.bincount(lab,weights=h,minlength=Kb)          # est haplotype freq
        true=np.bincount(lab,minlength=Kb)/F                 # true (equimolar)
        Kbs.append(Kb); nk.append(ncol)
        tvd.append(0.5*np.sum(np.abs(est-true)))             # total variation distance (scale-fair)
        worst.append(np.max(np.abs(est-true)))               # worst single haplotype (absolute)
        rel.append(np.mean(np.abs(est-true)/true))           # mean per-class RELATIVE error
    Kbs=np.array(Kbs); tvd=np.array(tvd); worst=np.array(worst); rel=np.array(rel); nk=np.array(nk)
    print(f"r2={r2:4s} {len(Kbs):6d} {np.median(Kbs):6.0f} {np.percentile(Kbs,90):6.0f} "
          f"{Kbs.max():6.0f} {np.median(tvd):8.3f} {np.percentile(tvd,90):8.3f} "
          f"{np.mean(worst):8.2e} {np.median(rel):7.2f} {np.median(nk):6.0f}")
print()
print("Kb  = # distinct haplotypes resolvable per block ('not 231' -- collinearity collapsed).")
print("TVD = total-variation distance between est & true HAPLOBLOCK freq distributions")
print("      (0.5*sum|est-true|; scale-fair -- both sum to 1). Lower = better recovery.")
print("relerr = median per-class |est-true|/true (fractional error of a haploblock's freq).")
print("Sweet spot: Kb << 231 (collinearity won) AND low TVD (haploblock freqs recovered).")
