"""Haploblock reframe with SIMILARITY-THRESHOLD clustering (k-mer native).
Per block: pairwise Hamming-fraction distance between the 231 founders' k-mer
presence signatures -> hierarchical (complete-linkage) clustering at threshold eps
-> haploblocks. eps=0 == exact identity. Sweeps eps x r2 block-size.
Scores each haploblock's frequency vs equimolar truth (TVD, scale-fair).
Post-hoc valid: haploblock freq = sum of member founders' block-EM h.
"""
import numpy as np, scipy.sparse as sp
from scipy.spatial.distance import pdist
from scipy.cluster.hierarchy import linkage, fcluster
ROOT="/global/scratch/users/tbellg/kmate"
OUT=f"{ROOT}/analysis/grenenet_gea/hap_blocks/bench_g0_231"
KMPRE=f"{ROOT}/benchmarks/p231/data/kmer_pa_p231_filt2inv/kmer_pa_Chr1"
F=231
EPS=[0.0,0.01,0.02,0.05]                 # k-mer disagreement fraction to merge founders

km=np.load(f"{KMPRE}.meta.npz",allow_pickle=True)
bpos=((km["bubble_start"].astype(np.int64)+km["bubble_end"].astype(np.int64))//2)
K=sp.load_npz(f"{KMPRE}.kmer_pa.npz").astype(bool).tocsc()
order=np.argsort(bpos); spos=bpos[order]

def block_matrix(a,b):
    lo=np.searchsorted(spos,a,"left"); hi=np.searchsorted(spos,b,"right")
    cols=order[lo:hi]
    if len(cols)==0: return None
    return np.asarray(K[:,cols].todense(),dtype=np.uint8)     # 231 x nk

def labels_at_eps(M,eps):
    """founder haploblock labels merging founders with < eps hamming-frac distance."""
    if eps<=0:                                                # exact identity (fast path)
        _,lab=np.unique([hash(M[i].tobytes()) for i in range(F)],return_inverse=True)
        return lab
    d=pdist(M,metric="hamming")                               # fraction of k-mers differing
    Z=linkage(d,method="complete")
    return fcluster(Z,t=eps,criterion="distance")-1

print(f"{'r2':5s} {'eps':5s} {'blocks':>6s} {'Kb_med':>6s} {'Kb_p90':>6s} {'TVD_med':>8s} {'TVD_p90':>8s}")
for r2 in ["0.10","0.20","0.30","0.40"]:
    z=np.load(f"{OUT}/alltype_r2_{r2}.h_blocks_per_chrom.npz",allow_pickle=True)
    hb=z["Chr1_h_blocks"]; bs=z["Chr1_block_start"]; be=z["Chr1_block_end"]
    # precompute block matrices + finite-h once
    mats=[]
    for j in range(len(hb)):
        h=hb[j]
        if not np.isfinite(h).all(): mats.append(None); continue
        mats.append((h/h.sum(), block_matrix(int(bs[j]),int(be[j]))))
    for eps in EPS:
        Kbs=[]; tvd=[]
        for mj in mats:
            if mj is None or mj[1] is None: continue
            h,M=mj
            lab=labels_at_eps(M,eps); Kb=lab.max()+1
            est=np.bincount(lab,weights=h,minlength=Kb)
            true=np.bincount(lab,minlength=Kb)/F
            Kbs.append(Kb); tvd.append(0.5*np.sum(np.abs(est-true)))
        Kbs=np.array(Kbs); tvd=np.array(tvd)
        print(f"{r2:5s} {eps:5.2f} {len(Kbs):6d} {np.median(Kbs):6.0f} {np.percentile(Kbs,90):6.0f} "
              f"{np.median(tvd):8.3f} {np.percentile(tvd,90):8.3f}",flush=True)
print()
print("eps=0 is exact identity. Higher eps merges near-identical founders -> fewer")
print("haploblocks (Kb down). Watch TVD: clustering should keep/lower it while cutting Kb.")
