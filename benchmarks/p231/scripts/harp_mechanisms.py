import numpy as np, scipy.sparse as sp, sys
sys.path.insert(0,"src")
import pandas as pd
from kmate.em_solver import solve_em
from kmate.block_em import BlockSpec, assign_kmers_to_blocks
from kmate.kmer_count import count_kmers_in_fasta
ROOT="/global/scratch/users/tbellg/kmate"
pfx=f"{ROOT}/data/kmer_pa_231_arch3_filt2inv/kmer_pa_Chr1"
Ksp=sp.load_npz(pfx+".kmer_pa.npz").tocsc(); F=Ksp.shape[0]
m=np.load(pfx+".meta.npz",allow_pickle=True); ki=m["kmer_index"]; fo=np.asarray(m["founders"]).astype(str)
bid=np.asarray(m["bubble_id"]).astype(np.int64); bchrom=np.asarray(m["bubble_chrom"]).astype(str)
bstart=np.asarray(m["bubble_start"]).astype(np.int64); bend=np.asarray(m["bubble_end"]).astype(np.int64)
bt=pd.read_csv(f"{ROOT}/benchmarks/p231/data/ld_blocks_r2_0.10_Chr1.tsv",sep="\t")
blocks=[BlockSpec("Chr1",int(s),int(e)) for s,e in zip(bt.start_pos,bt.end_pos)]
kmer_block=assign_kmers_to_blocks(bid,bchrom,bstart,bend,blocks)
vp=sp.load_npz(f"{ROOT}/panel/arch3/chr1/var_pa_231_arch3_chr1.var_pa.npz").astype(np.float32).tocsr()
vc=sp.load_npz(f"{ROOT}/panel/arch3/chr1/var_pa_231_arch3_chr1.var_called.npz").astype(np.float32).tocsr()
effn=lambda h:1.0/np.sum(h**2)
def af_mae(h,truth):
    af=np.asarray((h@vp)/np.maximum(h@vc,1e-9)).ravel(); ok=np.isfinite(af)&np.isfinite(truth)
    return float(np.mean(np.abs(af[ok]-truth[ok])))
def poisson_ll(h,c,Knz):
    mu=np.asarray(h@Knz).ravel(); lam=c.sum()/max(mu.sum(),1e-12); r=np.maximum(lam*mu,1e-12)
    return float(np.sum(c*np.log(r)-r))

for pool,sub in [("n231_g0","cov10_n231_g0_s42_hotspots_p231_chr1"),("n50_g0","cov10_n50_g0_s42_hotspots_p231_chr1")]:
    cd=count_kmers_in_fasta([f"{ROOT}/benchmarks/p231/sims/{sub}/reads/r1.fq",f"{ROOT}/benchmarks/p231/sims/{sub}/reads/r2.fq"],list(ki),k=31,threads=8,hash_size="3G")
    counts=np.array([cd[k] for k in ki],dtype=np.float64)
    truth=pd.read_csv(f"{ROOT}/benchmarks/p231/sims/{sub}/recomb_truth_raw.tsv.gz",sep="\t",usecols=["truth_af"])["truth_af"].values.astype(float)
    pw=pd.read_csv(f"{ROOT}/benchmarks/p231/sims/{sub}/pool_weights.tsv",sep="\t"); present=set(pw[pw.weight>0].founder.astype(str))
    absent=np.array([f not in present for f in fo])
    idx=np.flatnonzero(kmer_block==7); nz=idx[counts[idx]>0]
    Knz=np.ascontiguousarray(np.asarray(Ksp[:,nz].todense(),np.float32)); c=counts[nz].astype(np.float32)
    kfw=np.asarray(Ksp[:,idx].sum(1)).ravel().astype(np.float32)
    print(f"\n######## {pool}  blk7 centromere (n_nz={len(nz):,}, true present={len(present)}) ########",flush=True)
    # baseline (uniform init)
    h0,_=solve_em(c,Knz,1.0,max_iter=300,tol=1e-8,normalize="per_founder",kfw=kfw)
    print(f"  baseline (uniform init)   eff_n={effn(h0):6.1f} AF_MAE={af_mae(h0,truth):.4f} mass_absent={h0[absent].sum():.3f} LL={poisson_ll(h0,c,Knz):.0f}",flush=True)
    # HARP mechanism 1: multi-random-start (symmetric Dirichlet alpha=1), keep best LL
    rng=np.random.default_rng(0); best=None; effs=[]
    for s in range(8):
        hi=rng.dirichlet(np.ones(F)).astype(np.float32)
        hs,_=solve_em(c,Knz,1.0,h_init=hi,max_iter=300,tol=1e-8,normalize="per_founder",kfw=kfw)
        ll=poisson_ll(hs,c,Knz); effs.append(effn(hs))
        if best is None or ll>best[0]: best=(ll,hs)
    hb=best[1]
    print(f"  multi-start best-of-8     eff_n={effn(hb):6.1f} AF_MAE={af_mae(hb,truth):.4f} mass_absent={hb[absent].sum():.3f} LL={best[0]:.0f}  (start eff_n range {min(effs):.0f}-{max(effs):.0f})",flush=True)
    # HARP mechanism 2: em_min_freq_cutoff on baseline
    for cut in (1e-3,5e-3):
        hc=np.where(h0>=cut,h0,0.0); hc=hc/hc.sum()
        print(f"  cutoff {cut:.0e}            eff_n={effn(hc):6.1f} AF_MAE={af_mae(hc,truth):.4f} mass_absent={hc[absent].sum():.3f}",flush=True)
print("DONE-HARP")
