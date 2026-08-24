"""
POC: does LOCAL-window EM + identifiability-weighted aggregation recover the
absorbed founders that GLOBAL (chrom-wide) EM collapses?

Ideal noiseless test on Chr1: c = K^T @ h_true, h_true = uniform 1/231.
If global EM on ideal counts collapses the absorbed founders while local-window
EM recovers them, that isolates the estimator geometry (no read noise, no hapFIRE).

Copies the EM math from src/kmate/em_solver.solve_em (sparse-aware). No edits to src.
"""
import numpy as np, scipy.sparse as sp, zipfile, io, sys, time

DDIR='/global/scratch/users/tbellg/kmate/data/kmer_pa_231_arch3_filt2inv'
OUT='/global/scratch/users/tbellg/kmate/analysis/grenenet_gea/qc/seedmix_validation/fix_local'

def em(K, counts, h_init=None, max_iter=300, tol=1e-7):
    """K: F x Ncol sparse CSR. counts: (Ncol,). Returns h (F,). MLE EM (no prior)."""
    F=K.shape[0]
    counts=counts.astype(np.float32)
    h=np.full(F,1.0/F,np.float32) if h_init is None else h_init.astype(np.float32).copy()
    total_c=counts.sum()
    Kt=K.T.tocsr()  # for denom = Kt @ h  (Ncol,)
    if total_c<=0: return h.astype(np.float64),0
    for it in range(max_iter):
        denom=np.maximum(Kt.dot(h),np.float32(1e-7))
        cw=counts/denom
        em_term=h*K.dot(cw)
        h_new=em_term/max(em_term.sum(),np.float32(1e-12))
        d=np.linalg.norm(h_new-h); h=h_new
        if d<tol: break
    return h.astype(np.float64),it+1

def main():
    t0=time.time()
    print("loading Chr1 panel...",flush=True)
    K=sp.load_npz(f'{DDIR}/kmer_pa_Chr1.kmer_pa.npz').astype(np.float32)  # 231 x Ncol
    F,N=K.shape
    z=zipfile.ZipFile(f'{DDIR}/kmer_pa_Chr1.meta.npz')
    ld=lambda n: np.load(io.BytesIO(z.read(n)),allow_pickle=True)
    founders=ld('founders.npy').astype(str)
    bid=ld('bubble_id.npy'); bstart=ld('bubble_start.npy')
    pos=bstart[bid].astype(np.int64)  # per-kmer genomic position (cols sorted ascending)
    target=1.0/F
    fmap={f:i for i,f in enumerate(founders)}
    absorbed=['9977','9985','9941','10013','9507','9978','9748','9966']
    aidx=[fmap[a] for a in absorbed]
    print(f"K={K.shape} nnz={K.nnz} load {time.time()-t0:.0f}s",flush=True)

    # ideal noiseless counts from uniform truth
    h_true=np.full(F,target)
    c=np.asarray(K.T.dot(h_true)).ravel().astype(np.float32)

    # ---- GLOBAL EM (chrom-wide) ----
    tg=time.time()
    h_glob,itg=em(K,c)
    print(f"GLOBAL EM {itg} iters {time.time()-tg:.0f}s",flush=True)

    # ---- LOCAL window EM, multiple window sizes ----
    ac_full=np.asarray(K.sum(0)).ravel()  # per-kmer allele count
    results={'founders':founders,'target':target,'h_true':h_true,
             'h_global':h_glob,'absorbed':np.array(absorbed)}
    span=pos.max()-pos.min()
    for nwin in [20,100,500]:
        edges=np.linspace(pos.min(),pos.max()+1,nwin+1)
        # column boundaries (pos sorted ascending)
        bnds=np.searchsorted(pos,edges)
        Hwin=np.full((nwin,F),np.nan)
        # per-window per-founder identifiability weight: sum over window kmers the
        # founder carries of 1/ac  (distinctive-kmer mass); higher = better resolved
        Wid=np.zeros((nwin,F)); Wkf=np.zeros((nwin,F))
        inv_ac=1.0/ac_full
        for w in range(nwin):
            a,b=bnds[w],bnds[w+1]
            if b-a<10:  # too few kmers
                continue
            Kw=K[:,a:b]
            cw=c[a:b]
            hw,_=em(Kw,cw,max_iter=200)
            Hwin[w]=hw
            # identifiability: for each founder, sum of 1/ac over carried kmers in window
            invw=inv_ac[a:b]
            Wid[w]=np.asarray(Kw.dot(invw)).ravel()   # F-vector
            Wkf[w]=np.asarray(Kw.sum(1)).ravel()       # k-mer count per founder in window
        valid=~np.isnan(Hwin[:,0])
        Hw=Hwin[valid]; Wid_v=Wid[valid]; Wkf_v=Wkf[valid]
        # aggregation schemes (per founder, over windows)
        agg_unif=np.nanmean(Hw,0)
        def wmean(H,W):
            Wn=W/np.maximum(W.sum(0,keepdims=True),1e-12)
            return (Wn*H).sum(0)
        agg_kf=wmean(Hw,Wkf_v)     # ~ what global effectively does
        agg_id=wmean(Hw,Wid_v)     # identifiability-weighted
        results[f'nwin{nwin}_Hwin']=Hwin
        results[f'nwin{nwin}_agg_unif']=agg_unif
        results[f'nwin{nwin}_agg_kf']=agg_kf
        results[f'nwin{nwin}_agg_id']=agg_id
        print(f"\n=== nwin={nwin} ({valid.sum()} valid windows, ~{span/nwin/1e3:.0f}kb each) ===")
        def report(name,h):
            frac=h/target
            nbad=int((h<0.5*target).sum())
            l1=np.abs(h-target).sum()
            absorbed_frac=frac[aidx]
            print(f"  {name:10s} nbad(<0.5x)={nbad:3d}  L1={l1:.4f}  absorbed x/target:"
                  +" ".join(f"{x:.2f}" for x in absorbed_frac))
        report('GLOBAL',h_glob)
        report('unifmean',agg_unif)
        report('kf-wt',agg_kf)
        report('id-wt',agg_id)
    np.savez_compressed(f'{OUT}/poc_local_em_chr1.npz',**results)
    print(f"\nsaved. total {time.time()-t0:.0f}s")

if __name__=='__main__':
    main()
