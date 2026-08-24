#!/usr/bin/env python
"""Definitive mechanism: does kMate under-call exactly the founders with the LEAST PRIVATE
k-mers on the PRODUCTION panel (data/kmer_pa_231_arch3_filt2inv)?

For each founder, per chromosome:
  mean_ac  = avg # founders sharing each of its k-mers = (K @ ac)/Kf   (HIGH = non-private)
  priv     = # k-mers it alone carries (ac==1)                          (LOW  = non-identifiable)
  Kf       = its total k-mer count
Then correlate with the kMate seed-mix p0 (genome-wide) and per-chrom h. Prediction:
absorbed founders (low p0) have HIGH mean_ac / LOW priv; and a founder is absorbed on the
chromosomes where its k-mers are least private (per-(founder,chrom) test).

Reuses diag_9761.py's loading. Env: kmate. Heavy (5 x ~500MB k-mer matrices) -> background.
-> analysis/grenenet_gea/qc/seedmix_validation/kmer_identifiability.npz + prints.
"""
import numpy as np, glob, os, sys, time
from scipy.sparse import load_npz
from scipy import stats
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import lib

PANEL="data/kmer_pa_231_arch3_filt2inv/kmer_pa"
CH=[f"Chr{i}" for i in range(1,6)]; N=231
OUT="analysis/grenenet_gea/seedmix_validation"; os.makedirs(OUT,exist_ok=True)
t0=time.time()

# founder order of the k-mer panel (cheap: npz loads one key lazily)
kf_founders=np.load(f"{PANEL}_Chr1.meta.npz",allow_pickle=True)["founders"].astype(str)

# kMate seed-mix h aligned to the k-mer panel founder order
reps=sorted({os.path.basename(p).split("_Chr")[0] for p in glob.glob(f"{lib.SEEDMIX}/*_Chr1.h_per_chrom.npz")})
z0=np.load(f"{lib.SEEDMIX}/{reps[0]}_Chr1.h_per_chrom.npz",allow_pickle=True)
sm_founders=z0["founders"].astype(str)
pos_in_sm={f:i for i,f in enumerate(sm_founders)}
sel=np.array([pos_in_sm[f] for f in kf_founders])     # reorder seedmix -> k-mer panel order
Hk=np.zeros((len(reps),5,N))
for ri,s in enumerate(reps):
    for ci,c in enumerate(CH):
        Hk[ri,ci]=np.load(f"{lib.SEEDMIX}/{s}_{c}.h_per_chrom.npz",allow_pickle=True)[c].astype(float)[sel]
kmate=Hk.mean(1).mean(0); h_chrom=Hk.mean(0)          # (N,), (5,N) in k-mer panel order

mean_ac=np.zeros((5,N)); priv=np.zeros((5,N)); Kf=np.zeros((5,N))
for ci,c in enumerate(CH):
    K=load_npz(f"{PANEL}_{c}.kmer_pa.npz").astype(np.float32).tocsr()
    ac=np.asarray(K.sum(0)).ravel(); kf=np.asarray(K.sum(1)).ravel()
    mean_ac[ci]=(K @ ac)/np.maximum(kf,1)
    priv[ci]=np.asarray(K[:,ac==1].sum(1)).ravel()
    Kf[ci]=kf
    print(f"{c}: K={K.shape} nnz={K.nnz:,} kmers={K.shape[1]:,} singletons(ac==1)={int((ac==1).sum()):,} "
          f"({time.time()-t0:.0f}s)",flush=True)
    del K

# genome-wide per founder: k-mer-weighted mean commonness + total private
mac_gw=(mean_ac*Kf).sum(0)/Kf.sum(0); priv_gw=priv.sum(0); Kf_gw=Kf.sum(0)
low=kmate<=1e-3; mid=(kmate>1e-3)&(kmate<1/462); norm=kmate>=1/462

print("\n=== genome-wide k-mer identifiability vs kMate under-calling ===")
print(f"{'group':26} {'n':>4} {'med_mean_ac':>12} {'med_priv_kmers':>15} {'med_Kf':>10}")
for lab,m in [("kMate p0<=1e-3 (dropped)",low),("kMate 1e-3..1/462",mid),("kMate>=1/462 normal",norm)]:
    print(f"{lab:26} {m.sum():4d} {np.median(mac_gw[m]):12.1f} {np.median(priv_gw[m]):15.0f} {np.median(Kf_gw[m]):10.0f}")
print(f"\nSpearman(kMate p0, mean_ac)     = {stats.spearmanr(kmate,mac_gw).correlation:+.3f}  (POS => absorbed founders have MORE-common k-mers)")
print(f"Spearman(kMate p0, priv_kmers)  = {stats.spearmanr(kmate,priv_gw).correlation:+.3f}  (POS => absorbed founders have FEWER private k-mers)")
print(f"Spearman(kMate p0, Kf)          = {stats.spearmanr(kmate,Kf_gw).correlation:+.3f}")

print("\n=== the 19 low founders (k-mer panel order): mechanism ===")
print(f"{'founder':>8} {'kMate_p0':>9} {'mean_ac':>8} {'priv_kmers':>11} {'rank_mac':>9}  (rank 1 = most-common k-mers)")
mac_rank=stats.rankdata(-mac_gw)   # 1 = highest mean_ac
for i in np.where(low)[0][np.argsort(kmate[low])]:
    print(f"{kf_founders[i]:>8} {kmate[i]:.2e} {mac_gw[i]:8.1f} {priv_gw[i]:11.0f} {int(mac_rank[i]):9d}/231")

print("\n=== per-(founder,chrom): is a founder absorbed on the chroms where its k-mers are least private? ===")
ch_h=h_chrom.ravel(); ch_mac=mean_ac.ravel(); ch_priv=priv.ravel()
print(f"Spearman(per-chrom h, per-chrom mean_ac) all cells   = {stats.spearmanr(ch_h,ch_mac).correlation:+.3f}")
print(f"Spearman(per-chrom h, per-chrom priv)    all cells   = {stats.spearmanr(ch_h,ch_priv).correlation:+.3f}")
lm=np.repeat(low.reshape(1,-1),5,0).ravel()
print(f"Spearman(per-chrom h, per-chrom mean_ac) low founders= {stats.spearmanr(ch_h[lm],ch_mac[lm]).correlation:+.3f}")

np.savez(f"{OUT}/kmer_identifiability.npz",founders=kf_founders,kmate=kmate,h_chrom=h_chrom,
         mean_ac=mean_ac,priv=priv,Kf=Kf,mac_gw=mac_gw,priv_gw=priv_gw,low=low)
print(f"\nwrote {OUT}/kmer_identifiability.npz ({time.time()-t0:.0f}s)")
