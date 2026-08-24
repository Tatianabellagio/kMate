#!/usr/bin/env python
"""Why does kMate under-call ~19-49 seed-mix founders (vs hapFIRE ~1/231)?

Mechanism test: kMate's mixture EM resolves a founder through its PRIVATE / distinguishing
variants (proxy for private k-mers). A founder with few private markers sits on a flat
likelihood ridge and the non-negative EM regularizes it toward 0. This is a MARGINAL
identifiability property (how much unique info the founder carries), not pairwise-twin
similarity (which we already ruled out on the variant GRM).

Computes, per founder, genome-wide and per-chrom:
  * priv1/2/3  = # well-called variants where the founder is one of <=1/2/3 alt-carriers
                 (private / near-private markers -> the EM's leverage to detect it)
  * n_carry    = total alt-carrier variants (baseline)
and tests these against the kMate seed-mix p0 and per-chrom h.

Also emits the per-chrom x per-sample h matrix for the low founders ("are they always low?").
Env: kmate. Compute node. -> analysis/grenenet_selection/qc/seedmix_validation/identifiability.npz + prints.
"""
import numpy as np, glob, os, sys
import scipy.sparse as sp
from scipy import stats
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import lib

CH=[f"Chr{i}" for i in range(1,6)]; N=231
OUT="analysis/grenenet_selection/seedmix_validation"; os.makedirs(OUT,exist_ok=True)

# ---- kMate seed-mix h: reps x chroms x founders ----
reps=sorted({os.path.basename(p).split("_Chr")[0] for p in glob.glob(f"{lib.SEEDMIX}/*_Chr1.h_per_chrom.npz")})
founders=None; Hk=np.zeros((len(reps),5,N))
for ri,s in enumerate(reps):
    for ci,c in enumerate(CH):
        z=np.load(f"{lib.SEEDMIX}/{s}_{c}.h_per_chrom.npz",allow_pickle=True)
        if founders is None: founders=z["founders"].astype(str)
        Hk[ri,ci]=z[c].astype(float)
kmate=Hk.mean(1).mean(0)               # (N,) genome-wide p0
h_chrom=Hk.mean(0)                     # (5,N) rep-averaged per-chrom h

# ---- private-marker counts from var_pa, per chrom ----
panel_f=np.load(f"panel/arch3/chr1/var_pa_231_arch3_chr1.meta.npz",allow_pickle=True)["founders"].astype(str)
order=np.array([list(panel_f).index(f) for f in founders])   # align panel -> seedmix founder order
priv1=np.zeros((5,N)); priv2=np.zeros((5,N)); priv3=np.zeros((5,N)); ncarry=np.zeros((5,N))
for ci,c in enumerate(CH):
    cl=c.lower(); base=f"panel/arch3/{cl}/var_pa_231_arch3_{cl}"
    vp=sp.load_npz(f"{base}.var_pa.npz").tocsr()[order].tocsc()
    vc=sp.load_npz(f"{base}.var_called.npz").tocsr()[order].tocsc()
    n_alt=np.asarray(vp.sum(0)).ravel(); n_cal=np.asarray(vc.sum(0)).ravel()
    wc=n_cal>=0.9*N                       # well-called variants only
    ncarry[ci]=np.asarray(vp[:,wc].sum(1)).ravel()
    for k,arr in ((1,priv1),(2,priv2),(3,priv3)):
        cols=wc&(n_alt<=k)&(n_alt>=1)
        arr[ci]=np.asarray(vp[:,cols].sum(1)).ravel()
    print(f"{c}: variants={vp.shape[1]:,} wellcalled={int(wc.sum()):,} "
          f"singletons(n_alt==1)={int(((n_alt==1)&wc).sum()):,}",flush=True)

P1=priv1.sum(0); P2=priv2.sum(0); P3=priv3.sum(0); NC=ncarry.sum(0)   # genome-wide per founder
low=kmate<=1e-3; mid=(kmate>1e-3)&(kmate<1/462); norm=kmate>=1/462

print("\n=== genome-wide private-marker count vs kMate under-calling ===")
print(f"{'group':26} {'n':>4} {'med_priv1':>10} {'med_priv2':>10} {'med_priv3':>10} {'med_ncarry':>11}")
for lab,m in [("kMate p0<=1e-3 (dropped)",low),("kMate p0 1e-3..1/462",mid),("kMate p0>=1/462 normal",norm)]:
    print(f"{lab:26} {m.sum():4d} {np.median(P1[m]):10.0f} {np.median(P2[m]):10.0f} "
          f"{np.median(P3[m]):10.0f} {np.median(NC[m]):11.0f}")
for nm,v in (("priv1",P1),("priv2",P2),("priv3",P3),("ncarry",NC)):
    print(f"Spearman(kMate p0, {nm:6}) = {stats.spearmanr(kmate,v).correlation:+.3f}")

print("\n=== the 19 low founders: genome-wide private markers + is it always low per-chrom? ===")
print(f"{'founder':>8} {'kMate_p0':>9} {'priv1':>6} {'priv2':>6} | per-chrom h  (Chr1..5)   | #chr with h<1e-3")
for i in np.where(low)[0][np.argsort(kmate[low])]:
    hc=h_chrom[:,i]; nlow=int((hc<1e-3).sum())
    print(f"{founders[i]:>8} {kmate[i]:.2e} {P1[i]:6.0f} {P2[i]:6.0f} | "
          + " ".join(f"{v:.1e}" for v in hc) + f"  | {nlow}/5")

# ---- per (founder,chrom) cell test: does per-chrom h track per-chrom private count? ----
Pc=priv1  # (5,N) per-chrom private singletons
cell_h=h_chrom.ravel(); cell_p=Pc.ravel()
print(f"\n=== per-(founder,chrom) cells (n={cell_h.size}) ===")
print(f"Spearman(per-chrom h, per-chrom priv1) all cells = {stats.spearmanr(cell_h,cell_p).correlation:+.3f}")
# within the low founders only
lm=np.repeat(low,1).reshape(1,-1).repeat(5,0).ravel()
print(f"Spearman within low-founder cells               = {stats.spearmanr(cell_h[lm],cell_p[lm]).correlation:+.3f}")

np.savez(f"{OUT}/identifiability.npz",founders=founders,kmate=kmate,h_chrom=h_chrom,
         priv1=P1,priv2=P2,priv3=P3,ncarry=NC,priv1_chrom=priv1,low=low)
print(f"\nwrote {OUT}/identifiability.npz")
