#!/usr/bin/env python3
"""Validate the rebuilt p231 cn_full against the production cn_full_231_v3qc_v3_filt2.
If cn_full is representation-invariant (consensus-derived from the same panel), the
two should be nearly identical. Reports: founder-axis match, k-mer set overlap,
per-founder k-mer-count correlation, and carrier-pattern agreement on shared k-mers."""
import numpy as np
from scipy.sparse import load_npz
ROOT="/carnegie/nobackup/scratch/tbellagio/hapfire_sv"
P231=f"{ROOT}/control_p231/data/cn_full_p231_filt2"
PROD=f"{ROOT}/poolfreq/data/cn_full_231_v3qc_v3_filt2"

def load(d):
    cn=load_npz(f"{d}/cn_Chr1.cn.npz").tocsc()
    m=np.load(f"{d}/cn_Chr1.meta.npz",allow_pickle=True)
    return cn, np.asarray(m["founders"]).astype(str), np.asarray(m["kmer_index"]).astype(str)

cnA,foA,kmA=load(P231); cnB,foB,kmB=load(PROD)
print(f"p231 rebuilt : F={len(foA)} K={cnA.shape[1]:,} nnz={cnA.nnz:,}")
print(f"production   : F={len(foB)} K={cnB.shape[1]:,} nnz={cnB.nnz:,}")
print(f"\nfounders identical order: {np.array_equal(foA,foB)}")

sA,sB=set(kmA),set(kmB); inter=sA&sB
print(f"\nk-mer sets: p231={len(sA):,} prod={len(sB):,} shared={len(inter):,} "
      f"({100*len(inter)/max(len(sA),len(sB)):.2f}% of larger)")
print(f"  p231-only={len(sA-sB):,}  prod-only={len(sB-sA):,}")

# per-founder k-mer count (carrier count) correlation on the FULL matrices
kcA=np.asarray(cnA.sum(axis=1)).flatten(); kcB=np.asarray(cnB.sum(axis=1)).flatten()
if np.array_equal(foA,foB):
    r=np.corrcoef(kcA,kcB)[0,1]
    print(f"\nper-founder total k-mer count: corr(p231,prod)={r:.6f}")
    print(f"  median |Δ|={np.median(np.abs(kcA-kcB)):.0f}  median count prod={np.median(kcB):.0f} "
          f"(rel {100*np.median(np.abs(kcA-kcB))/np.median(kcB):.2f}%)")

# carrier-pattern agreement on shared k-mers
if np.array_equal(foA,foB) and inter:
    iA={k:i for i,k in enumerate(kmA)}; iB={k:i for i,k in enumerate(kmB)}
    shared=np.array(sorted(inter))
    samp=shared if len(shared)<=200000 else np.random.default_rng(0).choice(shared,200000,replace=False)
    colA=cnA[:,[iA[k] for k in samp]].toarray().astype(bool)
    colB=cnB[:,[iB[k] for k in samp]].toarray().astype(bool)
    agree=(colA==colB).mean()
    identical_cols=(colA==colB).all(axis=0).mean()
    print(f"\ncarrier-pattern agreement on {len(samp):,} shared k-mers:")
    print(f"  per-cell agreement: {100*agree:.4f}%")
    print(f"  fully-identical carrier columns: {100*identical_cols:.4f}%")
print("\n=> If overlap ~100% and carrier agreement ~100%, cn_full is representation-"
      "invariant; the rebuilt and production cn_full are equivalent (consistency proven).")
