#!/usr/bin/env python3
"""AUDIT: is the 0.0485 global noiseless floor genuine non-identifiability, or
the normalized-EM's estimator bias?

The production EM update (em_solver.solve_em) is h_f ∝ h_f Σ_k ω_k c_k K_fk/μ_k,
normalized to the simplex. That maximizes L(h)=Σ_k ω_k c_k log μ_k(h) s.t. Σh=1
(a multinomial/cross-entropy objective — NO Poisson −μ term). Its stationary
condition on the support is a_f := Σ_k ω_k K_fk = const ∀f. So h_true is a fixed
point ONLY IF every founder has the same ω-weighted k-mer count a_f.

Decisive test (memory-careful, matvecs only — no F×K temporaries):
  1. a_f = kmer_pa @ ω  → spread (CoV) across the true support.
  2. ||EMstep(h_true) − h_true|| on NOISELESS data (c_k=μ_k(h_true)). If ≫0,
     h_true is NOT a fixed point ⇒ the floor is (at least partly) estimator bias,
     NOT irreducible non-identifiability.
  3. ||PNMFstep(h_true) − h_true|| where PNMF divides by per-founder a_f
     (the Poisson-NMF update). Should be ~0 (h_true IS its fixed point).
  4. Run BOTH updates to convergence from uniform on noiseless data; report
     ||ĥ−h_true||. If PNMF → ~0 but production stays ~0.0485 ⇒ floor is the
     normalization choice (fixable), not collinearity.
"""
import sys, time
from pathlib import Path
import numpy as np, pandas as pd
from scipy import sparse

ROOT = Path("/global/scratch/users/tbellg/kmate")
PRE = ROOT / "benchmarks/p80/data/kmer_pa_p80_filt2"

t = time.time()
kp = sparse.load_npz(PRE / "kmer_pa_Chr1.kmer_pa.npz")
kmer_pa = np.asarray(kp.todense(), dtype=np.float32) if sparse.issparse(kp) \
    else np.asarray(kp, dtype=np.float32)
del kp
m = np.load(PRE / "kmer_pa_Chr1.meta.npz", allow_pickle=True)
founders = np.asarray(m["founders"]).astype(str)
bid = np.asarray(m["bubble_id"])
omega = (1.0 / np.bincount(bid)[bid]).astype(np.float32)
w = pd.read_csv(ROOT / "benchmarks/p80/sims/cov10_n50_g0_s42_hotspots_p80_chr1/pool_weights.tsv", sep="\t")
wmap = {str(f): float(x) for f, x in zip(w["founder"], w["weight"])}
h_true = np.array([wmap.get(str(f), 0.0) for f in founders], dtype=np.float64)
h_true /= h_true.sum()
supp = h_true > 0
F, K = kmer_pa.shape
print(f"load {time.time()-t:.0f}s  F={F} K={K:,} support={int(supp.sum())}", flush=True)

# 1. a_f = Σ_k ω_k K_fk  (matvec, no F×K temp)
a = kmer_pa @ omega                      # (F,)
a_s = a[supp]
print(f"\n[1] a_f = ω-weighted k-mer count per founder:")
print(f"    over support: mean={a_s.mean():.1f}  CoV={a_s.std()/a_s.mean():.3f}  "
      f"min={a_s.min():.0f}  max={a_s.max():.0f}  (max/min={a_s.max()/a_s.min():.2f})")
print(f"    ⇒ if CoV≈0 the normalized-EM is ~unbiased; large CoV ⇒ h_true not a fixed point")

mu = (h_true @ kmer_pa).astype(np.float64)   # noiseless counts, λ=1  (matvec)

def em_step(h):                          # production normalized update
    denom = np.maximum(h @ kmer_pa, 1e-12)
    cw = (omega * mu) / denom
    em = h * (kmer_pa @ cw)
    return em / em.sum()

def pnmf_step(h):                        # Poisson-NMF: divide by per-founder a_f
    denom = np.maximum(h @ kmer_pa, 1e-12)
    cw = (omega * mu) / denom
    num = h * (kmer_pa @ cw)
    hn = num / np.maximum(a, 1e-12)
    return hn / hn.sum()

# 2,3. one step from h_true
r_em = np.linalg.norm(em_step(h_true) - h_true)
r_pn = np.linalg.norm(pnmf_step(h_true) - h_true)
print(f"\n[2] ||EMstep(h_true)-h_true||   = {r_em:.4f}   (production normalized)")
print(f"[3] ||PNMFstep(h_true)-h_true|| = {r_pn:.4f}   (per-founder normalized)")
print(f"    ⇒ EM step ≫0 means h_true is NOT a fixed point of the production EM")

# 4. converge both from uniform on noiseless data
def run(step, n=4000, tol=1e-12):
    h = np.full(F, 1.0 / F)
    for it in range(n):
        hn = step(h)
        d = np.linalg.norm(hn - h); h = hn
        if d < tol:
            break
    return h, it + 1, d

t = time.time()
h_em, it_em, d_em = run(em_step)
h_pn, it_pn, d_pn = run(pnmf_step)
e_em = np.linalg.norm(h_em - h_true)
e_pn = np.linalg.norm(h_pn - h_true)
print(f"\n[4] noiseless convergence from uniform ({time.time()-t:.0f}s):")
print(f"    production EM : ||ĥ-h_true||={e_em:.4f}  (it {it_em}, dh {d_em:.1e}, "
      f"mass_off_supp={h_em[~supp].sum():.4f})")
print(f"    Poisson-NMF   : ||ĥ-h_true||={e_pn:.4f}  (it {it_pn}, dh {d_pn:.1e}, "
      f"mass_off_supp={h_pn[~supp].sum():.4f})")
print(f"\nVERDICT: if PNMF≈0 and production≈0.0485 ⇒ the global floor is the EM "
      f"NORMALIZATION (estimator bias, fixable), not irreducible non-identifiability.")
print(f"         if BOTH ≈0.0485 ⇒ genuine collinear non-identifiability (doc correct).")
