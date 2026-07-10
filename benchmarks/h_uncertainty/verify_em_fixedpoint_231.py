#!/usr/bin/env python3
"""Caveat check on the REAL 231 heterogeneous panel (78 cactus + 153 PG).

On p80 (homogeneous) ω=1/m_b nearly equalized the per-founder ω-weighted k-mer
count a_f (CoV 1.6%), so the normalized-EM's estimator bias was tiny (0.0023 vs
the 0.0485 non-identifiability floor). The 231 panel has the cactus/PG k-mer
richness asymmetry that ω only PARTIALLY cancels, so a_f spread — and thus the
estimator-bias slice — should be larger. Measure it.

Matvecs only (no F×K temporaries). 231 kmer_pa Chr1 dense ≈ 10 GB.
"""
import sys, time
import numpy as np, pandas as pd
from scipy import sparse
from pathlib import Path

ROOT = Path("/global/scratch/users/tbellg/kmate")
PRE = ROOT / "data/kmer_pa_231_arch3_filt2inv"
POOL = "cov10_n50_g0_s42_hotspots_p231_chr1"

t = time.time()
kp = sparse.load_npz(PRE / "kmer_pa_Chr1.kmer_pa.npz")
K = np.asarray(kp.todense(), dtype=np.float32) if sparse.issparse(kp) else np.asarray(kp, np.float32)
del kp
m = np.load(PRE / "kmer_pa_Chr1.meta.npz", allow_pickle=True)
fnd = np.asarray(m["founders"]).astype(str)
bid = np.asarray(m["bubble_id"])
omega = (1.0 / np.bincount(bid)[bid]).astype(np.float32)
F, Kn = K.shape
print(f"load {time.time()-t:.0f}s  F={F}  K={Kn:,}", flush=True)

a = (K @ omega).astype(np.float64)        # a_f = Σ_k ω_k K_fk  (h-independent)
print(f"\n[a_f over ALL {F} founders]  CoV={a.std()/a.mean():.3f}  "
      f"min={a.min():.0f}  max={a.max():.0f}  max/min={a.max()/a.min():.2f}", flush=True)
print(f"   (p80 reference: CoV 0.016, max/min 1.08)", flush=True)

# one-step residual at a real p231 pool truth
try:
    w = pd.read_csv(ROOT / f"benchmarks/p231/sims/{POOL}/pool_weights.tsv", sep="\t")
    wm = {str(f): float(x) for f, x in zip(w["founder"], w["weight"])}
    h = np.array([wm.get(str(f), 0.0) for f in fnd], dtype=np.float32); h /= h.sum()
    supp = h > 0
    print(f"\npool={POOL}  support={int(supp.sum())}", flush=True)
    asu = a[supp]
    print(f"[a_f over the {int(supp.sum())} support founders] CoV={asu.std()/asu.mean():.3f}  "
          f"max/min={asu.max()/asu.min():.2f}", flush=True)
    mu = (h @ K).astype(np.float32)
    den = np.maximum(h @ K, np.float32(1e-12)); cw = (omega * mu) / den
    e = h * (K @ cw); h1 = e / e.sum()
    r_em = np.linalg.norm(h1.astype(np.float64) - h.astype(np.float64))
    hp = (h.astype(np.float64) * a); hp /= hp.sum()
    print(f"[one-step] ||EMstep(h_true)-h_true|| = {r_em:.4f}  "
          f"(p80 was 0.0023);  predicted h∝h_true·a_f = {np.linalg.norm(hp-h.astype(np.float64)):.4f}", flush=True)
    print("\nINTERPRET: larger a_f CoV / one-step residual than p80 ⇒ on the real panel the "
          "estimator-NORMALIZATION bias is a bigger slice of the global floor (still added to, "
          "not replacing, the non-identifiability manifold).", flush=True)
except FileNotFoundError:
    print(f"(no weights for {POOL}; a_f CoV over all founders is the key number)", flush=True)
