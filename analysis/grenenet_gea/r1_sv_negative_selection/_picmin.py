#!/usr/bin/env python
"""REAL PicMin (Booker et al. 2024) applied to GrENE-net (user 2026-07-03). Lineages = SITES; per-
lineage per-locus statistic = the drift-controlled parallelism |z| (among-plot mean-slope z at that
site, from parallelism.npz). PicMin:
  1. per site, convert each locus |z| to an empirical p vs the SNP |z| distribution in its p0-bin
     (frequency-controlled; small p = strong parallel selection at that site).
  2. per locus, sort its N site p-values; the a-th order statistic ~ Beta(a, N-a+1) under H0.
     stat = min_a F_Beta(p_(a); a, N-a+1)  (the most-significant order = repeated across >=a sites).
  3. calibrate that min-over-orders statistic against a uniform null -> per-locus PicMin p-value; BH-FDR.
Then the SV-vs-SNP contrast: fraction of SVs that are repeated-adaptation loci vs frequency-matched SNPs.

Env: kmate. Reads analysis/grenenet_gea/sv_adaptive/parallelism.npz.
"""
import os, sys
import numpy as np
from scipy.stats import beta
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import lib

G = f"{lib.GEA}/sv_adaptive"
z = np.load(f"{G}/parallelism.npz")
zsv, zin, zsn = z["z_sv"], z["z_indel"], z["z_snp"]
p0sv, p0in, p0sn = z["p0_sv"], z["p0_indel"], z["p0_snp"]
isdel = z["isdel_sv"]
N = zsv.shape[1]; rng = np.random.default_rng(0)
print(f"lineages(sites) N={N}; SV={zsv.shape[0]} indel={zin.shape[0]} SNP={zsn.shape[0]}")

p0q = np.quantile(p0sn, np.linspace(0, 1, 11)); p0q[-1] += 1e-9
binf = lambda p: np.clip(np.digitize(p, p0q[1:-1]), 0, 9)
bsv, bin_, bsn = binf(p0sv), binf(p0in), binf(p0sn)


def emp_p(zmat, bins):
    """per site & p0-bin, empirical p of each locus vs the SNP |z| reference (small p = high |z|)."""
    P = np.empty_like(zmat, dtype=np.float64)
    for s in range(N):
        zs_ref = zsn[:, s]
        for b in range(10):
            ref = np.sort(zs_ref[bsn == b]); m = bins == b
            if ref.size == 0:
                P[m, s] = 0.5; continue
            cnt = ref.size - np.searchsorted(ref, zmat[m, s], side="left")   # #ref >= val
            P[m, s] = (1 + cnt) / (ref.size + 1)
    return P


def picmin_stat(P):
    Ps = np.sort(P, axis=1); a = np.arange(1, N + 1)
    F = beta.cdf(Ps, a[None, :], (N - a + 1)[None, :])
    return F.min(axis=1)                      # small = repeated across the best number of sites


print("empirical per-site p-values ...")
Psv, Pin, Psn = emp_p(zsv, bsv), emp_p(zin, bin_), emp_p(zsn, bsn)
print("PicMin order statistics + uniform-null calibration ...")
stat_sv, stat_in, stat_sn = picmin_stat(Psv), picmin_stat(Pin), picmin_stat(Psn)
nulls = np.sort(picmin_stat(rng.uniform(size=(300000, N))))
pv = lambda st: (1 + np.searchsorted(nulls, st, side="right")) / (nulls.size + 1)
p_sv, p_in, p_sn = pv(stat_sv), pv(stat_in), pv(stat_sn)


def bh(p, q=0.1):
    o = np.argsort(p); n = len(p); thr = q * (np.arange(1, n + 1)) / n
    below = p[o] <= thr; k = np.where(below)[0].max() + 1 if below.any() else 0
    sig = np.zeros(n, bool); sig[o[:k]] = True; return sig


def matchfrac(sig_ref, p0_ref, p0_tgt, val="sig"):
    """frequency-match ref(SNP) to target p0 dist, return matched fraction significant."""
    e = np.quantile(np.concatenate([p0_tgt, p0_ref]), np.linspace(0, 1, 26)); e[-1] += 1e-9
    tb = np.digitize(p0_tgt, e[1:-1]); rb = np.digitize(p0_ref, e[1:-1]); acc = []
    for k in range(25):
        pool = sig_ref[rb == k]; nn = int(round((tb == k).mean() * 20000))
        if nn and pool.size: acc.append(rng.choice(pool, nn, replace=True))
    return np.concatenate(acc).mean()


for q in (0.10, 0.05):
    sig_sv = bh(p_sv, q); sig_sn = bh(p_sn, q); sig_in = bh(p_in, q)
    m_sn = matchfrac(sig_sn, p0sn, p0sv)
    print(f"\n=== PicMin repeated-adaptation loci (BH q<{q}) ===")
    print(f"  SV:    {sig_sv.mean()*100:.2f}% ({sig_sv.sum()}/{len(sig_sv)})   "
          f"vs freq-matched SNP {m_sn*100:.2f}%   fold={sig_sv.mean()/max(m_sn,1e-9):.2f}")
    print(f"  indel: {sig_in.mean()*100:.2f}%   SNP(raw): {sig_sn.mean()*100:.2f}%")
    print(f"  SV insertions {sig_sv[~isdel].mean()*100:.2f}%  vs deletions {sig_sv[isdel].mean()*100:.2f}%")

np.savez(f"{G}/picmin.npz", p_sv=p_sv, p_indel=p_in, p_snp=p_sn, p0_sv=p0sv, p0_snp=p0sn,
         p0_indel=p0in, isdel_sv=isdel, stat_sv=stat_sv, stat_snp=stat_sn)
print("\n[wrote] picmin.npz")
