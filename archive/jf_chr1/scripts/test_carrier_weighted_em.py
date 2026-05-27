"""Test: weight each k-mer by its carrier count (ac_k) in the EM.

Rationale: standard EM amplifies low-carrier-count k-mers' updates by 1/ac_k
(via the denominator). Weighting by ac_k cancels this amplification, putting
every k-mer on equal-evidence footing — analogous to freqk's per-allele
normalization but applied per k-mer.

Implementation: replace c[k] with c[k] * ac_k before standard EM.
"""
import numpy as np
import pandas as pd
import scipy.sparse as sp
from pathlib import Path
import time

ROOT = Path('/carnegie/nobackup/scratch/tbellagio/hapfire_sv')


def em_sparse(counts, cn_sparse, max_iter=200, tol=1e-7, h_init=None):
    F, K = cn_sparse.shape
    counts = counts.astype(np.float32)
    h = np.full(F, 1.0/F, dtype=np.float32) if h_init is None else h_init.astype(np.float32)
    cn = cn_sparse.tocsr()
    cn_T = cn.T.tocsr()
    for it in range(max_iter):
        denom = cn_T @ h
        np.maximum(denom, np.float32(1e-7), out=denom)
        cw = counts / denom
        em_term = h * (cn @ cw)
        total_c = counts.sum()
        h_new = em_term / max(float(total_c), 1e-12)
        h_new = h_new / h_new.sum()
        delta = float(np.linalg.norm(h_new - h))
        h = h_new
        if delta < tol:
            break
    return h.astype(np.float64), {'iterations': it+1, 'final_delta': delta, 'converged': delta<tol}


def report(name, h, h_truth, is_cactus, iters=None):
    l1 = float(np.abs(h - h_truth).sum())
    l2 = float(((h - h_truth)**2).sum()**0.5)
    r = float(np.corrcoef(h, h_truth)[0,1]) if h.std()>0 else float('nan')
    return dict(
        name=name, iters=iters, L1=l1, L2=l2, r=r,
        cactus_mean=float(h[is_cactus].mean()),
        PG_mean=float(h[~is_cactus].mean()),
    )


# Load (same as v2)
print('Loading...', flush=True)
t0 = time.time()
cn_full = sp.load_npz(ROOT/'poolfreq/data/cn_full_231_v3/cn_Chr1.cn.npz').tocsr().astype(np.float32)
m1 = np.load(ROOT/'poolfreq/data/cn_full_231_v3/cn_Chr1.meta.npz', allow_pickle=True)
founders = list(np.asarray(m1['founders']).astype(str))
counts_all = np.load(ROOT/'poolfreq/data/seedmix_S1_counts_genomewide.npz')['counts']
n_chr1 = cn_full.shape[1]
counts_chr1 = counts_all[:n_chr1].astype(np.float32)

nz = counts_chr1 > 0
c = counts_chr1[nz]
cn = cn_full[:, nz].tocsr()
ac_k = np.asarray(cn.sum(axis=0)).ravel().astype(np.float32)
print(f'cn={cn.shape}, nnz={cn.nnz:,}, n_singletons={int((ac_k==1).sum()):,}  ({time.time()-t0:.0f}s)', flush=True)

cactus_ids = set((ROOT/'pangenie_genotyping/data/merged/cactus_overlap_80.txt').read_text().split())
is_cactus = np.array([f in cactus_ids for f in founders])

recipe = pd.read_csv(ROOT/'data/seedmix_recipe_normalized.tsv', sep='\t', dtype={'ID':str})
h_truth = np.array([dict(zip(recipe['ID'], recipe['seed_prop'])).get(f, 1/231) for f in founders])

h_v3_saved = np.load(ROOT/'poolfreq/results/seedmix_231_v3_perchrom/SEEDMIX_S1.h_per_chrom.npz', allow_pickle=True)['Chr1']

# === Tests ===
results = []

# A) baseline
print('\nA) baseline EM ...', flush=True)
t0=time.time()
hA, infoA = em_sparse(c, cn, max_iter=200)
print(f'  done ({time.time()-t0:.0f}s, iters={infoA["iterations"]}, delta={infoA["final_delta"]:.2e})', flush=True)
results.append(report('A_baseline', hA, h_truth, is_cactus, iters=infoA['iterations']))

# H) carrier-count weighted: c_new = c * ac_k
print('\nH) carrier-count-weighted EM (c_new = c * ac_k) ...', flush=True)
c_weighted = c * ac_k
print(f'  c sum: original={c.sum():.0f}, weighted={c_weighted.sum():.0f}  (ratio={c_weighted.sum()/c.sum():.2f})', flush=True)
t0=time.time()
hH, infoH = em_sparse(c_weighted, cn, max_iter=200)
print(f'  done ({time.time()-t0:.0f}s, iters={infoH["iterations"]}, delta={infoH["final_delta"]:.2e})', flush=True)
results.append(report('H_carrier_weighted', hH, h_truth, is_cactus, iters=infoH['iterations']))

# H2) sqrt(ac_k) weighting — milder version
print('\nH2) sqrt(ac_k) weighted EM (c_new = c * sqrt(ac_k)) ...', flush=True)
c_sqrt = c * np.sqrt(ac_k)
t0=time.time()
hH2, infoH2 = em_sparse(c_sqrt, cn, max_iter=200)
print(f'  done ({time.time()-t0:.0f}s, iters={infoH2["iterations"]}, delta={infoH2["final_delta"]:.2e})', flush=True)
results.append(report('H2_sqrt_ac_weighted', hH2, h_truth, is_cactus, iters=infoH2['iterations']))

# H3) cap update: replace ac_k weighting with min(ac_k, 50)
print('\nH3) capped ac_k=min(ac, 50) weighted EM ...', flush=True)
c_capped = c * np.minimum(ac_k, 50)
t0=time.time()
hH3, infoH3 = em_sparse(c_capped, cn, max_iter=200)
print(f'  done ({time.time()-t0:.0f}s, iters={infoH3["iterations"]}, delta={infoH3["final_delta"]:.2e})', flush=True)
results.append(report('H3_cap_ac_50', hH3, h_truth, is_cactus, iters=infoH3['iterations']))

# Saved + uniform refs
results.append(report('saved_h_v3', h_v3_saved, h_truth, is_cactus))
results.append(report('uniform_1_231', np.full(231, 1/231), h_truth, is_cactus))

print('\n=== Summary: distance from recipe truth ===', flush=True)
df = pd.DataFrame(results)[['name','iters','L1','L2','r','cactus_mean','PG_mean']]
df['cactus_bias_pct'] = (df['cactus_mean']/h_truth[is_cactus].mean() - 1) * 100
df['PG_bias_pct']     = (df['PG_mean']/h_truth[~is_cactus].mean() - 1) * 100
print(df.to_string(index=False, float_format='%.4f'), flush=True)

print(f'\nRecipe-truth means:  cactus={h_truth[is_cactus].mean():.5f}  PG={h_truth[~is_cactus].mean():.5f}', flush=True)

# Save
np.savez(ROOT/'jf_chr1/carrier_weighted_results.npz',
         hA=hA, hH=hH, hH2=hH2, hH3=hH3,
         h_truth=h_truth, h_v3_saved=h_v3_saved,
         founders=np.array(founders), is_cactus=is_cactus, ac_k=ac_k)
print(f'\nsaved {ROOT}/jf_chr1/carrier_weighted_results.npz', flush=True)
