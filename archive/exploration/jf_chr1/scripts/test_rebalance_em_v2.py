"""Test column/per-bubble normalizations of cn_full for cactus_em on SEEDMIX_S1.

Strategies tested:
  A) baseline: standard EM (matches production)
  B) row-normalize: cn[f,k]/K_f                         (BALANCING_KMERS opt 1; already shown bad)
  C) column-normalize by carrier count: cn[:,k]/Σf cn[f,k]
  D) per-bubble: cn[:,k] / n_kmers_in_bubble(k)          (freqk-style)
  E) drop singletons: drop k-mers with carrier_count == 1
  F) drop very-common k-mers: drop k-mers with carrier_count >= 50
  G) drop both ends: keep k-mers with 2 <= carrier_count < 50
"""
import numpy as np
import pandas as pd
import scipy.sparse as sp
from pathlib import Path
import time

ROOT = Path(__file__).resolve().parents[2]


def em_sparse(counts, cn_sparse, max_iter=200, tol=1e-7, h_init=None):
    """Sparse EM. cn_sparse: (F,K) csr float32. counts: (K,) float32."""
    F, K = cn_sparse.shape
    counts = counts.astype(np.float32)
    h = np.full(F, 1.0/F, dtype=np.float32) if h_init is None else h_init.astype(np.float32)
    cn_T = cn_sparse.T.tocsr()
    cn = cn_sparse.tocsr()
    history = []
    for it in range(max_iter):
        denom = cn_T @ h
        np.maximum(denom, np.float32(1e-7), out=denom)
        cw = counts / denom
        em_term = h * (cn @ cw)
        total_c = counts.sum()
        h_new = em_term / max(float(total_c), 1e-12)
        h_new = h_new / h_new.sum()
        delta = float(np.linalg.norm(h_new - h))
        history.append(delta)
        h = h_new
        if delta < tol:
            break
    return h.astype(np.float64), {'iterations': it+1, 'converged': delta < tol, 'final_delta': delta}


def report(name, h, h_truth, is_cactus):
    l1 = float(np.abs(h - h_truth).sum())
    l2 = float(((h - h_truth)**2).sum()**0.5)
    r = float(np.corrcoef(h, h_truth)[0,1]) if h.std()>0 else float('nan')
    return dict(
        name=name, L1=l1, L2=l2, r=r,
        cactus_mean=float(h[is_cactus].mean()),
        PG_mean=float(h[~is_cactus].mean()),
    )


# Load
print('Loading...', flush=True)
t0 = time.time()
cn_full = sp.load_npz(ROOT/'poolfreq/data/cn_full_231_v3/cn_Chr1.cn.npz').tocsr().astype(np.float32)
m1 = np.load(ROOT/'poolfreq/data/cn_full_231_v3/cn_Chr1.meta.npz', allow_pickle=True)
founders = list(np.asarray(m1['founders']).astype(str))
bubble_id = np.asarray(m1['bubble_id']).astype(np.int64)
counts_all = np.load(ROOT/'poolfreq/data/seedmix_S1_counts_genomewide.npz')['counts']
n_chr1 = cn_full.shape[1]
counts_chr1 = counts_all[:n_chr1].astype(np.float32)
print(f'cn={cn_full.shape}, nnz={cn_full.nnz:,}; loaded in {time.time()-t0:.0f}s', flush=True)

# Filter to nonzero counts
nz = counts_chr1 > 0
c = counts_chr1[nz]
cn = cn_full[:, nz].tocsr()
bid_nz = bubble_id[nz]
print(f'filtered: cn={cn.shape}, nnz={cn.nnz:,}', flush=True)

# Helpers
K_f = np.asarray(cn_full.sum(axis=1)).ravel().astype(np.float32)
ac_k = np.asarray(cn.sum(axis=0)).ravel().astype(np.float32)  # per-kmer carrier count
print(f'ac_k: min={ac_k.min():.0f} median={np.median(ac_k):.0f} max={ac_k.max():.0f} '
      f'(% singletons={(ac_k==1).mean()*100:.1f}, %>=50={(ac_k>=50).mean()*100:.1f})', flush=True)

cactus_ids = set((ROOT/'pangenie_genotyping/data/merged/cactus_overlap_80.txt').read_text().split())
is_cactus = np.array([f in cactus_ids for f in founders])

recipe = pd.read_csv(ROOT/'data/seedmix_recipe_normalized.tsv', sep='\t', dtype={'ID':str})
h_truth = np.array([dict(zip(recipe['ID'], recipe['seed_prop'])).get(f, 1/231) for f in founders])

h_v3_saved = np.load(ROOT/'poolfreq/results/seedmix_231_v3_perchrom/SEEDMIX_S1.h_per_chrom.npz', allow_pickle=True)['Chr1']

# Per-bubble: n k-mers per bubble (on the filtered set)
ubids, n_per_bub = np.unique(bid_nz, return_counts=True)
nkmers_per_bub = np.zeros(int(ubids.max())+1, dtype=np.float32)
nkmers_per_bub[ubids] = n_per_bub
n_per_kmer = nkmers_per_bub[bid_nz]
print(f'bubbles: {len(ubids):,}  median n_k_per_bub={np.median(n_per_bub):.0f}  '
      f'mean={n_per_bub.mean():.1f}  max={n_per_bub.max():.0f}', flush=True)

# === Run all scenarios ===
results = []

# Baseline
print('\nA) baseline EM ...', flush=True)
t0=time.time()
hA, info = em_sparse(c, cn, max_iter=200)
print(f'  iters={info["iterations"]} delta={info["final_delta"]:.2e}  ({time.time()-t0:.0f}s)', flush=True)
results.append({**report('A_baseline', hA, h_truth, is_cactus), 'iters': info['iterations']})

# Row normalize (per-founder)
print('\nB) row-normalize ...', flush=True)
cn_row = (sp.diags(1.0/K_f) @ cn).astype(np.float32).tocsr()
t0=time.time()
hB, info = em_sparse(c, cn_row, max_iter=200)
print(f'  iters={info["iterations"]} delta={info["final_delta"]:.2e}  ({time.time()-t0:.0f}s)', flush=True)
results.append({**report('B_row_normalize', hB, h_truth, is_cactus), 'iters': info['iterations']})

# Column normalize (per-kmer by carrier count)
print('\nC) column-normalize by carrier count ...', flush=True)
inv_ac = 1.0 / np.maximum(ac_k, 1)
cn_col = (cn.multiply(inv_ac)).astype(np.float32).tocsr()
t0=time.time()
hC, info = em_sparse(c, cn_col, max_iter=200)
print(f'  iters={info["iterations"]} delta={info["final_delta"]:.2e}  ({time.time()-t0:.0f}s)', flush=True)
results.append({**report('C_col_normalize', hC, h_truth, is_cactus), 'iters': info['iterations']})

# Per-bubble (freqk style)
print('\nD) per-bubble normalize ...', flush=True)
inv_npk = 1.0 / np.maximum(n_per_kmer, 1)
cn_bub = (cn.multiply(inv_npk)).astype(np.float32).tocsr()
t0=time.time()
hD, info = em_sparse(c, cn_bub, max_iter=200)
print(f'  iters={info["iterations"]} delta={info["final_delta"]:.2e}  ({time.time()-t0:.0f}s)', flush=True)
results.append({**report('D_per_bubble', hD, h_truth, is_cactus), 'iters': info['iterations']})

# Drop singletons
print('\nE) drop singletons (ac_k==1) ...', flush=True)
keep_E = ac_k > 1
print(f'  kept {int(keep_E.sum()):,} of {len(ac_k):,} k-mers', flush=True)
cn_E = cn[:, keep_E].tocsr()
c_E = c[keep_E]
t0=time.time()
hE, info = em_sparse(c_E, cn_E, max_iter=200)
print(f'  iters={info["iterations"]} delta={info["final_delta"]:.2e}  ({time.time()-t0:.0f}s)', flush=True)
results.append({**report('E_drop_singletons', hE, h_truth, is_cactus), 'iters': info['iterations']})

# Drop very common (ac_k >= 50)
print('\nF) drop common (ac_k >= 50) ...', flush=True)
keep_F = ac_k < 50
print(f'  kept {int(keep_F.sum()):,} of {len(ac_k):,} k-mers', flush=True)
cn_F = cn[:, keep_F].tocsr()
c_F = c[keep_F]
t0=time.time()
hF, info = em_sparse(c_F, cn_F, max_iter=200)
print(f'  iters={info["iterations"]} delta={info["final_delta"]:.2e}  ({time.time()-t0:.0f}s)', flush=True)
results.append({**report('F_drop_common', hF, h_truth, is_cactus), 'iters': info['iterations']})

# Keep middle (2 <= ac_k < 50)
print('\nG) keep middle band (2 <= ac_k < 50) ...', flush=True)
keep_G = (ac_k >= 2) & (ac_k < 50)
print(f'  kept {int(keep_G.sum()):,} of {len(ac_k):,} k-mers', flush=True)
cn_G = cn[:, keep_G].tocsr()
c_G = c[keep_G]
t0=time.time()
hG, info = em_sparse(c_G, cn_G, max_iter=200)
print(f'  iters={info["iterations"]} delta={info["final_delta"]:.2e}  ({time.time()-t0:.0f}s)', flush=True)
results.append({**report('G_keep_2_to_49', hG, h_truth, is_cactus), 'iters': info['iterations']})

# Saved + uniform refs
results.append({**report('saved_h_v3', h_v3_saved, h_truth, is_cactus), 'iters': None})
results.append({**report('uniform_1_231', np.full(231, 1/231), h_truth, is_cactus), 'iters': None})

# Print
print('\n=== Summary: distance from recipe truth (smaller = closer) ===', flush=True)
df = pd.DataFrame(results)[['name','iters','L1','L2','r','cactus_mean','PG_mean']]
df['cactus_vs_truth'] = df['cactus_mean'] / h_truth[is_cactus].mean() - 1
df['PG_vs_truth']     = df['PG_mean']     / h_truth[~is_cactus].mean() - 1
print(df.to_string(index=False, float_format='%.4f'), flush=True)

print(f'\nRecipe-truth means:  cactus={h_truth[is_cactus].mean():.5f}  PG={h_truth[~is_cactus].mean():.5f}', flush=True)

# Save
np.savez(ROOT/'jf_chr1/rebalance_v2_results.npz',
         hA=hA, hB=hB, hC=hC, hD=hD, hE=hE, hF=hF, hG=hG, h_truth=h_truth,
         h_v3_saved=h_v3_saved, K_f=K_f, founders=np.array(founders), is_cactus=is_cactus)
print(f'saved {ROOT}/jf_chr1/rebalance_v2_results.npz', flush=True)
