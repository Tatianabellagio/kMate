"""Test row-normalized EM on SEEDMIX_S1 Chr1 vs standard EM."""
import numpy as np
import pandas as pd
import scipy.sparse as sp
from pathlib import Path
import time, sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'poolfreq/src'))


def em_sparse(counts, cn_sparse, max_iter=80, tol=1e-7, h_init=None):
    """EM that keeps cn sparse to save memory.
    cn_sparse: (F, K) sparse, dtype float32
    counts: (K,) float32
    Returns h (F,) and info dict.
    """
    F, K = cn_sparse.shape
    counts = counts.astype(np.float32)
    h = np.full(F, 1.0/F, dtype=np.float32) if h_init is None else h_init.astype(np.float32)
    history = []
    cn_T = cn_sparse.T.tocsr()
    cn = cn_sparse.tocsr()
    for it in range(max_iter):
        # denom_k = h @ cn[:, k]   (K-vector)
        # equivalent to cn.T @ h
        denom = cn_T @ h
        denom = np.maximum(denom, np.float32(1e-7))
        cw = counts / denom
        # em_term[f] = h[f] * sum_k cn[f,k] * cw[k]
        em_term = h * (cn @ cw)
        total_c = counts.sum()
        h_new = em_term / max(float(total_c), 1e-12)
        h_new = h_new / h_new.sum()
        delta = float(np.linalg.norm(h_new - h))
        history.append(delta)
        h = h_new
        if delta < tol:
            break
    return h.astype(np.float64), {'iterations': it+1, 'history': history, 'converged': delta < tol}


# Load
print('Loading cn_full Chr1...', flush=True)
t0=time.time()
cn_full = sp.load_npz(ROOT/'poolfreq/data/cn_full_231_v3/cn_Chr1.cn.npz').tocsr()
print(f'  cn_full: {cn_full.shape}, nnz={cn_full.nnz:,}  ({time.time()-t0:.0f}s)', flush=True)

m1 = np.load(ROOT/'poolfreq/data/cn_full_231_v3/cn_Chr1.meta.npz', allow_pickle=True)
founders = list(np.asarray(m1['founders']).astype(str))
n_chr1 = len(m1['kmer_index'])

counts_all = np.load(ROOT/'poolfreq/data/seedmix_S1_counts_genomewide.npz')['counts']
counts_chr1 = counts_all[:n_chr1].astype(np.float32)
print(f'  Chr1 counts: nz={int((counts_chr1>0).sum()):,}', flush=True)

# Filter to nonzero
nz_mask = counts_chr1 > 0
c_nz = counts_chr1[nz_mask]
cn_nz = cn_full[:, nz_mask].astype(np.float32).tocsr()
print(f'  filtered: cn={cn_nz.shape}, nnz={cn_nz.nnz:,}', flush=True)

K_f = np.asarray(cn_full.sum(axis=1)).ravel().astype(np.float32)
cactus_ids = set((ROOT/'pangenie_genotyping/data/merged/cactus_overlap_80.txt').read_text().split())
is_cactus = np.array([f in cactus_ids for f in founders])
print(f'  K_f cactus mean={K_f[is_cactus].mean():,.0f}  PG mean={K_f[~is_cactus].mean():,.0f}', flush=True)

# STANDARD EM
print('\n=== STANDARD EM ===', flush=True)
t0=time.time()
h_std, info1 = em_sparse(c_nz, cn_nz, max_iter=80)
print(f'  done in {time.time()-t0:.0f}s, iters={info1["iterations"]}, converged={info1["converged"]}', flush=True)

h_npz = np.load(ROOT/'poolfreq/results/seedmix_231_v3_perchrom/SEEDMIX_S1.h_per_chrom.npz', allow_pickle=True)
h_v3 = np.asarray(h_npz['Chr1']).astype(np.float64)
print(f'  std vs saved h_v3: r={np.corrcoef(h_std, h_v3)[0,1]:.4f}  max|diff|={np.abs(h_std-h_v3).max():.5f}', flush=True)

# ROW-NORMALIZED EM
print('\n=== ROW-NORMALIZED EM (cn_norm = cn / K_f[:, None]) ===', flush=True)
# Use sp.diags to scale rows efficiently — sparse-friendly
D_inv = sp.diags(1.0 / K_f.astype(np.float64))
cn_norm = (D_inv @ cn_nz).astype(np.float32).tocsr()
print(f'  cn_norm: nnz={cn_norm.nnz:,}', flush=True)

t0=time.time()
h_norm, info2 = em_sparse(c_nz, cn_norm, max_iter=80)
print(f'  done in {time.time()-t0:.0f}s, iters={info2["iterations"]}, converged={info2["converged"]}', flush=True)

# Convert to literal mass: h_mass = h_norm * K_f, renormalize
h_mass = h_norm * K_f
h_mass = h_mass / h_mass.sum()

# COMPARE TO RECIPE TRUTH
recipe = pd.read_csv(ROOT/'data/seedmix_recipe_normalized.tsv', sep='\t', dtype={'ID':str})
recipe_map = dict(zip(recipe['ID'], recipe['seed_prop']))
h_truth = np.array([recipe_map.get(f, 1/231) for f in founders])

print('\n=== Comparison: distance from recipe truth ===', flush=True)
print(f'{"method":<28s}  {"L1":>8s}  {"L2":>8s}  {"r":>8s}  {"cactus mean":>12s}  {"PG mean":>10s}', flush=True)
for name, h in [('uniform 1/231',           np.full(231, 1/231)),
                ('h_v3 (saved, original)',  h_v3),
                ('h_std (recomputed EM)',   h_std),
                ('h_norm (row-normalized)', h_norm),
                ('h_mass (= h_norm * K_f)', h_mass),
                ('recipe (self check)',     h_truth)]:
    l1 = float(np.abs(h - h_truth).sum())
    l2 = float(((h - h_truth)**2).sum()**0.5)
    r = float(np.corrcoef(h, h_truth)[0,1]) if h.std()>0 else float('nan')
    print(f'  {name:<28s}  {l1:>8.4f}  {l2:>8.4f}  {r:>8.4f}  {h[is_cactus].mean():>12.5f}  {h[~is_cactus].mean():>10.5f}', flush=True)

print(f'\nRecipe-truth means: cactus={h_truth[is_cactus].mean():.5f}, PG={h_truth[~is_cactus].mean():.5f}', flush=True)

# Save results
np.savez(ROOT/'jf_chr1/rebalance_results.npz',
         h_std=h_std, h_norm=h_norm, h_mass=h_mass, h_truth=h_truth,
         h_v3_saved=h_v3, K_f=K_f, founders=np.array(founders), is_cactus=is_cactus)
print(f'\nsaved {ROOT}/jf_chr1/rebalance_results.npz', flush=True)
