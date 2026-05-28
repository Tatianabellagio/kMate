"""Run ONE mitigation approach on cov10_n200_g1 reads and save the h vector.

Designed for SLURM array: pass SLURM_ARRAY_TASK_ID via $1 to pick a config from
APPROACHES below. Each task loads cn_full_v3 Chr1, applies its approach, runs
the EM, and writes <out_dir>/h_<task>.npz with {h, label, meta}.

Aggregation is a separate step (aggregate_three_approaches.py) that reads all
h_<task>.npz files and computes per-founder + per-record AF stats.
"""
from __future__ import annotations
import gc
import json
import os
import sys
import time
from pathlib import Path

import numpy as np
from scipy.sparse import load_npz

ROOT = Path('/carnegie/nobackup/scratch/tbellagio/kmate')
sys.path.insert(0, str(ROOT / 'src'))
from em_solver import solve_em
from kmer_count import count_kmers_in_fasta

SIM = ROOT / 'sims/visor_freqk/pool_sweep_82_recomb' / 'cov10_n200_g1_s42_hotspots_p231_chr1'
CN_PREFIX = ROOT / 'data/cn_full_231_v3/cn'
SCRATCH = ROOT / 'scratch'
OUT_DIR = SCRATCH / 'three_approaches'
OUT_DIR.mkdir(exist_ok=True)
COUNTS_CACHE = SCRATCH / 'counts_chr1_cov10_n200_g1.npz'

# Approach config: (label, approach_type, params)
APPROACHES = [
    ('baseline',                 'baseline',  {}),
    ('stratified_alpha=0.0',     'stratified', {'alpha': 0.0}),
    ('stratified_alpha=0.1',     'stratified', {'alpha': 0.1}),
    ('stratified_alpha=0.25',    'stratified', {'alpha': 0.25}),
    ('stratified_alpha=0.5',     'stratified', {'alpha': 0.5}),
    ('stratified_alpha=0.75',    'stratified', {'alpha': 0.75}),
    ('stratified_alpha=1.5',     'stratified', {'alpha': 1.5}),  # upweight private
    ('stratified_alpha=2.0',     'stratified', {'alpha': 2.0}),  # upweight private
    ('perbubble_cap_kmax=2',     'cap',        {'k_max': 2}),
    ('perbubble_cap_kmax=4',     'cap',        {'k_max': 4}),
    ('perbubble_cap_kmax=8',     'cap',        {'k_max': 8}),
    ('perbubble_rate_norm',      'rate_norm',  {}),
]


def log(*a, **kw): print(*a, **kw, flush=True)


def load_cn_and_counts():
    """Load cn_full_v3 Chr1 (densify) and the cached k-mer counts."""
    t0 = time.time()
    log('Loading cn_full_v3 Chr1...')
    cn = load_npz(str(CN_PREFIX) + '_Chr1.cn.npz')
    meta = np.load(str(CN_PREFIX) + '_Chr1.meta.npz', allow_pickle=True)
    F, K = cn.shape
    log(f'  cn: {F} × {K:,}, nnz={cn.nnz:,}')
    bubble_id = np.asarray(meta['bubble_id'], dtype=np.int64)

    if COUNTS_CACHE.exists():
        log(f'Loading counts from {COUNTS_CACHE}...')
        counts = np.load(COUNTS_CACHE)['counts']
    else:
        kmer_index = list(meta['kmer_index'])
        r1 = SIM / 'reads/r1.fq'; r2 = SIM / 'reads/r2.fq'
        log(f'  Counting k-mers (~5-10 min)...')
        cd = count_kmers_in_fasta([str(r1), str(r2)], kmer_index, k=31,
                                   threads=8, hash_size='3G')
        counts = np.array([cd[km] for km in kmer_index], dtype=np.int64)
        np.savez(COUNTS_CACHE, counts=counts)
    log(f'  counts: nonzero {(counts>0).sum():,}/{K:,}')

    log('Densifying cn for fast EM...')
    cn_dense = np.asarray(cn.todense()).astype(np.float32)
    del cn; gc.collect()
    log(f'  cn_dense: {cn_dense.nbytes/1e9:.1f} GB, load {time.time()-t0:.0f}s')
    return cn_dense, counts, bubble_id


def run_baseline(cn, counts):
    nz = counts > 0
    cn_em = np.ascontiguousarray(cn[:, nz])
    c_em = counts[nz].astype(np.float32)
    h, info = solve_em(c_em, cn_em, 0.0, max_iter=200, tol=1e-7)
    return h, info


def run_stratified(cn, counts, alpha):
    """Weight private k-mer counts by alpha (relative to shared, which stay 1)."""
    ac_k = cn.sum(axis=0).astype(np.int32)
    w = np.where(ac_k == 1, alpha, 1.0).astype(np.float32)
    counts_w = counts.astype(np.float32) * w
    nz = counts_w > 1e-9
    cn_em = np.ascontiguousarray(cn[:, nz])
    c_em = counts_w[nz]
    log(f'  alpha={alpha}: weighted {(ac_k==1).sum():,} private k-mers; '
        f'cn_em shape after nz: {cn_em.shape}')
    h, info = solve_em(c_em, cn_em, 0.0, max_iter=200, tol=1e-7)
    return h, info


def run_perbubble_cap(cn_dense, counts, bubble_id, k_max, rng_seed=42):
    """Cap each founder to <= k_max k-mers per bubble. Random drop excess."""
    F, K = cn_dense.shape
    n_bubbles = int(bubble_id.max() + 1)
    # Group columns by bubble for fast iteration
    order = np.argsort(bubble_id, kind='stable')
    n_per_b = np.bincount(bubble_id, minlength=n_bubbles)
    offsets = np.concatenate([[0], np.cumsum(n_per_b)])
    cn_capped = cn_dense.copy()
    rng = np.random.default_rng(rng_seed)
    n_capped = 0
    t = time.time()
    for b in range(n_bubbles):
        cols = order[offsets[b]:offsets[b+1]]
        if len(cols) == 0:
            continue
        sub = cn_capped[:, cols]
        cnt = sub.sum(axis=1)
        over = np.where(cnt > k_max)[0]
        if over.size == 0:
            continue
        for f in over:
            carry_idx = np.where(sub[f] > 0)[0]
            drop = rng.choice(carry_idx, size=carry_idx.size - k_max, replace=False)
            cn_capped[f, cols[drop]] = 0
            n_capped += drop.size
    log(f'  capped {n_capped:,} (founder, k-mer) pairs in {time.time()-t:.0f}s')
    nz = counts > 0
    cn_em = np.ascontiguousarray(cn_capped[:, nz])
    c_em = counts[nz].astype(np.float32)
    del cn_capped; gc.collect()
    h, info = solve_em(c_em, cn_em, 0.0, max_iter=200, tol=1e-7)
    return h, info


def run_rate_norm(cn, counts, bubble_id):
    """Per-bubble rate normalization: counts /= n_kmers_in_bubble."""
    n_per_b = np.bincount(bubble_id).astype(np.float32)
    w = 1.0 / n_per_b[bubble_id]
    counts_w = counts.astype(np.float32) * w
    nz = counts_w > 1e-9
    cn_em = np.ascontiguousarray(cn[:, nz])
    c_em = counts_w[nz]
    log(f'  rate_norm: median n_kmers/bubble = {int(np.median(n_per_b))}')
    h, info = solve_em(c_em, cn_em, 0.0, max_iter=200, tol=1e-7)
    return h, info


def main():
    task = int(sys.argv[1]) - 1   # SLURM array index is 1-based
    label, approach, params = APPROACHES[task]
    out_path = OUT_DIR / f'h_task{task:02d}_{label}.npz'
    log(f'=== TASK {task} ({label}, {approach}, {params}) ===')

    cn, counts, bubble_id = load_cn_and_counts()
    F, K = cn.shape

    t = time.time()
    if approach == 'baseline':
        h, info = run_baseline(cn, counts)
    elif approach == 'stratified':
        h, info = run_stratified(cn, counts, **params)
    elif approach == 'cap':
        h, info = run_perbubble_cap(cn, counts, bubble_id, **params)
    elif approach == 'rate_norm':
        h, info = run_rate_norm(cn, counts, bubble_id)
    else:
        raise ValueError(f'unknown approach: {approach!r}')
    em_wall = time.time() - t
    log(f'  EM done: iters={info["iterations"]}, wall={em_wall:.0f}s, '
        f'eff_n={1/np.sum(h**2):.1f}')

    np.savez(out_path,
             h=h.astype(np.float64),
             label=label, approach=approach,
             params=json.dumps(params),
             em_iters=info['iterations'],
             em_wall_sec=em_wall)
    log(f'  saved {out_path}')


if __name__ == '__main__':
    main()
