"""Run baseline + 3 imbalance-mitigation EM variants on the same reads.

Approaches:
  1. Stratified EM: scale counts at PRIVATE k-mers (ac_k=1) by alpha.
     alpha=0 → drop private (equivalent to filt2). alpha=1 → no change.
     Sweep alpha ∈ {0.0, 0.1, 0.25, 0.5, 0.75, 1.0}.
  2. Per-bubble per-founder cap: post-process cn_full so no founder carries
     more than k_max k-mers in any single bubble. Sweep k_max ∈ {2, 4, 8}.
  3. Per-bubble rate normalization: scale counts by 1/n_kmers_in_bubble
     so every bubble contributes ~1 unit of evidence per carrier.

Reads cov10_n200_g1 simulated reads. Runs global EM (Chr1 only) for each
variant. Saves h vectors + per-founder bias stats. Truth is pool_weights.tsv.
"""
from __future__ import annotations
import json
import sys
import time
from pathlib import Path

import numpy as np
from scipy.sparse import load_npz, csr_matrix, save_npz

ROOT = Path('/carnegie/nobackup/scratch/tbellagio/kmate')
sys.path.insert(0, str(ROOT / 'poolfreq/src'))
from em_solver import solve_em
from kmer_count import count_kmers_in_fasta

SIM = ROOT / 'sims/visor_freqk/pool_sweep_82_recomb' / 'cov10_n200_g1_s42_hotspots_p231_chr1'
CN_PREFIX = ROOT / 'poolfreq/data/cn_full_231_v3/cn'
SAMPLES_FILE = ROOT / 'data/vcf_samples_231.txt'
SPLIT_FILE = ROOT / 'data/founder_split_cactus_pg.json'
SCRATCH = ROOT / 'scratch'
SCRATCH.mkdir(exist_ok=True)


def log(*args, **kw):
    print(*args, **kw, flush=True)


def load_data():
    """Load cn_full_v3 Chr1, count reads, return (cn, counts, meta, K_f, ac_k, bubble_id)."""
    log('Loading cn_full_v3 Chr1...')
    cn = load_npz(str(CN_PREFIX) + '_Chr1.cn.npz')
    meta = np.load(str(CN_PREFIX) + '_Chr1.meta.npz', allow_pickle=True)
    F, K = cn.shape
    log(f'  cn: {F} × {K:,}, nnz={cn.nnz:,}')
    kmer_index = list(meta['kmer_index'])
    bubble_id = np.asarray(meta['bubble_id'], dtype=np.int64)
    log(f'  n_bubbles: {bubble_id.max()+1:,}')

    counts_cache = SCRATCH / 'counts_chr1_cov10_n200_g1.npz'
    if counts_cache.exists():
        log(f'Loading cached counts from {counts_cache}...')
        counts = np.load(counts_cache)['counts']
    else:
        r1 = SIM / 'reads/r1.fq'
        r2 = SIM / 'reads/r2.fq'
        log(f'Counting k-mers from {r1.name} + {r2.name} (~10 min)...')
        t = time.time()
        cd = count_kmers_in_fasta([str(r1), str(r2)], kmer_index, k=31,
                                   threads=8, hash_size='3G')
        log(f'  jellyfish done in {time.time()-t:.0f}s')
        counts = np.array([cd[km] for km in kmer_index], dtype=np.int64)
        np.savez(counts_cache, counts=counts)
        log(f'  cached at {counts_cache}')
    log(f'  counts: nonzero {(counts>0).sum():,}/{K:,}, sum={counts.sum():,}')

    log('Densifying cn for fast EM...')
    cn_dense = np.asarray(cn.todense()).astype(np.float32)
    log(f'  cn_dense: {cn_dense.nbytes/1e9:.1f} GB')

    K_f = cn_dense.sum(axis=1)
    ac_k = cn_dense.sum(axis=0).astype(np.int32)
    return cn_dense, counts, meta, K_f, ac_k, bubble_id, F, K


def load_truth():
    """Truth h vector: founder_id → expected pool fraction."""
    samples = [s.strip() for s in open(SAMPLES_FILE)]
    pw = {}
    with open(SIM / 'pool_weights.tsv') as f:
        header = next(f)
        for line in f:
            founder, count, weight = line.strip().split('\t')
            pw[founder] = float(weight)
    h_true = np.array([pw.get(s, 0.0) for s in samples], dtype=np.float64)
    h_true = h_true / h_true.sum()  # normalize defensively
    log(f'  truth: {(h_true > 0).sum()} founders with nonzero weight, '
        f'sum={h_true.sum():.4f}, max={h_true.max():.4f}, min_nz={h_true[h_true>0].min():.4f}')
    return h_true, samples


def evaluate_h(h_est, h_true, samples, is_cactus, is_pg, label):
    """Per-founder bias + summary stats by class."""
    h_est = np.asarray(h_est, dtype=np.float64)
    h_true = np.asarray(h_true, dtype=np.float64)
    bias = h_est - h_true  # raw bias
    # Per-founder over-credit ratio (h_est / h_true) for nonzero truth founders
    nz = h_true > 1e-9
    ratio = np.full_like(h_est, np.nan)
    ratio[nz] = h_est[nz] / h_true[nz]

    mae = np.mean(np.abs(bias))
    rmse = np.sqrt(np.mean(bias ** 2))
    sum_cac = h_est[is_cactus].sum()
    sum_pg = h_est[is_pg].sum()
    true_sum_cac = h_true[is_cactus].sum()
    true_sum_pg = h_true[is_pg].sum()

    # h_est mass on cactus vs PG, divided by truth: a clean "class-level over-credit"
    over_cactus = sum_cac / max(true_sum_cac, 1e-9)
    over_pg = sum_pg / max(true_sum_pg, 1e-9)

    eff_n = 1.0 / np.sum(h_est ** 2)
    log(f'\n=== {label} ===')
    log(f'  per-founder MAE = {mae:.6f}, RMSE = {rmse:.6f}, eff_n_founders = {eff_n:.1f}')
    log(f'  Σh on cactus 80: est {sum_cac:.4f} vs truth {true_sum_cac:.4f} → '
        f'{over_cactus:.3f}× ({100*(over_cactus-1):+.1f}%)')
    log(f'  Σh on PG    151: est {sum_pg:.4f} vs truth {true_sum_pg:.4f} → '
        f'{over_pg:.3f}× ({100*(over_pg-1):+.1f}%)')
    return dict(label=label, mae=mae, rmse=rmse, eff_n=eff_n,
                over_cactus=float(over_cactus), over_pg=float(over_pg),
                h=h_est.tolist())


def run_baseline(cn, counts):
    """Standard EM, no modifications."""
    log('\n>>> BASELINE: standard EM')
    nz = counts > 0
    cn_em = np.ascontiguousarray(cn[:, nz])
    c_em = counts[nz].astype(np.float32)
    t = time.time()
    h, info = solve_em(c_em, cn_em, coverage=0.0, max_iter=200, tol=1e-7)
    log(f'  converged in {info["iterations"]} iters [{time.time()-t:.0f}s]')
    return h


def run_stratified(cn, counts, ac_k, alpha):
    """Approach 1: scale counts at private k-mers by alpha. alpha=0 → drop private."""
    log(f'\n>>> APPROACH 1: stratified EM, alpha={alpha}')
    w = np.where(ac_k == 1, alpha, 1.0).astype(np.float32)
    counts_w = counts.astype(np.float32) * w
    nz = counts_w > 0
    cn_em = np.ascontiguousarray(cn[:, nz])
    c_em = counts_w[nz]
    n_priv = int((ac_k[nz] == 1).sum())
    log(f'  applied weight alpha={alpha} on {(ac_k==1).sum():,} private k-mers; '
        f'{n_priv:,} of those have nonzero counts')
    t = time.time()
    h, info = solve_em(c_em, cn_em, coverage=0.0, max_iter=200, tol=1e-7)
    log(f'  converged in {info["iterations"]} iters [{time.time()-t:.0f}s]')
    return h


def run_perbubble_cap(cn_dense, counts, bubble_id, k_max, rng_seed=42):
    """Approach 2: cap each founder to <= k_max k-mers per bubble. Random drop."""
    log(f'\n>>> APPROACH 2: per-bubble per-founder cap, k_max={k_max}')
    F, K = cn_dense.shape
    # Precompute per-bubble column lists
    n_bubbles = int(bubble_id.max() + 1)
    log(f'  Building per-bubble column lookup (n_bubbles={n_bubbles:,})...')
    order = np.argsort(bubble_id, kind='stable')
    # Build offsets via bincount
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
        sub = cn_capped[:, cols]  # view
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

    # Now run EM with capped cn
    nz = counts > 0
    cn_em = np.ascontiguousarray(cn_capped[:, nz])
    c_em = counts[nz].astype(np.float32)
    t = time.time()
    h, info = solve_em(c_em, cn_em, coverage=0.0, max_iter=200, tol=1e-7)
    log(f'  EM converged in {info["iterations"]} iters [{time.time()-t:.0f}s]')
    return h


def run_perbubble_rate(cn, counts, bubble_id):
    """Approach 3: per-bubble rate normalization. counts_k → counts_k / n_kmers_in_b."""
    log(f'\n>>> APPROACH 3: per-bubble rate norm (counts /= n_kmers_in_bubble)')
    n_per_b = np.bincount(bubble_id).astype(np.float32)
    w = 1.0 / n_per_b[bubble_id]
    counts_w = counts.astype(np.float32) * w
    log(f'  median n_kmers/bubble = {int(np.median(n_per_b))}; '
        f'effective per-bubble evidence weight = ~1 unit')
    nz = counts_w > 1e-9
    cn_em = np.ascontiguousarray(cn[:, nz])
    c_em = counts_w[nz]
    t = time.time()
    h, info = solve_em(c_em, cn_em, coverage=0.0, max_iter=200, tol=1e-7)
    log(f'  converged in {info["iterations"]} iters [{time.time()-t:.0f}s]')
    return h


def main():
    cn, counts, meta, K_f, ac_k, bubble_id, F, K = load_data()
    h_true, samples = load_truth()

    # Class membership
    split = json.load(open(SPLIT_FILE))
    cactus_set = set(split['cactus'])
    pg_set = set(split['PG'])
    is_cactus = np.array([s in cactus_set for s in samples], dtype=bool)
    is_pg = np.array([s in pg_set for s in samples], dtype=bool)

    results = []

    # Baseline
    h_base = run_baseline(cn, counts)
    results.append(evaluate_h(h_base, h_true, samples, is_cactus, is_pg, 'baseline'))

    # Approach 1: stratified
    for alpha in [0.0, 0.1, 0.25, 0.5, 0.75]:
        h_a = run_stratified(cn, counts, ac_k, alpha=alpha)
        results.append(evaluate_h(h_a, h_true, samples, is_cactus, is_pg,
                                   f'stratified_alpha={alpha}'))

    # Approach 2: per-bubble per-founder cap
    for k_max in [2, 4, 8]:
        h_c = run_perbubble_cap(cn, counts, bubble_id, k_max=k_max)
        results.append(evaluate_h(h_c, h_true, samples, is_cactus, is_pg,
                                   f'perbubble_cap_kmax={k_max}'))

    # Approach 3: per-bubble rate normalization
    h_r = run_perbubble_rate(cn, counts, bubble_id)
    results.append(evaluate_h(h_r, h_true, samples, is_cactus, is_pg,
                               'perbubble_rate_norm'))

    # Save all results
    out_json = SCRATCH / 'three_approaches_results.json'
    with open(out_json, 'w') as f:
        json.dump({'results': results, 'samples': samples,
                   'h_true': h_true.tolist(),
                   'is_cactus': is_cactus.tolist(),
                   'is_pg': is_pg.tolist()}, f)
    log(f'\nSaved {out_json}')

    # Summary table
    log('\n=== SUMMARY ===')
    log(f'{"label":40s}  {"MAE":>10s}  {"RMSE":>10s}  {"eff_n":>8s}  '
        f'{"over_cac":>10s}  {"over_pg":>10s}')
    for r in results:
        log(f'{r["label"]:40s}  {r["mae"]:>10.6f}  {r["rmse"]:>10.6f}  '
            f'{r["eff_n"]:>8.1f}  {r["over_cactus"]:>9.3f}×  {r["over_pg"]:>9.3f}×')


if __name__ == '__main__':
    main()
