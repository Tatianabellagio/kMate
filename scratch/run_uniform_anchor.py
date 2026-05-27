"""Run cactus_em with a Dirichlet anchor toward uniform 1/F (not toward h_global).

This is the under-explored lever for the over-concentration finding:
hapFIRE has eff_n=187/231 on SEEDMIX_S1, cactus_em has eff_n=128/231 —
cactus_em concentrates on too few founders. A uniform prior pulls every
founder toward 1/F.

Uses solve_em's existing prior_h + prior_weight infrastructure with
prior_h = uniform = ones(F)/F.

Each SLURM array task picks a different prior_weight from CONFIGS.
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

ROOT = Path('/carnegie/nobackup/scratch/tbellagio/hapfire_sv')
sys.path.insert(0, str(ROOT / 'poolfreq/src'))
from em_solver import solve_em

SIM = ROOT / 'sims/visor_freqk/pool_sweep_82_recomb' / 'cov10_n200_g1_s42_hotspots_p231_chr1'
CN_PREFIX = ROOT / 'poolfreq/data/cn_full_231_v3/cn'
COUNTS_CACHE = ROOT / 'scratch/counts_chr1_cov10_n200_g1.npz'
OUT_DIR = ROOT / 'scratch/three_approaches'
OUT_DIR.mkdir(parents=True, exist_ok=True)

# (label, prior_weight). prior_h = uniform 1/F.
# Calibration: total counts ~17.5M, F=231 → per-founder data weight ~76K.
# prior_weight λ multiplies (total_c * 1/F) per founder. So λ=0.3 → ~30% prior;
# λ=1 → ~equal prior+data; λ=10 → prior dominates.
CONFIGS = [
    ('uniform_anchor_lambda=0.1',  0.1),
    ('uniform_anchor_lambda=0.3',  0.3),
    ('uniform_anchor_lambda=0.5',  0.5),
    ('uniform_anchor_lambda=1.0',  1.0),
    ('uniform_anchor_lambda=2.0',  2.0),
    ('uniform_anchor_lambda=5.0',  5.0),
]


def log(*a, **kw): print(*a, **kw, flush=True)


def main():
    task = int(sys.argv[1]) - 1
    label, lam = CONFIGS[task]
    out_path = OUT_DIR / f'h_task20_uniform_lambda={lam}.npz'
    log(f'=== TASK {task} ({label}, prior_weight={lam}) ===')

    t0 = time.time()
    log('Loading cn_full_v3 Chr1...')
    cn = load_npz(str(CN_PREFIX) + '_Chr1.cn.npz')
    F, K = cn.shape
    log(f'  cn: {F} × {K:,}')
    log(f'Loading counts from {COUNTS_CACHE}...')
    counts = np.load(COUNTS_CACHE)['counts']
    log('Densifying...')
    cn_dense = np.asarray(cn.todense()).astype(np.float32)
    del cn; gc.collect()
    log(f'  cn_dense: {cn_dense.nbytes/1e9:.1f} GB, load {time.time()-t0:.0f}s')

    nz = counts > 0
    cn_em = np.ascontiguousarray(cn_dense[:, nz])
    c_em = counts[nz].astype(np.float32)
    del cn_dense; gc.collect()

    # uniform prior
    prior_h = np.full(F, 1.0 / F, dtype=np.float32)

    t = time.time()
    h, info = solve_em(c_em, cn_em, 0.0, max_iter=200, tol=1e-7,
                        prior_h=prior_h, prior_weight=lam)
    em_wall = time.time() - t
    log(f'  EM done: iters={info["iterations"]}, wall={em_wall:.0f}s, '
        f'eff_n={1/np.sum(h**2):.1f}')

    np.savez(out_path,
             h=h.astype(np.float64),
             label=label, approach='uniform_anchor',
             params=json.dumps({'prior_weight': lam}),
             em_iters=info['iterations'],
             em_wall_sec=em_wall)
    log(f'  saved {out_path}')


if __name__ == '__main__':
    main()
