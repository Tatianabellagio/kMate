"""Plot per-founder h-bias by class for all three approaches + AF-level eval.

Inputs: scratch/three_approaches_results.json (from three_approaches.py),
truth_v3 at recomb_truth_v3.tsv.gz in the sim dir.

Outputs:
  scratch/three_approaches_h_bias.png  — per-founder h-bias, faceted by approach
  scratch/three_approaches_af_eval.tsv — per-approach per-record AF accuracy
"""
from __future__ import annotations
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.sparse import load_npz

ROOT = Path('/carnegie/nobackup/scratch/tbellagio/hapfire_sv')
SIM = ROOT / 'sims/visor_freqk/pool_sweep_82_recomb' / 'cov10_n200_g1_s42_hotspots_p231_chr1'
CN_VAR = ROOT / 'poolfreq/data/cn_var_231_v3.cn_var.npz'
CN_VAR_META = ROOT / 'poolfreq/data/cn_var_231_v3.meta.npz'
SCRATCH = ROOT / 'scratch'


def log(*a, **kw): print(*a, **kw, flush=True)


def main():
    with open(SCRATCH / 'three_approaches_results.json') as f:
        data = json.load(f)
    results = data['results']
    h_true = np.array(data['h_true'])
    is_cactus = np.array(data['is_cactus'])
    is_pg = np.array(data['is_pg'])
    samples = data['samples']

    # ---- Plot 1: per-founder h_est vs h_true, faceted by approach ----
    n = len(results)
    ncols = 4
    nrows = (n + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(4.5 * ncols, 4.2 * nrows),
                              sharex=True, sharey=True)
    axes = axes.flatten()
    xmax = max(h_true.max(), 0.012)
    for ax, r in zip(axes, results):
        h_est = np.array(r['h'])
        ax.scatter(h_true[is_cactus], h_est[is_cactus],
                   s=15, alpha=0.55, color='#d62728', label='cactus 80', edgecolors='none')
        ax.scatter(h_true[is_pg], h_est[is_pg],
                   s=15, alpha=0.55, color='#1f77b4', label='PG 151', edgecolors='none')
        ax.plot([0, xmax], [0, xmax], 'k--', lw=0.7, alpha=0.6)
        ax.set_xlim(0, xmax); ax.set_ylim(0, xmax * 1.5)
        ax.set_title(f'{r["label"]}\nMAE={r["mae"]:.5f}  cac={r["over_cactus"]:.2f}×  '
                     f'pg={r["over_pg"]:.2f}×', fontsize=9)
        if ax is axes[0]:
            ax.legend(fontsize=8)
    for ax in axes[n:]:
        ax.axis('off')
    fig.supxlabel('h_true (per founder)')
    fig.supylabel('h_est')
    fig.suptitle('Per-founder h_est vs h_true on cov10_n200_g1 (Chr1)', y=1.001)
    fig.tight_layout()
    out_png = SCRATCH / 'three_approaches_h_bias.png'
    fig.savefig(out_png, dpi=120, bbox_inches='tight')
    log(f'Wrote {out_png}')

    # ---- Plot 2: class over-credit bar chart ----
    fig, ax = plt.subplots(figsize=(11, 5))
    labels = [r['label'] for r in results]
    over_cac = [r['over_cactus'] for r in results]
    over_pg = [r['over_pg'] for r in results]
    x = np.arange(len(labels))
    w = 0.4
    ax.bar(x - w/2, over_cac, w, color='#d62728', label='cactus 80 mass / truth', edgecolor='black')
    ax.bar(x + w/2, over_pg, w, color='#1f77b4', label='PG 151 mass / truth', edgecolor='black')
    ax.axhline(1.0, color='k', linestyle='--', lw=0.8)
    ax.set_ylabel('Σh_est on class / Σh_true on class')
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=30, ha='right', fontsize=9)
    ax.legend(fontsize=9, loc='upper right')
    ax.set_title('Class-level mass over-credit (1.0 = unbiased)')
    fig.tight_layout()
    out_png2 = SCRATCH / 'three_approaches_class_over_credit.png'
    fig.savefig(out_png2, dpi=120, bbox_inches='tight')
    log(f'Wrote {out_png2}')

    # ---- Per-record AF accuracy ----
    log('\nLoading truth_v3 + cn_var_v3 for AF eval...')
    truth_df = pd.read_csv(SIM / 'recomb_truth_v3.tsv.gz', sep='\t', compression='gzip')
    truth = truth_df.truth_af.to_numpy().astype(np.float32)
    log(f'  truth records: {len(truth):,}')

    cn_var = load_npz(CN_VAR)
    cn_var_meta = np.load(CN_VAR_META, allow_pickle=True)
    rec_chrom = np.asarray(cn_var_meta['chrom']).astype(str)
    chr1_idx = np.where(rec_chrom == 'Chr1')[0]
    log(f'  cn_var_v3 Chr1: {len(chr1_idx):,} records')
    cn_var_chr1 = cn_var[:, chr1_idx]

    # Filter to polymorphic + nonzero subset on Chr1 only
    if len(chr1_idx) != len(truth):
        log(f'  WARN: cn_var Chr1 ({len(chr1_idx)}) != truth ({len(truth)}); slicing truth to first n')
        n_min = min(len(chr1_idx), len(truth))
        truth = truth[:n_min]
        cn_var_chr1 = cn_var_chr1[:, :n_min]

    poly_mask = (truth > 0) & (truth < 1)
    log(f'  polymorphic records: {poly_mask.sum():,}')
    truth_poly = truth[poly_mask]

    rows = []
    for r in results:
        h_est = np.array(r['h'], dtype=np.float64)
        af_pred = (h_est @ cn_var_chr1.toarray()).astype(np.float32)
        af_pred_poly = af_pred[poly_mask]
        err = af_pred_poly - truth_poly
        mae = float(np.mean(np.abs(err)))
        rmse = float(np.sqrt(np.mean(err ** 2)))
        r_pearson = float(np.corrcoef(truth_poly, af_pred_poly)[0, 1])
        slope = float(np.polyfit(truth_poly, af_pred_poly, 1)[0])
        out_frac = float(np.mean(np.abs(err) > 0.1))
        rows.append(dict(label=r['label'], mae=mae, rmse=rmse, r2=r_pearson**2,
                         slope=slope, outlier_frac=out_frac))
        log(f'  {r["label"]:40s}  MAE={mae:.4f}  RMSE={rmse:.4f}  R²={r_pearson**2:.4f}  '
            f'slope={slope:+.3f}  out={100*out_frac:.2f}%')

    df = pd.DataFrame(rows)
    out_tsv = SCRATCH / 'three_approaches_af_eval.tsv'
    df.to_csv(out_tsv, sep='\t', index=False)
    log(f'\nWrote {out_tsv}')


if __name__ == '__main__':
    main()
