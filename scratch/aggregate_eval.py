"""Aggregate 12 per-approach h files and run h-bias + AF accuracy eval.

Reads scratch/three_approaches/h_task*_*.npz, joins with cov10_n200_g1 truth,
splits by cactus/PG, computes per-record AF via cn_var_v3 projection.

Outputs:
  scratch/three_approaches_results.json        — all h vectors + class-bias stats
  scratch/three_approaches_af_eval.tsv         — per-approach AF metrics
  scratch/three_approaches_h_bias.png          — per-founder scatter, faceted by approach
  scratch/three_approaches_class_over_credit.png — class-level mass bar chart
"""
from __future__ import annotations
import glob
import json
import sys
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
H_DIR = SCRATCH / 'three_approaches'
SAMPLES_FILE = ROOT / 'data/vcf_samples_231.txt'
SPLIT_FILE = ROOT / 'data/founder_split_cactus_pg.json'


def log(*a, **kw): print(*a, **kw, flush=True)


def load_truth_and_class():
    samples = [s.strip() for s in open(SAMPLES_FILE)]
    pw = {}
    with open(SIM / 'pool_weights.tsv') as f:
        next(f)
        for line in f:
            founder, count, weight = line.strip().split('\t')
            pw[founder] = float(weight)
    h_true = np.array([pw.get(s, 0.0) for s in samples], dtype=np.float64)
    h_true = h_true / h_true.sum()
    split = json.load(open(SPLIT_FILE))
    cactus_set = set(split['cactus']); pg_set = set(split['PG'])
    is_cactus = np.array([s in cactus_set for s in samples], dtype=bool)
    is_pg = np.array([s in pg_set for s in samples], dtype=bool)
    return h_true, samples, is_cactus, is_pg


def class_stats(h_est, h_true, is_cactus, is_pg, label):
    bias = h_est - h_true
    mae = float(np.mean(np.abs(bias)))
    rmse = float(np.sqrt(np.mean(bias ** 2)))
    sum_cac = float(h_est[is_cactus].sum())
    sum_pg = float(h_est[is_pg].sum())
    true_cac = float(h_true[is_cactus].sum())
    true_pg = float(h_true[is_pg].sum())
    eff_n = float(1.0 / np.sum(h_est ** 2))
    return dict(label=label, mae=mae, rmse=rmse, eff_n=eff_n,
                sum_cactus_est=sum_cac, sum_pg_est=sum_pg,
                sum_cactus_true=true_cac, sum_pg_true=true_pg,
                over_cactus=sum_cac / max(true_cac, 1e-9),
                over_pg=sum_pg / max(true_pg, 1e-9))


def main():
    h_true, samples, is_cactus, is_pg = load_truth_and_class()
    log(f'Truth: {(h_true>0).sum()} nonzero founders, '
        f'cactus_mass={h_true[is_cactus].sum():.4f}, pg_mass={h_true[is_pg].sum():.4f}')

    # Collect all per-approach h files
    files = sorted(glob.glob(str(H_DIR / 'h_task*.npz')))
    log(f'Found {len(files)} approach result files')
    results = []
    for f in files:
        z = np.load(f, allow_pickle=True)
        label = str(z['label'])
        h = z['h']
        params = z['params'].item() if z['params'].size == 1 else str(z['params'])
        em_iters = int(z['em_iters'])
        em_wall = float(z['em_wall_sec'])
        stats = class_stats(h, h_true, is_cactus, is_pg, label)
        stats['em_iters'] = em_iters
        stats['em_wall_sec'] = em_wall
        stats['h'] = h.tolist()
        stats['params'] = params
        results.append(stats)

    # AF-level evaluation: project h through cn_var_v3 (Chr1) and compare to truth_v3
    log('\nLoading cn_var_v3 + truth_v3 for AF eval...')
    truth_df = pd.read_csv(SIM / 'recomb_truth_v3.tsv.gz', sep='\t', compression='gzip')
    truth = truth_df.truth_af.to_numpy().astype(np.float32)
    cn_var = load_npz(CN_VAR)
    cn_var_meta = np.load(CN_VAR_META, allow_pickle=True)
    rec_chrom = np.asarray(cn_var_meta['chrom']).astype(str)
    chr1_idx = np.where(rec_chrom == 'Chr1')[0]
    cn_var_chr1 = cn_var[:, chr1_idx]
    log(f'  cn_var_v3 Chr1: {cn_var_chr1.shape[1]:,} records; truth: {len(truth):,}')
    n_min = min(cn_var_chr1.shape[1], len(truth))
    truth = truth[:n_min]
    cn_var_chr1 = cn_var_chr1[:, :n_min]
    poly_mask = (truth > 0) & (truth < 1)
    truth_poly = truth[poly_mask]
    log(f'  polymorphic records: {poly_mask.sum():,}')

    # Project each approach's h through cn_var; AF metrics on polymorphic records
    af_rows = []
    cv_dense = cn_var_chr1.toarray().astype(np.float32)
    for r in results:
        h = np.array(r['h'], dtype=np.float32)
        af = (h @ cv_dense).astype(np.float32)
        af_poly = af[poly_mask]
        err = af_poly - truth_poly
        mae = float(np.mean(np.abs(err)))
        rmse = float(np.sqrt(np.mean(err ** 2)))
        r2 = float(np.corrcoef(truth_poly, af_poly)[0, 1] ** 2)
        slope = float(np.polyfit(truth_poly, af_poly, 1)[0])
        outlier_frac = float(np.mean(np.abs(err) > 0.1))
        af_rows.append(dict(label=r['label'], af_mae=mae, af_rmse=rmse,
                             af_r2=r2, af_slope=slope, af_outlier_frac=outlier_frac))
        r['af_mae'] = mae; r['af_rmse'] = rmse; r['af_r2'] = r2
        r['af_slope'] = slope; r['af_outlier_frac'] = outlier_frac

    # Include hapFIRE as a reference (row-aligned to cn_var_v3)
    hf_path = SIM / 'hapfire_v3_win10kb_perblock_proj.tsv'
    if hf_path.exists():
        hf_df = pd.read_csv(hf_path, sep='\t')
        if len(hf_df) >= n_min:
            hf_chr1 = hf_df.alt_freq.to_numpy()[chr1_idx[:n_min]] if len(hf_df) > n_min else hf_df.alt_freq.to_numpy()[:n_min]
            hf_poly = hf_chr1[poly_mask]
            finite = np.isfinite(hf_poly)
            if finite.any():
                err = hf_poly[finite] - truth_poly[finite]
                mae = float(np.mean(np.abs(err)))
                rmse = float(np.sqrt(np.mean(err ** 2)))
                r2 = float(np.corrcoef(truth_poly[finite], hf_poly[finite])[0, 1] ** 2)
                slope = float(np.polyfit(truth_poly[finite], hf_poly[finite], 1)[0])
                outlier_frac = float(np.mean(np.abs(err) > 0.1))
                af_rows.append(dict(label='[REF] hapfire_v3 10kb', af_mae=mae, af_rmse=rmse,
                                     af_r2=r2, af_slope=slope, af_outlier_frac=outlier_frac))
                log(f'  hapfire_v3 reference: MAE={mae:.5f}, R²={r2:.4f}, '
                    f'slope={slope:+.3f}, out={100*outlier_frac:.3f}%')

    # Save full results JSON
    out_json = SCRATCH / 'three_approaches_results.json'
    with open(out_json, 'w') as f:
        json.dump({'results': results, 'samples': samples,
                   'h_true': h_true.tolist(),
                   'is_cactus': is_cactus.tolist(),
                   'is_pg': is_pg.tolist()}, f)
    log(f'Wrote {out_json}')

    # Save AF eval TSV
    af_df = pd.DataFrame(af_rows)
    af_df.to_csv(SCRATCH / 'three_approaches_af_eval.tsv', sep='\t', index=False)

    # Headline table (sorted by AF MAE)
    log('\n' + '='*100)
    log(f'{"label":30s}  {"h_MAE":>9s}  {"over_cac":>8s}  {"over_pg":>8s}  '
        f'{"AF_MAE":>8s}  {"AF_R²":>7s}  {"AF_slope":>9s}  {"AF_out%":>8s}')
    log('-'*100)
    df = pd.DataFrame(results)
    df_sorted = df.sort_values('af_mae')
    for _, r in df_sorted.iterrows():
        log(f'{r["label"]:30s}  {r["mae"]:>9.6f}  {r["over_cactus"]:>7.3f}×  '
            f'{r["over_pg"]:>7.3f}×  {r["af_mae"]:>8.5f}  {r["af_r2"]:>7.4f}  '
            f'{r["af_slope"]:>+9.4f}  {100*r["af_outlier_frac"]:>7.3f}%')

    # Plot 1: per-founder scatter, faceted by approach
    n = len(results)
    ncols = 4
    nrows = (n + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(4.5 * ncols, 4.2 * nrows),
                              sharex=True, sharey=True)
    axes = axes.flatten() if hasattr(axes, 'flatten') else [axes]
    xmax = max(float(h_true.max()), 0.012)
    for ax, r in zip(axes, results):
        h_est = np.array(r['h'])
        ax.scatter(h_true[is_cactus], h_est[is_cactus],
                   s=15, alpha=0.55, color='#d62728', label='cactus 80', edgecolors='none')
        ax.scatter(h_true[is_pg], h_est[is_pg],
                   s=15, alpha=0.55, color='#1f77b4', label='PG 151', edgecolors='none')
        ax.plot([0, xmax], [0, xmax], 'k--', lw=0.7, alpha=0.6)
        ax.set_xlim(0, xmax); ax.set_ylim(0, xmax * 1.5)
        ax.set_title(f'{r["label"]}\nh_MAE={r["mae"]:.5f}  AF_MAE={r["af_mae"]:.4f}',
                     fontsize=9)
        if ax is axes[0]:
            ax.legend(fontsize=8)
    for ax in axes[len(results):]:
        ax.axis('off')
    fig.supxlabel('h_true per founder'); fig.supylabel('h_est')
    fig.suptitle('cov10_n200_g1 per-founder h_est vs h_true', y=1.001)
    fig.tight_layout()
    fig.savefig(SCRATCH / 'three_approaches_h_bias.png', dpi=120, bbox_inches='tight')
    log(f'Wrote {SCRATCH / "three_approaches_h_bias.png"}')

    # Plot 2: class over-credit bar chart
    fig, ax = plt.subplots(figsize=(13, 5))
    labels = [r['label'] for r in results]
    x = np.arange(len(labels)); w = 0.4
    ax.bar(x - w/2, [r['over_cactus'] for r in results], w,
           color='#d62728', label='Σh on cactus / truth', edgecolor='black')
    ax.bar(x + w/2, [r['over_pg'] for r in results], w,
           color='#1f77b4', label='Σh on PG / truth', edgecolor='black')
    ax.axhline(1.0, color='k', linestyle='--', lw=0.8)
    ax.set_xticks(x); ax.set_xticklabels(labels, rotation=30, ha='right', fontsize=9)
    ax.set_ylabel('Class mass ratio (est/truth); 1.0 = unbiased')
    ax.legend(); ax.set_title('Class-level over/under-credit by approach')
    fig.tight_layout()
    fig.savefig(SCRATCH / 'three_approaches_class_over_credit.png', dpi=120, bbox_inches='tight')
    log(f'Wrote {SCRATCH / "three_approaches_class_over_credit.png"}')


if __name__ == '__main__':
    main()
