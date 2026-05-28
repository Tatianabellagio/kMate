"""Cross-replicate analysis of SEEDMIX S1-S8: cactus_em (v3, dedup'd reads) vs hapFIRE.

Both methods estimate the same h per replicate (231 founders). The 8 replicates
are technical replicates of the same recipe (intended uniform 1/231). Cross-rep
variance is the noise floor; cross-method disagreement above the noise floor is
real methodological bias.

Inputs:
  cactus_em h: scratch/seedmix_v3_dedup/SEEDMIX_S{1..8}.h_per_chrom.npz (Chr1-5 each)
  hapFIRE h:   /global/scratch/projects/fc_moilab/projects/grenenet-phase1/frequency/hapFIRE_frequencies/seed_mix/s{1..8}_ecotype_frequency.txt
  panel split: data/founder_split_cactus_pg.json
  founder order: data/vcf_samples_231.txt

Outputs:
  scratch/cross_replicate_seedmix.tsv      — per-founder per-method mean+std across reps
  scratch/cross_replicate_metrics.tsv      — per-replicate, per-method summary stats
  scratch/cross_replicate_plots/*.png      — diagnostic figures
"""
from __future__ import annotations
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path('/carnegie/nobackup/scratch/tbellagio/kmate')
SAMPLES_FILE = ROOT / 'data/vcf_samples_231.txt'
SPLIT_FILE = ROOT / 'data/founder_split_cactus_pg.json'
CEM_DIR = ROOT / 'scratch/seedmix_v3_dedup'
HF_DIR = Path('/global/scratch/projects/fc_moilab/projects/grenenet-phase1/frequency/hapFIRE_frequencies/seed_mix')
OUT_DIR = ROOT / 'scratch/cross_replicate_plots'
OUT_DIR.mkdir(exist_ok=True)
TRUTH = 1.0 / 231


def log(*a, **kw): print(*a, **kw, flush=True)


def load_cactus_em(s, founders):
    """Load cactus_em h for replicate s; average per-chrom h vectors."""
    f = CEM_DIR / f'SEEDMIX_S{s}.h_per_chrom.npz'
    z = np.load(f, allow_pickle=True)
    h_chroms = []
    for chrom in ['Chr1', 'Chr2', 'Chr3', 'Chr4', 'Chr5']:
        if chrom in z:
            h_chroms.append(z[chrom].astype(np.float64))
    h_mean = np.mean(h_chroms, axis=0)
    h_mean = h_mean / h_mean.sum()  # defensive renormalize
    return h_mean


def load_hapfire(s, founders):
    """Load hapFIRE h for replicate s; ensure same founder order as cactus_em."""
    f = HF_DIR / f's{s}_ecotype_frequency.txt'
    df = pd.read_csv(f, sep='\t', header=None, names=['founder', 'h'])
    df['founder'] = df['founder'].astype(str)
    # Reorder to match `founders`
    fmap = dict(zip(df['founder'], df['h']))
    h = np.array([fmap.get(f, np.nan) for f in founders], dtype=np.float64)
    # Renormalize (hapFIRE outputs should sum to 1, but defensive)
    if np.isnan(h).any():
        log(f'  WARN S{s}: {np.isnan(h).sum()} founders missing in hapFIRE output')
        h[np.isnan(h)] = 0.0
    h = h / h.sum()
    return h


def per_replicate_stats(h, label, rep, is_cactus, is_pg):
    eff_n = float(1.0 / np.sum(h ** 2))
    l1 = float(np.sum(np.abs(h - TRUTH)))
    rmse = float(np.sqrt(np.mean((h - TRUTH) ** 2)))
    cv = float(h.std() / h.mean())
    sum_cac = float(h[is_cactus].sum())
    sum_pg = float(h[is_pg].sum())
    # truth class mass (uniform): cactus 80/231=0.3463, PG 151/231=0.6537
    return dict(method=label, replicate=rep,
                eff_n=eff_n, l1_to_uniform=l1, rmse_to_uniform=rmse, cv=cv,
                sum_cactus=sum_cac, sum_pg=sum_pg,
                over_cactus=sum_cac / 0.3463, over_pg=sum_pg / 0.6537)


def main():
    founders = [l.strip() for l in open(SAMPLES_FILE)]
    split = json.load(open(SPLIT_FILE))
    is_cactus = np.array([f in set(split['cactus']) for f in founders])
    is_pg = np.array([f in set(split['PG']) for f in founders])
    log(f'231 founders: {is_cactus.sum()} cactus + {is_pg.sum()} PG')

    # Load all 16 h vectors (8 reps × 2 methods)
    cem = np.zeros((8, 231))
    hf = np.zeros((8, 231))
    missing = []
    for s in range(1, 9):
        cem_f = CEM_DIR / f'SEEDMIX_S{s}.h_per_chrom.npz'
        if not cem_f.exists():
            missing.append(f'cactus_em S{s}')
            continue
        cem[s-1] = load_cactus_em(s, founders)
        hf[s-1] = load_hapfire(s, founders)
        log(f'  loaded S{s}: cem sum={cem[s-1].sum():.4f}, hapfire sum={hf[s-1].sum():.4f}')
    if missing:
        log(f'WARN missing: {missing}')

    # Per-replicate metrics
    rows = []
    for s in range(1, 9):
        rows.append(per_replicate_stats(cem[s-1], 'cactus_em', s, is_cactus, is_pg))
        rows.append(per_replicate_stats(hf[s-1], 'hapFIRE', s, is_cactus, is_pg))
    metrics = pd.DataFrame(rows)
    metrics.to_csv(ROOT / 'scratch/cross_replicate_metrics.tsv', sep='\t', index=False)

    log('\n=== Per-replicate metrics ===')
    log(metrics.to_string(index=False))

    log('\n=== Summary by method (mean ± std across 8 reps) ===')
    summary = metrics.groupby('method').agg(['mean', 'std']).round(4)
    log(summary.to_string())

    # Per-founder mean+std across reps
    cem_mean = cem.mean(axis=0); cem_std = cem.std(axis=0)
    hf_mean = hf.mean(axis=0);   hf_std = hf.std(axis=0)
    per_founder = pd.DataFrame({
        'founder': founders,
        'panel': np.where(is_cactus, 'cactus', np.where(is_pg, 'PG', 'unknown')),
        'cem_mean': cem_mean, 'cem_std': cem_std,
        'hf_mean': hf_mean,   'hf_std': hf_std,
        'agreement_pearson_r_across_reps': [
            float(np.corrcoef(cem[:, i], hf[:, i])[0, 1]) if cem[:, i].std() > 0 and hf[:, i].std() > 0
            else np.nan for i in range(231)],
        'systematic_diff_cem_minus_hf': cem_mean - hf_mean,
    })
    per_founder.to_csv(ROOT / 'scratch/cross_replicate_seedmix.tsv', sep='\t', index=False)

    # === Key cross-method comparison: are mean h_cem and mean h_hf in agreement? ===
    r_cem_hf_mean = float(np.corrcoef(cem_mean, hf_mean)[0, 1])
    l1_cem_vs_hf_mean = float(np.sum(np.abs(cem_mean - hf_mean)))
    log(f'\n=== Cross-method on mean (8-rep average h vectors) ===')
    log(f'  cactus_em mean h vs hapFIRE mean h: Pearson r = {r_cem_hf_mean:.4f}')
    log(f'  L1 distance between mean h vectors = {l1_cem_vs_hf_mean:.4f}')
    log(f'  L1 from cactus_em mean to uniform = {np.sum(np.abs(cem_mean - TRUTH)):.4f}')
    log(f'  L1 from hapFIRE   mean to uniform = {np.sum(np.abs(hf_mean - TRUTH)):.4f}')

    # Per-founder noise floor: std across reps. If two methods disagree by more
    # than sqrt(σ_cem² + σ_hf²), the disagreement is real, not noise.
    noise_floor = np.sqrt(cem_std ** 2 + hf_std ** 2)
    real_disagree = np.abs(cem_mean - hf_mean) > 2 * noise_floor
    log(f'\n  founders with |cem_mean - hf_mean| > 2·noise_floor: {real_disagree.sum()}/231 '
        f'({100*real_disagree.sum()/231:.1f}%)')
    log(f'    of which on cactus side: {(real_disagree & is_cactus).sum()}/80')
    log(f'    of which on PG side:     {(real_disagree & is_pg).sum()}/151')

    # === PLOTS ===
    # 1) eff_n + L1 boxplots per method
    fig, axes = plt.subplots(1, 2, figsize=(11, 5))
    methods = ['cactus_em', 'hapFIRE']
    eff_data = [metrics[metrics.method == m].eff_n.values for m in methods]
    axes[0].boxplot(eff_data, labels=methods)
    axes[0].axhline(231, color='k', linestyle='--', lw=0.8, label='truth (uniform=231)')
    axes[0].set_ylabel('eff_n_founders')
    axes[0].set_title('Effective # founders across 8 SEEDMIX reps')
    axes[0].legend()
    axes[0].grid(alpha=0.3, axis='y')

    l1_data = [metrics[metrics.method == m].l1_to_uniform.values for m in methods]
    axes[1].boxplot(l1_data, labels=methods)
    axes[1].set_ylabel('L1 distance from uniform 1/231')
    axes[1].set_title('Per-replicate L1 to uniform expectation')
    axes[1].grid(alpha=0.3, axis='y')
    fig.tight_layout()
    fig.savefig(OUT_DIR / 'concentration_per_replicate.png', dpi=120, bbox_inches='tight')
    log(f'Wrote {OUT_DIR / "concentration_per_replicate.png"}')

    # 2) Scatter: cem_mean vs hf_mean, colored by panel
    fig, ax = plt.subplots(figsize=(7, 7))
    ax.scatter(hf_mean[is_cactus], cem_mean[is_cactus], s=20, alpha=0.7,
                color='#d62728', label='cactus 80', edgecolors='none')
    ax.scatter(hf_mean[is_pg], cem_mean[is_pg], s=20, alpha=0.7,
                color='#1f77b4', label='PG 151', edgecolors='none')
    xmax = max(cem_mean.max(), hf_mean.max()) * 1.05
    ax.plot([0, xmax], [0, xmax], 'k--', lw=0.8, label='identity')
    ax.axhline(TRUTH, color='gray', linestyle=':', lw=0.5)
    ax.axvline(TRUTH, color='gray', linestyle=':', lw=0.5)
    ax.set_xlim(0, xmax); ax.set_ylim(0, xmax)
    ax.set_xlabel('hapFIRE mean h (across 8 reps)')
    ax.set_ylabel('cactus_em mean h (across 8 reps)')
    ax.set_title(f'Per-founder mean h (Pearson r = {r_cem_hf_mean:.3f})')
    ax.legend()
    fig.tight_layout()
    fig.savefig(OUT_DIR / 'cem_vs_hf_mean_scatter.png', dpi=120, bbox_inches='tight')
    log(f'Wrote {OUT_DIR / "cem_vs_hf_mean_scatter.png"}')

    # 3) Per-founder std (noise floor) per method
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.scatter(cem_std[is_cactus], hf_std[is_cactus], s=20, alpha=0.7,
                color='#d62728', label='cactus 80', edgecolors='none')
    ax.scatter(cem_std[is_pg], hf_std[is_pg], s=20, alpha=0.7,
                color='#1f77b4', label='PG 151', edgecolors='none')
    smax = max(cem_std.max(), hf_std.max()) * 1.05
    ax.plot([0, smax], [0, smax], 'k--', lw=0.8, label='identity')
    ax.set_xlim(0, smax); ax.set_ylim(0, smax)
    ax.set_xlabel('cactus_em per-founder σ across 8 reps')
    ax.set_ylabel('hapFIRE per-founder σ across 8 reps')
    ax.set_title('Per-founder noise floor (each point is one founder)')
    ax.legend()
    fig.tight_layout()
    fig.savefig(OUT_DIR / 'noise_floor_per_founder.png', dpi=120, bbox_inches='tight')
    log(f'Wrote {OUT_DIR / "noise_floor_per_founder.png"}')

    # 4) Rank-h plot like user's earlier figure, on the 8-rep MEAN
    fig, ax = plt.subplots(figsize=(11, 5))
    order = np.argsort(-cem_mean)  # descending by cem
    rank = np.arange(231)
    for cls, mask, color, label in [
        ('cactus', is_cactus, '#9ecae1', 'cactus_em (cactus 80)'),
        ('PG', is_pg, '#3182bd', 'cactus_em (PG 151)'),
    ]:
        idx = np.where(mask)[0]
        order_of_idx = np.array([np.where(order == i)[0][0] for i in idx])
        ax.scatter(order_of_idx, cem_mean[idx], s=15, alpha=0.5, color=color, label=label)
    ax.scatter(rank, hf_mean[order], s=15, alpha=0.5, color='#2ca02c', label='hapFIRE')
    ax.axhline(TRUTH, color='gray', linestyle='--', lw=0.8, label=f'1/231 ≈ {TRUTH:.5f}')
    ax.set_xlabel('Founder rank (sorted by cactus_em mean)')
    ax.set_ylabel('Mean h across 8 SEEDMIX reps')
    ax.set_title(f'Per-founder mean h across 8 dedup\'d SEEDMIX reps; cactus_em eff_n_mean='
                  f'{1/np.sum(cem_mean**2):.0f}, hapFIRE eff_n_mean={1/np.sum(hf_mean**2):.0f}')
    ax.legend()
    fig.tight_layout()
    fig.savefig(OUT_DIR / 'rank_h_mean_8rep.png', dpi=120, bbox_inches='tight')
    log(f'Wrote {OUT_DIR / "rank_h_mean_8rep.png"}')

    log('\nDone.')


if __name__ == '__main__':
    main()
