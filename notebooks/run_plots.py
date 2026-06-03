"""Generate all plots from results_comparison.ipynb as PNG files into ../plots/

Run with: python run_plots.py
"""
import os, sys, glob, json
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.sparse import load_npz

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / 'poolfreq' / 'data'
RES  = ROOT / 'poolfreq' / 'results'
PLOTS = ROOT / 'plots'
PLOTS.mkdir(exist_ok=True)
plt.rcParams.update({'figure.dpi': 110, 'savefig.dpi': 200, 'font.size': 10,
                     'axes.spines.top': False, 'axes.spines.right': False})


def metrics(y_pred, y_true):
    y_pred = np.asarray(y_pred); y_true = np.asarray(y_true)
    rmse = float(np.sqrt(np.mean((y_pred-y_true)**2)))
    if y_true.std()>1e-12:
        r2 = 1 - ((y_true-y_pred)**2).sum() / ((y_true-y_true.mean())**2).sum()
    else: r2 = float('nan')
    r = float(np.corrcoef(y_pred, y_true)[0,1])
    return r2, rmse, r


def project_recipe(cn_var, founders, panel_map_path):
    if panel_map_path:
        pm = pd.read_csv(panel_map_path, sep='\t')
        m = dict(zip(pm.Assembly_ID.astype(str), pm.Accession_ID.astype(str)))
        ids = [m.get(str(f)) for f in founders]
    else:
        ids = [str(f) for f in founders]
    rec = pd.read_csv(ROOT/'data'/'seedmix_recipe_normalized.tsv', sep='\t')
    rd = dict(zip(rec.ID.astype(str), rec.seed_prop))
    h = np.array([rd.get(i, 0.0) if i else 0.0 for i in ids])
    if h.sum() > 0: h = h / h.sum()
    af = (h.astype(np.float32) @ cn_var.toarray()).astype(np.float64)
    return af, h


# === Plot 1: cactus capped vs nocap SV size distribution ===
print("[1/6] cactus capped vs nocap SV sizes")
labels  = ['snp+1bp', '1-50bp', '50-200bp', '200bp-1kb', '1-10kb', '10-100kb', '100kb-1Mb']
capped  = [3608497, 797061, 17383, 14152, 11885, 1413, 20]
nocap   = [3608172, 797090, 17394, 14151, 11878, 1412, 20]
fig, axes = plt.subplots(1, 2, figsize=(11, 4))
x = np.arange(len(labels)); w = 0.4
axes[0].bar(x-w/2, capped, w, label='capped (Mar 19)', color='#1f77b4')
axes[0].bar(x+w/2, nocap,  w, label='nocap (Apr 29)',  color='#ff7f0e')
axes[0].set_yscale('log')
axes[0].set_xticks(x); axes[0].set_xticklabels(labels, rotation=30, ha='right')
axes[0].set_ylabel('# variants (log)'); axes[0].set_title('Cactus VCF: SV size distribution')
axes[0].legend()
delta = np.array(nocap) - np.array(capped)
axes[1].bar(x, delta, color=['#1f77b4' if d<0 else '#ff7f0e' for d in delta])
axes[1].axhline(0, color='black', linewidth=0.5)
axes[1].set_xticks(x); axes[1].set_xticklabels(labels, rotation=30, ha='right')
axes[1].set_ylabel('Δ records (nocap - capped)')
axes[1].set_title(f'Per-bin difference (total Δ = {sum(delta):+d})')
plt.tight_layout()
plt.savefig(PLOTS / 'cactus_capped_vs_nocap_sv_sizes.png', bbox_inches='tight')
plt.close()


# === Plot 2: SEEDMIX panel × mode comparison ===
print("[2/6] SEEDMIX panel × mode")
panels = {
    '82-founder': dict(cn_var=DATA/'cn_var_82.cn_var.npz',
                       meta=DATA/'cn_var_82.meta.npz',
                       panel_map=ROOT/'data'/'sv_panel_to_accession_id.tsv',
                       results_dir=RES/'seedmix_82'),
    '231-founder': dict(cn_var=DATA/'cn_var_231.cn_var.npz',
                        meta=DATA/'cn_var_231.meta.npz',
                        panel_map=None,
                        results_dir=RES/'seedmix_231'),
}

rows = []
for panel_name, p in panels.items():
    cn_var = load_npz(str(p['cn_var']))
    meta = np.load(p['meta'], allow_pickle=True)
    af_truth, _ = project_recipe(cn_var, list(meta['founders']), p['panel_map'])
    for tsv in sorted(p['results_dir'].glob('SEEDMIX_S*.tsv')):
        df = pd.read_csv(tsv, sep='\t')
        if len(df) != len(af_truth): continue
        af_pred = df.alt_freq.to_numpy()
        r2, rmse, r = metrics(af_pred, af_truth)
        row = {'panel': panel_name, 'sample': tsv.stem, 'R2_raw': r2, 'r_raw': r, 'RMSE_raw': rmse}
        if panel_name == '231-founder':
            a, b = np.polyfit(af_truth, af_pred, 1)
            af_cal = (af_pred - b) / a
            r2c, rmsec, rc = metrics(af_cal, af_truth)
            row.update({'R2_cal': r2c, 'r_cal': rc, 'RMSE_cal': rmsec, 'slope': a, 'intercept': b})
        rows.append(row)
summary = pd.DataFrame(rows)

fig, ax = plt.subplots(figsize=(8, 4))
groups = ['82-founder\nraw', '231-founder\nraw', '231-founder\ncalibrated']
vals = [
    summary.loc[summary.panel=='82-founder', 'R2_raw'].dropna().values,
    summary.loc[summary.panel=='231-founder', 'R2_raw'].dropna().values,
    summary.loc[summary.panel=='231-founder', 'R2_cal'].dropna().values,
]
colors = ['#1f77b4', '#d62728', '#2ca02c']
for i, (label, v, c) in enumerate(zip(groups, vals, colors)):
    ax.bar(i, v.mean(), color=c, alpha=0.6, edgecolor='black')
    ax.scatter([i]*len(v), v, color=c, edgecolor='black', s=40, zorder=3)
    ax.text(i, v.mean()+0.015, f'{v.mean():.3f}±{v.std():.3f}', ha='center', fontsize=9)
ax.set_xticks(range(3)); ax.set_xticklabels(groups)
ax.set_ylabel('R² (per-record alt_freq vs recipe truth)')
ax.set_ylim(0, 1.02)
ax.set_title('SEEDMIX × 8 replicates — recipe-truth recovery')
plt.tight_layout()
plt.savefig(PLOTS / 'seedmix_panel_mode_comparison.png', bbox_inches='tight')
plt.close()


# === Plot 3: 231-founder calibration scatter ===
print("[3/6] 231-founder calibration scatter")
p = panels['231-founder']
cn_var = load_npz(str(p['cn_var']))
meta = np.load(p['meta'], allow_pickle=True)
af_truth, _ = project_recipe(cn_var, list(meta['founders']), p['panel_map'])
df = pd.read_csv(p['results_dir']/'SEEDMIX_S1.tsv', sep='\t')
af_pred = df.alt_freq.to_numpy()
a, b = np.polyfit(af_truth, af_pred, 1)
af_cal = (af_pred - b) / a
rng = np.random.default_rng(0)
idx = rng.choice(len(af_truth), 50_000, replace=False)
fig, axes = plt.subplots(1, 2, figsize=(11, 5), sharex=True, sharey=True)
for ax, y, title in [(axes[0], af_pred[idx], f'Raw  (slope = {a:.3f})'),
                      (axes[1], af_cal[idx], 'Calibrated')]:
    ax.scatter(af_truth[idx], y, s=2, alpha=0.15, color='#1f77b4')
    ax.plot([0, 1], [0, 1], 'r-', linewidth=1, alpha=0.7)
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    ax.set_xlabel('truth alt_freq (recipe @ cn_var_231)')
    ax.set_title(title)
    r2, rmse, r = metrics(y, af_truth[idx])
    ax.text(0.05, 0.92, f'R² = {r2:.3f}\nRMSE = {rmse:.3f}', transform=ax.transAxes,
            fontsize=11, va='top', bbox=dict(facecolor='white', edgecolor='gray', boxstyle='round,pad=0.4'))
axes[0].set_ylabel('predicted alt_freq')
fig.suptitle('231-founder window mode: SEEDMIX_S1, scale-bias correction')
plt.tight_layout()
plt.savefig(PLOTS / 'calibration_scatter_231_seedmix_S1.png', bbox_inches='tight')
plt.close()


# === Plot 4: AC-stratified accuracy ===
print("[4/6] AC-stratified accuracy")
ac_bins = [('AC=1', 1, 1), ('AC 2-4', 2, 4), ('AC 5-10', 5, 10), ('AC>10', 11, 999_999)]
stratified = []
for panel_name, p in panels.items():
    cn_var = load_npz(str(p['cn_var']))
    meta = np.load(p['meta'], allow_pickle=True)
    af_truth, _ = project_recipe(cn_var, list(meta['founders']), p['panel_map'])
    ac = np.asarray(cn_var.sum(axis=0)).flatten()
    for tsv in sorted(p['results_dir'].glob('SEEDMIX_S*.tsv')):
        df = pd.read_csv(tsv, sep='\t')
        if len(df) != len(af_truth): continue
        af_pred = df.alt_freq.to_numpy()
        if panel_name == '231-founder':
            a, b = np.polyfit(af_truth, af_pred, 1)
            af_use = (af_pred - b) / a
            label = '231-founder calibrated'
        else:
            af_use = af_pred
            label = '82-founder raw'
        for bn, lo, hi in ac_bins:
            m = (ac>=lo) & (ac<=hi)
            if m.sum() < 100: continue
            r2, rmse, r = metrics(af_use[m], af_truth[m])
            stratified.append({'panel': label, 'sample': tsv.stem, 'ac_bin': bn,
                              'r2': r2, 'rmse': rmse, 'r': r})
df_strat = pd.DataFrame(stratified)

fig, axes = plt.subplots(1, 2, figsize=(12, 4))
x = np.arange(len(ac_bins)); labels_ac = [b[0] for b in ac_bins]
panels_to_plot = ['82-founder raw', '231-founder calibrated']
colors = ['#1f77b4', '#2ca02c']
for ax, metric in zip(axes, ['r2', 'r']):
    w = 0.35
    for i, (panel, color) in enumerate(zip(panels_to_plot, colors)):
        sub = df_strat[df_strat.panel == panel]
        means = [sub[sub.ac_bin==b][metric].mean() for b in labels_ac]
        stds  = [sub[sub.ac_bin==b][metric].std()  for b in labels_ac]
        ax.bar(x + i*w - w/2, means, w, yerr=stds, capsize=3, color=color, label=panel, alpha=0.7, edgecolor='black')
    ax.set_xticks(x); ax.set_xticklabels(labels_ac)
    ax.legend()
axes[0].set_ylabel('R²'); axes[0].set_title('AC-stratified R² (8-rep mean ± std)')
axes[0].set_ylim(-3, 1.05); axes[0].axhline(0, color='black', linewidth=0.4)
axes[1].set_ylabel('Pearson r'); axes[1].set_title('AC-stratified Pearson r')
axes[1].set_ylim(0, 1.05)
plt.tight_layout()
plt.savefig(PLOTS / 'ac_stratified_accuracy.png', bbox_inches='tight')
plt.close()


# === Plot 5: Tier 1 cactus_em sim ===
print("[5/6] Tier 1 cactus_em sim")
TIER1 = '/global/home/users/tbellg/scratch/visor_freqk/results/cactus_em/cactus_pos10mb/var/cov30/var_del_1kb_n10_f50_err001'
CHOSEN = ['100042','100043','100047','100053','100108','100130','100158','100208','100226','100282']
cn_var_82 = load_npz(str(DATA/'cn_var_82.cn_var.npz'))
meta82 = np.load(str(DATA/'cn_var_82.meta.npz'), allow_pickle=True)
founders_82 = list(meta82['founders'])
h_truth = np.zeros(82, dtype=np.float32)
for f in CHOSEN: h_truth[founders_82.index(f)] = 1/10
af_truth = (h_truth @ cn_var_82.toarray()).astype(np.float64)
chr1 = meta82['chrom'] == 'Chr1'
out_del = chr1 & ~((meta82['pos']>=10_000_000) & (meta82['pos']<10_001_000))
in_del  = chr1 &  ((meta82['pos']>=10_000_000) & (meta82['pos']<10_001_000))

fig, axes = plt.subplots(1, 2, figsize=(11, 5), sharex=True, sharey=True)
for ax, mode in zip(axes, ['global', 'window']):
    df = pd.read_csv(f'{TIER1}/cactus_em_{mode}.tsv', sep='\t')
    af = df.alt_freq.to_numpy()
    rng = np.random.default_rng(1)
    sub = rng.choice(np.flatnonzero(out_del), 30_000, replace=False)
    ax.scatter(af_truth[sub], af[sub], s=3, alpha=0.2, color='#1f77b4', label='outside DEL')
    ax.scatter(af_truth[in_del], af[in_del], s=15, alpha=0.6, color='#d62728', edgecolor='black', label='inside DEL (60 records)')
    ax.plot([0, 1], [0, 1], 'k-', lw=1, alpha=0.7)
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    ax.set_xlabel('truth alt_freq (uniform 1/10)')
    r2, rmse, r = metrics(af[out_del], af_truth[out_del])
    ax.set_title(f'cactus_em {mode} mode\nR²={r2:.4f}, RMSE={rmse:.4f}')
    ax.legend(loc='lower right', fontsize=8)
axes[0].set_ylabel('predicted alt_freq')
fig.suptitle('Tier 1: 10-founder pool + 1kb deletion @ Chr1:10M, f=0.5')
plt.tight_layout()
plt.savefig(PLOTS / 'tier1_cactus_em_global_vs_window.png', bbox_inches='tight')
plt.close()


# === Plot 6: Calibration slope stability ===
print("[6/6] Calibration slope stability")
df231 = summary[summary.panel == '231-founder'].dropna(subset=['slope']).reset_index(drop=True)
fig, ax = plt.subplots(figsize=(7, 4))
ax.scatter(range(len(df231)), df231.slope, s=80, color='#2ca02c', edgecolor='black')
ax.axhline(df231.slope.mean(), color='red', linestyle='--',
           label=f'mean = {df231.slope.mean():.4f} ± {df231.slope.std():.4f}')
ax.set_xticks(range(len(df231)))
ax.set_xticklabels([s.replace('SEEDMIX_', '') for s in df231['sample']])
ax.set_xlabel('SEEDMIX replicate'); ax.set_ylabel('calibration slope')
ax.set_title('231-founder calibration slope across SEEDMIX replicates\n(panel-intrinsic; std = 0.4%)')
ax.legend(); ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig(PLOTS / 'calibration_slope_stability.png', bbox_inches='tight')
plt.close()

print("\nDone. Plots in:", PLOTS)
for p in sorted(PLOTS.glob('*.png')):
    print(f"  {p.name}  ({p.stat().st_size/1024:.0f} kB)")
