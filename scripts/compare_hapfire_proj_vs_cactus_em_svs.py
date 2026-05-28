"""hapFIRE-projection vs cactus_em on site04 SVs (and SNPs as control).

Input pipeline for the hapFIRE-projection branch:
  - per-sample ecotype_frequency.txt (231-vec) from xwu's hapFIRE outputs
  - cn_var_231_v2 (231 founders × N_records) — same panel hapFIRE was run on
  - project: hapfire_proj_af[r] = h · cn_var_231_v2[:, r]

Compared against cactus_em (window-mode, 200 kb) per-record AFs from
results/site04_231_v2/.

Stratification:
  - by var_type (SNP, INS, DEL)
  - by SV size class for non-SNPs
  - by whether the SV's window has any hapFIRE-panel SNP nearby (proxy for LD-tagged)

Headline question: do hapFIRE-projection (genome-wide h) and cactus_em
(per-window h) produce different SV AFs? If they agree, the cactus_em-specific
work isn't needed.
"""
from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.sparse import load_npz


def metrics(t, p, label=''):
    t, p = np.asarray(t), np.asarray(p)
    err = p - t
    rmse = float(np.sqrt(np.mean(err**2)))
    r = float(np.corrcoef(t, p)[0, 1]) if np.std(t) > 0 and np.std(p) > 0 else float('nan')
    slope, intercept = np.polyfit(t, p, 1) if np.std(t) > 0 else (float('nan'), float('nan'))
    return {'label': label, 'n': len(t), 'MAE': float(np.mean(np.abs(err))), 'RMSE': rmse,
            'R2': r * r, 'slope': float(slope), 'intercept': float(intercept)}


def main():
    _root = Path(__file__).resolve().parents[2]
    ap = argparse.ArgumentParser()
    ap.add_argument('--cn-var', default=str(_root / 'data/cn_var_231_v2.cn_var.npz'))
    ap.add_argument('--cn-var-meta', default=str(_root / 'data/cn_var_231_v2.meta.npz'))
    ap.add_argument('--hapfire-ecotype-dir', default='/global/scratch/projects/fc_moilab/projects/grenenet-phase1/frequency/hapFIRE_frequencies/samples/ecotype_frequency')
    ap.add_argument('--cactus-em-dir', default=str(_root / 'results/site04_231_v2'))
    ap.add_argument('--samples', nargs='*', default=None,
                    help='sample IDs; default = all MLFH04* present in both dirs')
    ap.add_argument('--out', default=str(_root / 'results/hapfire_proj_vs_cem_svs_summary.tsv'))
    args = ap.parse_args()

    # Load panel
    print('Loading cn_var + meta...', flush=True)
    cn = load_npz(args.cn_var)            # 231 × N records (sparse 0/1)
    meta = np.load(args.cn_var_meta, allow_pickle=True)
    founders = [str(f) for f in meta['founders']]
    rec_chrom = np.asarray(meta['chrom']).astype(str)
    rec_pos = np.asarray(meta['pos']).astype(np.int64)
    ref_len = np.asarray(meta['ref_len']).astype(np.int32)
    alt_len = np.asarray(meta['alt_len']).astype(np.int32)
    N = len(rec_chrom)
    print(f'  cn_var: {cn.shape}, founders={len(founders)}, records={N:,}', flush=True)

    # Find samples present in both dirs
    cem_dir = Path(args.cactus_em_dir)
    hf_dir = Path(args.hapfire_ecotype_dir)
    if args.samples:
        samples = args.samples
    else:
        cem_samples = {p.stem for p in cem_dir.glob('MLFH04*.tsv')}
        hf_samples = {p.name.replace('_ecotype_frequency.txt', '')
                      for p in hf_dir.glob('MLFH04*_ecotype_frequency.txt')}
        samples = sorted(cem_samples & hf_samples)
    print(f'  samples to compare: {len(samples)} (e.g., {samples[:3]})', flush=True)

    # var_type stratum
    var_type = np.where((ref_len == 1) & (alt_len == 1), 'SNP',
               np.where(alt_len > ref_len, 'INS',
               np.where(ref_len > alt_len, 'DEL', 'COMPLEX')))
    sv_size = np.maximum(ref_len, alt_len) - np.minimum(ref_len, alt_len)  # bp difference
    is_sv = var_type != 'SNP'
    print(f'  records by class: SNP={int((var_type=="SNP").sum()):,}  '
          f'INS={int((var_type=="INS").sum()):,}  DEL={int((var_type=="DEL").sum()):,}  '
          f'COMPLEX={int((var_type=="COMPLEX").sum()):,}', flush=True)

    cn_csc = cn.tocsc().astype(np.float32)
    rows = []
    for sid in samples:
        # Load hapFIRE per-sample ecotype freq (231-vec, ordered by ecotype id)
        hf = pd.read_csv(hf_dir / f'{sid}_ecotype_frequency.txt', sep='\t',
                         header=None, names=['ecotype', 'freq'])
        hf['ecotype'] = hf['ecotype'].astype(str)
        # Re-order to match cn_var founder order
        hf_dict = dict(zip(hf['ecotype'], hf['freq'].astype(np.float32)))
        h = np.array([hf_dict.get(f, 0.0) for f in founders], dtype=np.float32)
        h_sum = float(h.sum())
        # Project: hapfire_af[r] = h @ cn[:, r]
        hapfire_af = (h @ cn_csc).A1.astype(np.float32) if hasattr((h @ cn_csc), 'A1') \
                     else np.asarray(h @ cn_csc).flatten().astype(np.float32)

        # Load cactus_em per-record AFs (row-aligned with cn_var)
        cem = pd.read_csv(cem_dir / f'{sid}.tsv', sep='\t')
        cem_af = cem['alt_freq'].to_numpy().astype(np.float32)

        # Subset to polymorphic records (drop both-zero or both-1 to avoid trivial agreement)
        polymask = (cem_af > 0.001) & (cem_af < 0.999) & np.isfinite(cem_af) & np.isfinite(hapfire_af)
        for stratum_name, stratum_mask in [
            ('all_poly', polymask),
            ('SNP', polymask & (var_type == 'SNP')),
            ('INS', polymask & (var_type == 'INS')),
            ('DEL', polymask & (var_type == 'DEL')),
            ('SV_all', polymask & is_sv),
            ('SV_small_<50bp', polymask & is_sv & (sv_size < 50)),
            ('SV_med_50-500bp', polymask & is_sv & (sv_size >= 50) & (sv_size < 500)),
            ('SV_big_500-5kb', polymask & is_sv & (sv_size >= 500) & (sv_size < 5000)),
            ('SV_huge_>=5kb', polymask & is_sv & (sv_size >= 5000)),
        ]:
            n = int(stratum_mask.sum())
            if n < 5:
                continue
            m = metrics(cem_af[stratum_mask], hapfire_af[stratum_mask], label=stratum_name)
            m['sample'] = sid
            m['h_sum'] = h_sum
            rows.append(m)

    out = pd.DataFrame(rows)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out, sep='\t', index=False)

    print(f'\n=== aggregate across {len(samples)} samples ===')
    agg = out.groupby('label').agg(
        n_samples=('sample', 'nunique'),
        mean_n=('n', 'mean'),
        mean_R2=('R2', 'mean'),
        std_R2=('R2', 'std'),
        mean_RMSE=('RMSE', 'mean'),
        mean_MAE=('MAE', 'mean'),
        mean_slope=('slope', 'mean'),
    ).round(4)
    print(agg)
    print(f'\nwrote {args.out}')


if __name__ == '__main__':
    main()
