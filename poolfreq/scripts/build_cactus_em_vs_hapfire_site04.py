#!/usr/bin/env python
"""Mirror of freqk_gr's build_hapfire_comparison: cactus_em vs hapFIRE on site04.

For each cactus_em output in results/site04_231_v2/MLFH04*.tsv:
  - Restrict to SNPs (ref_len=1, alt_len=1), drop duplicate (chrom,pos)
  - Match against hapfire's per-sample SNP frequency file
  - Compute per-sample r, R2, RMSE, MAE
Output: results/site04_231_v2/cactus_em_vs_hapfire_per_sample_corr.tsv (mirror of freqk_gr).
"""
from __future__ import annotations
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path('/carnegie/nobackup/scratch/tbellagio/hapfire_sv/poolfreq/results/site04_231_v2')
HAPFIRE_DIR = Path('/carnegie/nobackup/scratch/xwu/GrENE_net/hapFIRE_frequencies/samples/snp_frequency')
FREQK_CORR = Path('/carnegie/nobackup/scratch/tbellagio/freqk_gr/results/hapfire_vs_freqk_AF_per_sample_corr.tsv')


def main():
    samples = sorted([p.stem for p in ROOT.glob('MLFH04*.tsv')])
    print(f'cactus_em samples available: {len(samples)}')
    if not samples:
        print('No cactus_em outputs yet — run site04 array first.')
        return

    rows = []
    for sid in samples:
        cem_path = ROOT / f'{sid}.tsv'
        hf_path  = HAPFIRE_DIR / f'{sid}_snp_frequency.txt'
        if not hf_path.exists():
            print(f'  {sid}: no hapfire file, skipping')
            continue
        cem = pd.read_csv(cem_path, sep='\t')
        cem_snp = cem[(cem.ref_len==1) & (cem.alt_len==1)].drop_duplicates(['chrom','pos']).copy()
        cem_snp['chrom'] = cem_snp['chrom'].astype(str)

        hf = pd.read_csv(hf_path, sep='\t', header=None,
                         names=['chrom','pos','hapfire_p'],
                         dtype={'chrom':str, 'pos':int, 'hapfire_p':float})
        hf['chrom'] = 'Chr' + hf['chrom']

        m = cem_snp.merge(hf, on=['chrom','pos'], how='inner').dropna(subset=['alt_freq','hapfire_p'])
        if len(m) < 2:
            continue
        diff = m.alt_freq.to_numpy() - m.hapfire_p.to_numpy()
        r = float(np.corrcoef(m.alt_freq, m.hapfire_p)[0,1])
        rows.append({'sample_id': sid, 'n_snps': len(m),
                     'pearson_r': r, 'r_squared': r*r,
                     'rmse': float(np.sqrt(np.mean(diff**2))),
                     'mae': float(np.mean(np.abs(diff)))})
        print(f'  {sid}: n={len(m):,}  r={r:.4f}  R2={r*r:.4f}  RMSE={float(np.sqrt(np.mean(diff**2))):.4f}')

    out = ROOT / 'cactus_em_vs_hapfire_per_sample_corr.tsv'
    df_cem = pd.DataFrame(rows)
    df_cem.to_csv(out, sep='\t', index=False)
    print(f'\nWrote {out}  (n_samples={len(df_cem)})')

    # Side-by-side with freqk benchmark
    if FREQK_CORR.exists():
        fq = pd.read_csv(FREQK_CORR, sep='\t')
        fq_site04 = fq[fq.sample_id.str.startswith('MLFH04')]
        merged = df_cem.merge(fq_site04, on='sample_id', suffixes=('_cem','_freqk'))
        print(f'\nSummary across {len(merged)} site04 samples (mean ± std):')
        for col_cem, col_fq, label in [
            ('r_squared_cem','r_squared_freqk','R²'),
            ('rmse_cem','rmse_freqk','RMSE'),
            ('mae_cem','mae_freqk','MAE'),
            ('n_snps_cem','n_snps_freqk','n_snps'),
        ]:
            cem_mean = merged[col_cem].mean(); cem_std = merged[col_cem].std()
            fq_mean = merged[col_fq].mean();   fq_std = merged[col_fq].std()
            print(f'  {label:<8s}  cactus_em: {cem_mean:.4f} ± {cem_std:.4f}   freqk: {fq_mean:.4f} ± {fq_std:.4f}')

        out_merged = ROOT / 'cactus_em_vs_freqk_vs_hapfire_per_sample.tsv'
        merged.to_csv(out_merged, sep='\t', index=False)
        print(f'Wrote {out_merged}')


if __name__ == '__main__':
    main()
