"""hapFIRE-perblock-projection vs cactus_em on SVs, stratified by LD-tagging strength.

Hypothesis: cactus_em should win specifically on SVs that are POORLY TAGGED
by panel SNPs (low max_r²). For SVs in good LD with some panel SNP, hapFIRE's
per-block h gives a clean projection. For untagged SVs, cactus_em's k-mer-direct
evidence is supposed to add value.

Inputs:
  - per-sample hapFIRE per-block projection TSV (project_hapfire_perblock_to_records.py output)
  - per-sample cactus_em window-mode TSV (results/site04_231_v2/<sample>.tsv)
  - SV LD-tagging table: ld_sv_vs_xwu_snps_per_sv.tsv.gz with cols:
        chrom, sv_pos, sv_class, sv_size, sv_ac, max_r2_in_xwu, ...

Output: per-stratum metrics (hapFIRE vs cactus_em agreement) — disagreement is
the proxy for "would they make different SV-AF calls". A neat presentation:
across LD-tag bins (untagged / weak / strong), do the methods agree more or less?
"""
from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np
import pandas as pd


def metrics(t, p, label=''):
    t, p = np.asarray(t, dtype=np.float64), np.asarray(p, dtype=np.float64)
    finite = np.isfinite(t) & np.isfinite(p)
    t, p = t[finite], p[finite]
    if len(t) < 2:
        return {'label': label, 'n': int(finite.sum()), 'MAE': float('nan'),
                'RMSE': float('nan'), 'R2': float('nan'),
                'slope': float('nan'), 'intercept': float('nan')}
    err = p - t
    rmse = float(np.sqrt(np.mean(err**2)))
    r = float(np.corrcoef(t, p)[0, 1]) if np.std(t) > 0 and np.std(p) > 0 else float('nan')
    try:
        slope, intercept = np.polyfit(t, p, 1)
    except (np.linalg.LinAlgError, ValueError):
        slope, intercept = float('nan'), float('nan')
    return {'label': label, 'n': len(t), 'MAE': float(np.mean(np.abs(err))), 'RMSE': rmse,
            'R2': r * r, 'slope': float(slope), 'intercept': float(intercept)}


def main():
    _root = Path(__file__).resolve().parents[2]
    ap = argparse.ArgumentParser()
    ap.add_argument('--hapfire-proj-dir', default=str(_root / 'results/site04_hapfire_perblock_proj'))
    ap.add_argument('--cactus-em-dir', default=str(_root / 'results/site04_231_v2'))
    ap.add_argument('--ld-table', default='/global/scratch/users/tbellg/freqk_gr/ld_sv_snp/data/ld_sv_vs_xwu_snps_per_sv.tsv.gz')
    ap.add_argument('--samples', nargs='*', default=None,
                    help='sample IDs (default: all present in both dirs)')
    ap.add_argument('--out', default=str(_root / 'results/hapfire_perblock_vs_cem_by_ldtag.tsv'))
    args = ap.parse_args()

    print(f'Loading LD-tagging table from {args.ld_table}')
    ld = pd.read_csv(args.ld_table, sep='\t')
    ld['chrom'] = ld['chrom'].astype(str)
    print(f'  n_svs in table: {len(ld):,}')

    # Bin SVs by max LD-tag r²
    def ld_bin(r2):
        if pd.isna(r2):
            return 'no_snp_in_window'
        if r2 < 0.2:
            return '1.untagged_<0.2'
        if r2 < 0.5:
            return '2.weak_0.2-0.5'
        if r2 < 0.8:
            return '3.medium_0.5-0.8'
        return '4.strong_>=0.8'
    ld['ld_bin'] = ld['max_r2_in_xwu'].apply(ld_bin)
    print(f'  ld_bin counts:')
    for k, v in ld['ld_bin'].value_counts().sort_index().items():
        print(f'    {k}: {v:,}')

    hf_dir = Path(args.hapfire_proj_dir)
    cem_dir = Path(args.cactus_em_dir)
    if args.samples:
        samples = args.samples
    else:
        hf_samples = {p.stem for p in hf_dir.glob('MLFH04*.tsv')}
        cem_samples = {p.stem for p in cem_dir.glob('MLFH04*.tsv')}
        samples = sorted(hf_samples & cem_samples)
    print(f'  samples to compare: {len(samples)} (e.g., {samples[:3]})', flush=True)

    if not samples:
        raise SystemExit('no overlapping samples')

    rows = []
    for sid in samples:
        hf = pd.read_csv(hf_dir / f'{sid}.tsv', sep='\t')
        cem = pd.read_csv(cem_dir / f'{sid}.tsv', sep='\t')
        # Same row order (cn_var) — zip
        if len(hf) != len(cem):
            raise SystemExit(f'{sid}: row count mismatch ({len(hf)} vs {len(cem)})')
        df = hf[['chrom', 'pos', 'ref_len', 'alt_len']].copy()
        df['hapfire_af'] = hf['alt_freq'].to_numpy()
        df['cem_af'] = cem['alt_freq'].to_numpy()
        df['chrom'] = df['chrom'].astype(str)

        # Var type
        df['var_type'] = np.where((df.ref_len == 1) & (df.alt_len == 1), 'SNP',
                         np.where(df.alt_len > df.ref_len, 'INS',
                         np.where(df.ref_len > df.alt_len, 'DEL', 'COMPLEX')))
        # SVs only — join to LD-tag table on (chrom, pos)
        # ld_table's sv_pos may not match cn_var pos exactly; let's see
        sv_df = df[df.var_type != 'SNP'].rename(columns={'pos': 'sv_pos'})
        merged = sv_df.merge(ld[['chrom', 'sv_pos', 'sv_class', 'sv_size', 'sv_ac',
                                  'max_r2_in_xwu', 'ld_bin', 'n_xwu_tags_r2_0.2']],
                              on=['chrom', 'sv_pos'], how='inner')
        # Keep polymorphic only (drop both ~0 or ~1 to avoid trivial agreement)
        merged = merged[(merged.hapfire_af.between(0.001, 0.999, inclusive='neither')) |
                        (merged.cem_af.between(0.001, 0.999, inclusive='neither'))].copy()

        # Aggregate per-sample metrics by ld_bin
        for bin_label in sorted(merged['ld_bin'].unique()):
            sub = merged[merged.ld_bin == bin_label]
            if len(sub) < 5:
                continue
            m = metrics(sub.cem_af, sub.hapfire_af,
                        label=bin_label)
            m['sample'] = sid
            m['n_svs_in_bin'] = len(sub)
            rows.append(m)
        # Also overall SV agreement
        if len(merged) > 5:
            m = metrics(merged.cem_af, merged.hapfire_af, label='SV_all')
            m['sample'] = sid
            m['n_svs_in_bin'] = len(merged)
            rows.append(m)

    out_df = pd.DataFrame(rows)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    out_df.to_csv(args.out, sep='\t', index=False)

    # Aggregate
    print(f'\n=== aggregate across {len(samples)} samples (mean across samples per ld_bin) ===')
    agg = out_df.groupby('label').agg(
        n_samples=('sample', 'nunique'),
        mean_n_svs=('n_svs_in_bin', 'mean'),
        mean_R2=('R2', 'mean'),
        std_R2=('R2', 'std'),
        mean_RMSE=('RMSE', 'mean'),
        mean_MAE=('MAE', 'mean'),
        mean_slope=('slope', 'mean'),
    ).round(4)
    print(agg.to_string())
    print(f'\nwrote {args.out}')
    print(f'\ninterpretation:')
    print(f'  - if cem and hapfire agree at ALL ld_bins similarly → cem redundant')
    print(f'  - if agreement collapses in untagged_<0.2 bin → cem provides info hapfire cant')


if __name__ == '__main__':
    main()
