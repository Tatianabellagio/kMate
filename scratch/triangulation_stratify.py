"""
Stratified triangulation on SEEDMIX_S1 Chr1 single-ALT SNPs.
Loads the parquet from triangulation_audit.py and slices by three mechanisms.

For each stratum: report n, mean cem, mean hf, mean rec_lab, and the differences:
- hf - cem            (the headline disagreement)
- hf - rec_lab        (hapFIRE vs panel-side anchor)
- cem - rec_lab       (EM convergence error)
- rec_unif - rec_lab  (uniform-vs-lab anchor delta, mostly noise check)
"""
import numpy as np, pandas as pd

ROOT = '/carnegie/nobackup/scratch/tbellagio/hapfire_sv'
df = pd.read_parquet(f'{ROOT}/scratch/triangulation_S1_chr1.parquet')

m = (np.isfinite(df.cem) & np.isfinite(df.hf)
     & np.isfinite(df.rec_lab) & np.isfinite(df.rec_unif)
     & (df.cem.between(0,1)) & (df.hf.between(0,1))
     & (df.rec_lab.between(0,1)) & df.single_alt)
df = df[m].copy()
print(f'records analyzed: {len(df):,}')
df['hf_minus_cem']     = df.hf      - df.cem
df['hf_minus_reclab']  = df.hf      - df.rec_lab
df['cem_minus_reclab'] = df.cem     - df.rec_lab
df['recunif_minus_lab']= df.rec_unif- df.rec_lab

def stratify(df, col, bins, labels=None):
    if labels is None:
        labels = [f'[{a:g},{b:g})' for a,b in zip(bins[:-1], bins[1:])]
    df = df.copy()
    df['bin'] = pd.cut(df[col], bins=bins, labels=labels, include_lowest=True)
    out = df.groupby('bin', observed=True).agg(
        n=('hf','size'),
        mean_cem=('cem','mean'),
        mean_hf=('hf','mean'),
        mean_reclab=('rec_lab','mean'),
        hf_minus_cem=('hf_minus_cem','mean'),
        hf_minus_reclab=('hf_minus_reclab','mean'),
        cem_minus_reclab=('cem_minus_reclab','mean'),
        recunif_minus_lab=('recunif_minus_lab','mean'),
        mae_hf_cem=('hf_minus_cem', lambda x: float(np.mean(np.abs(x)))),
        out_pct=('hf_minus_cem', lambda x: float(100*np.mean(np.abs(x)>0.1))),
    )
    return out

def fmt(out):
    cols = ['n','mean_cem','mean_hf','mean_reclab',
            'hf_minus_cem','hf_minus_reclab','cem_minus_reclab',
            'recunif_minus_lab','mae_hf_cem','out_pct']
    s = out[cols].copy()
    for c in cols:
        if c == 'n':
            s[c] = s[c].map(lambda v: f'{int(v):>7,}')
        elif c == 'out_pct':
            s[c] = s[c].map(lambda v: f'{v:5.2f}%')
        else:
            s[c] = s[c].map(lambda v: f'{v:+.4f}')
    return s.to_string()

print('\n=== STRATIFY by F_MISSING in v3qc_v3 cn_var ===')
print(fmt(stratify(df, 'fmiss', [-0.0001, 0.01, 0.05, 0.1, 0.2, 0.5, 1.01])))

print('\n=== STRATIFY by pre-mask PG het rate (Mechanism 1: het-mask damage) ===')
df_m1 = df.dropna(subset=['pre_mask_het']).copy()
print(f'(records with het_rate available: {len(df_m1):,})')
print(fmt(stratify(df_m1, 'pre_mask_het', [-0.0001, 0.01, 0.05, 0.10, 0.20, 0.50, 1.01])))

print('\n=== STRATIFY by log2(carrier_pg / carrier_cactus) (Mechanism 2: bubble flattening / panel-side asymmetry) ===')
# Build asymmetry on AF scale (not log) too, but stratify on log to see direction
print(fmt(stratify(df, 'log2_pg_c', [-10, -2, -0.5, 0.5, 2, 10],
                   labels=['cactus-spec (<-2)','mild-cac (-2..-0.5)','balanced (|·|<.5)','mild-pg (0.5..2)','PG-spec (>2)'])))

print('\n=== STRATIFY by carrier_pg (panel-side PG carrier rate, intrinsic site frequency) ===')
print(fmt(stratify(df, 'carrier_pg', [-0.0001, 0.01, 0.05, 0.10, 0.25, 0.50, 1.01])))

# --- "Where does the concentration sit?" — look at sign of hf-cem by AF bin ---
print('\n=== Direction of hf-cem disagreement by AF level (rec_lab as x-axis) ===')
print(fmt(stratify(df, 'rec_lab', [-0.0001, 0.02, 0.05, 0.10, 0.25, 0.50, 0.99, 1.01],
                   labels=['<0.02','0.02-0.05','0.05-0.10','0.10-0.25','0.25-0.50','0.50-0.99','>=0.99'])))

# --- Outliers: where do |hf-cem|>0.1 sit on the AF axis, and what's their het_rate distribution? ---
print('\n=== OUTLIERS ONLY (|hf-cem|>0.1): joint distribution by panel-side AF and het_rate ===')
out = df[np.abs(df.hf_minus_cem)>0.1].copy()
print(f'n outliers: {len(out):,} ({100*len(out)/len(df):.2f}% of records)')
print(f'  mean rec_lab: {out.rec_lab.mean():.3f}  mean fmiss: {out.fmiss.mean():.3f}')
print(f'  mean pre_mask_het: {out.pre_mask_het.mean():.4f}  (vs population mean: {df.pre_mask_het.mean():.4f})')
print(f'  outliers by direction:')
print(f'    hf > cem (cactus-em under-predicts): {(out.hf_minus_cem>0).sum():,} ({100*(out.hf_minus_cem>0).mean():.1f}% of outliers)')
print(f'    hf < cem (cactus-em over-predicts):  {(out.hf_minus_cem<0).sum():,} ({100*(out.hf_minus_cem<0).mean():.1f}% of outliers)')
print(f'  mean hf for hf>cem outliers: {out[out.hf_minus_cem>0].hf.mean():.3f}, cem: {out[out.hf_minus_cem>0].cem.mean():.3f}, rec_lab: {out[out.hf_minus_cem>0].rec_lab.mean():.3f}')
print(f'  mean hf for hf<cem outliers: {out[out.hf_minus_cem<0].hf.mean():.3f}, cem: {out[out.hf_minus_cem<0].cem.mean():.3f}, rec_lab: {out[out.hf_minus_cem<0].rec_lab.mean():.3f}')

# Does rec_lab side with hapFIRE or cactus_em on each outlier?
print('\n=== AT OUTLIERS, WHO DOES THE PANEL-SIDE ANCHOR (rec_lab) AGREE WITH? ===')
side_with_hf  = (np.abs(out.rec_lab - out.hf)  < np.abs(out.rec_lab - out.cem))
side_with_cem = ~side_with_hf
print(f'  rec_lab closer to hapFIRE : {side_with_hf.sum():,} ({100*side_with_hf.mean():.1f}%)')
print(f'  rec_lab closer to cactus_em: {side_with_cem.sum():,} ({100*side_with_cem.mean():.1f}%)')
print(f'  Interpretation:')
print(f'    rec_lab is h_lab @ cn_var. It carries any cactus-panel construction error.')
print(f'    If it sides with cactus_em on outliers, the panel and EM agree — but hapFIRE')
print(f'    sees something different. Either panel-construction issue (mech 1+2) OR hapFIRE error.')
print(f'    If it sides with hapFIRE, the panel knows the truer answer and the EM is off.')
