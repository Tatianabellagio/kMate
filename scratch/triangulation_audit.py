"""
Triangulation audit: hapFIRE vs cactus_em vs lab-recipe-projection, SEEDMIX_S1 Chr1.

What this script does:
- Builds two "anchors" alongside the two estimators:
    A. recipe_lab_proj  =  h_lab @ cn_var  /  h_lab @ cn_var_called
       where h_lab is the lab-documented seed_prop. Anchor for *what cactus_em should
       output IF the EM converged perfectly and the cactus panel genotyped correctly*.
       NOT an external anchor — bakes in mechanisms 1 + 2 (het-mask, bubble flattening).
    B. recipe_uniform   =  h_unif @ cn_var / h_unif @ cn_var_called   (TRUTH_H=1/231)
       This is what the existing PIPELINE_STATE notebook calls 'recipe FIXED'.
       Equivalent to "AC/AN over called founders, ignoring lab mix."
- Estimators:
    cactus_em-AF: the production v3qc_v3 mixed-loose result
    hapFIRE-AF:  xwu's per-SNP AF (joined by position on Chr1)
- Stratifications:
    1. F_MISSING in v3qc_v3 cn_var (proxy for panel sparsity / het-mask aftermath)
    2. pre-mask PG het rate (Mechanism 1: site-level het-mask damage)
    3. carrier-rate asymmetry log2(PG/cactus) (Mechanism 2: bubble-flattening direction)
"""
import json
import numpy as np
import pandas as pd
import scipy.sparse as sp

ROOT = '/carnegie/nobackup/scratch/tbellagio/kmate'
SCRATCH = f'{ROOT}/scratch'
HF = '/global/scratch/projects/fc_moilab/projects/grenenet-phase1/frequency/hapFIRE_frequencies/seed_mix/s1_snp_frequency.txt'

# --- load cn_var v3qc_v3 + called mask + meta ---
meta = np.load(f'{ROOT}/data/cn_var_231_v3qc_v3.meta.npz', allow_pickle=True)
founders = list(np.asarray(meta['founders']).astype(str))
F = len(founders)
chrom = np.asarray(meta['chrom']).astype(str)
pos   = np.asarray(meta['pos']).astype(np.int64)
rlen  = np.asarray(meta['ref_len']).astype(np.int32)
alen  = np.asarray(meta['alt_len']).astype(np.int32)

chr1_snp = (chrom == 'Chr1') & (rlen == 1) & (alen == 1)
chr1_snp_idx = np.where(chr1_snp)[0]
chr1_snp_pos = pos[chr1_snp_idx]
print(f'Chr1 SNP records: {len(chr1_snp_idx):,}')

cn_var = sp.load_npz(f'{ROOT}/data/cn_var_231_v3qc_v3.cn_var.npz')
cn_var_called = sp.load_npz(f'{ROOT}/data/cn_var_231_v3qc_v3.cn_var_called.npz')

cv = cn_var[:, chr1_snp_idx].toarray().astype(np.float32)         # (F, N_chr1_snp)
cvc = cn_var_called[:, chr1_snp_idx].toarray().astype(np.float32) # (F, N_chr1_snp)

# --- lab recipe h ---
recipe_df = pd.read_csv(f'{ROOT}/data/seedmix_recipe_normalized.tsv', sep='\t', dtype={'ID':str})
recipe_map = dict(zip(recipe_df.ID.astype(str), recipe_df.seed_prop))
h_lab = np.array([recipe_map.get(f, 0.0) for f in founders], dtype=np.float64)
print(f'h_lab sum = {h_lab.sum():.6f}  min={h_lab.min():.5f}  max={h_lab.max():.5f}')

h_unif = np.full(F, 1.0/F, dtype=np.float64)

def project(h, cv, cvc):
    num = h @ cv
    den = h @ cvc
    af = np.where(den > 1e-12, num / np.maximum(den, 1e-12), np.nan)
    return af.astype(np.float32), den.astype(np.float32)

recipe_lab,    den_lab    = project(h_lab,  cv, cvc)
recipe_unif,   den_unif   = project(h_unif, cv, cvc)

# --- cactus_em production AF ---
cem_df = pd.read_csv(f'{SCRATCH}/v3qc_v3_mixedloose_chr1/SEEDMIX_S1.tsv', sep='\t')
cem_af_full = cem_df['alt_freq'].to_numpy(np.float32)
cem_af = cem_af_full[chr1_snp_idx]
print(f'cactus_em AF: n_finite={np.isfinite(cem_af).sum():,} / {len(cem_af):,}')

# --- hapFIRE AF (xwu) — position-joined on Chr1 ---
hf = pd.read_csv(HF, sep='\t', header=None, names=['chrom','pos','af'])
hf = hf[hf.chrom.astype(str) == '1']
hf_d = dict(zip(hf.pos.astype(np.int64), hf.af.astype(np.float32)))
hf_af = np.array([hf_d.get(p, np.nan) for p in chr1_snp_pos], dtype=np.float32)
print(f'hapFIRE aligned: n_finite={np.isfinite(hf_af).sum():,} / {len(hf_af):,}')

# Drop multi-allelic positions for the hapFIRE join (single-AF-per-pos in xwu file)
from collections import Counter
pos_counts = Counter(chr1_snp_pos.tolist())
single_alt = np.array([pos_counts[p] == 1 for p in chr1_snp_pos])
print(f'Single-ALT positions (for hapFIRE join): {single_alt.sum():,} / {len(chr1_snp_pos):,}')

# --- carrier asymmetry (mechanism 2 proxy) ---
with open(f'{ROOT}/data/founder_split_cactus_pg.json') as fh:
    split = json.load(fh)
cactus_set = set(map(str, split['cactus']))
pg_set     = set(map(str, split['PG']))
is_cactus = np.array([f in cactus_set for f in founders])
is_pg     = np.array([f in pg_set     for f in founders])
print(f'cactus founders: {is_cactus.sum()}, PG founders: {is_pg.sum()}')

carrier_cactus = cv[is_cactus].sum(axis=0) / np.maximum(cvc[is_cactus].sum(axis=0), 1e-9)
carrier_pg     = cv[is_pg].sum(axis=0)     / np.maximum(cvc[is_pg].sum(axis=0), 1e-9)
EPS = 1.0/F
log2_pg_c = np.log2((carrier_pg + EPS) / (carrier_cactus + EPS))

# --- pre-mask PG het rate (mechanism 1 proxy) ---
pg_het = np.load(f'{SCRATCH}/pg_het_rate_chr1.npz', allow_pickle=True)
het_pos = pg_het['pos'].astype(np.int64)
het_rate = pg_het['het_rate'].astype(np.float32)
het_d = dict(zip(het_pos, het_rate))
pre_mask_het = np.array([het_d.get(p, np.nan) for p in chr1_snp_pos], dtype=np.float32)
print(f'pre-mask het rate matched: {np.isfinite(pre_mask_het).sum():,} / {len(pre_mask_het):,}')

# --- F_MISSING in v3qc_v3 ---
fmiss = 1.0 - cvc.sum(axis=0) / F

# Save aligned array for downstream stratification
df = pd.DataFrame({
    'pos': chr1_snp_pos,
    'cem': cem_af,
    'hf':  hf_af,
    'rec_lab':  recipe_lab,
    'rec_unif': recipe_unif,
    'fmiss': fmiss.astype(np.float32),
    'pre_mask_het': pre_mask_het,
    'log2_pg_c': log2_pg_c.astype(np.float32),
    'carrier_cactus': carrier_cactus.astype(np.float32),
    'carrier_pg':     carrier_pg.astype(np.float32),
    'single_alt': single_alt,
})
df.to_parquet(f'{SCRATCH}/triangulation_S1_chr1.parquet')
print(f'\nSaved {SCRATCH}/triangulation_S1_chr1.parquet  ({len(df):,} rows)')

# --- TOP-LEVEL: pairwise MAE/signed-bias against the lab-recipe anchor ---
def metrics(x, y, mask=None):
    m = np.isfinite(x) & np.isfinite(y) & (x>=0)&(x<=1)&(y>=0)&(y<=1)
    if mask is not None: m &= mask
    if m.sum() < 5: return None
    xv, yv = x[m], y[m]
    return dict(n=int(m.sum()),
                r2=float(np.corrcoef(xv,yv)[0,1]**2),
                mae=float(np.mean(np.abs(yv-xv))),
                signed=float(np.mean(yv-xv)),
                slope=float(np.polyfit(xv,yv,1)[0]),
                out=float(np.mean(np.abs(yv-xv)>0.1)))

print('\n=== TRIANGULATION TABLE (SEEDMIX_S1 Chr1 single-ALT SNPs) ===')
print('All comparisons on the same record set (single-ALT positions, all anchors finite, all estimators finite).')
all_finite = (np.isfinite(cem_af) & np.isfinite(hf_af)
              & np.isfinite(recipe_lab) & np.isfinite(recipe_unif)
              & single_alt)
print(f'n records: {all_finite.sum():,}')

def fmt(d):
    return f"n={d['n']:>7,}  R²={d['r2']:.4f}  MAE={d['mae']:.4f}  signed(y-x)={d['signed']:+.4f}  slope={d['slope']:+.3f}  out>0.1={d['out']*100:5.2f}%"

print('\n-- vs LAB RECIPE projected through cactus panel (h_lab @ cn_var / h_lab @ cn_var_called) --')
print(f'  cactus_em vs rec_lab : {fmt(metrics(df.rec_lab.values,  df.cem.values, all_finite))}')
print(f'  hapFIRE   vs rec_lab : {fmt(metrics(df.rec_lab.values,  df.hf.values,  all_finite))}')
print('\n-- vs UNIFORM-h projection (the "recipe FIXED" of PIPELINE_STATE) --')
print(f'  cactus_em vs rec_unif: {fmt(metrics(df.rec_unif.values, df.cem.values, all_finite))}')
print(f'  hapFIRE   vs rec_unif: {fmt(metrics(df.rec_unif.values, df.hf.values,  all_finite))}')
print('\n-- pairwise --')
print(f'  cactus_em vs hapFIRE : {fmt(metrics(df.hf.values,       df.cem.values, all_finite))}')
print(f'  rec_lab   vs hapFIRE : {fmt(metrics(df.hf.values,       df.rec_lab.values, all_finite))}')
print(f'  rec_unif  vs hapFIRE : {fmt(metrics(df.hf.values,       df.rec_unif.values, all_finite))}')
print(f'  rec_lab   vs rec_unif: {fmt(metrics(df.rec_unif.values, df.rec_lab.values, all_finite))}')
