#!/usr/bin/env python3
"""
F2-honest patch: at cn_var SNP records, override cactus founders' GT from 0->1
where the RAW deconstruct VCF (biallelic-normed) says they actually carry the
SNP base on their haplotype path.

This addresses Mechanism 2 (bubble flattening) at its source: cactus founders
that traverse an SV branch but whose path still passes through the SNP-ALT base
get correctly counted as SNP-ALT carriers, instead of being defaulted to REF
by the vcfbub-l-0 multi-allelic decomposition.

Notes / approximations:
  - We re-project AF using the RECIPE-space (uniform h=1/231) formula, which
    gives the exact F2-honest answer when h is uniform. For non-uniform h (real
    EM), we add a multiplicative-correction approximation.
  - We only override 0->1 (add missing carriers). We do NOT override 1->0 — the
    current cn_var GT=1 founders are trusted.
"""
import numpy as np
import pandas as pd
import scipy.sparse as sp
from pathlib import Path
from scipy.stats import pearsonr

ROOT = Path('/carnegie/nobackup/scratch/tbellagio/hapfire_sv')
GT_TSV = ROOT / 'scratch/f2_honest/raw_chr1_biallelic_GT.tsv'

# ---------- 1) Build Assembly_ID -> Accession_ID mapping ----------
print('[1] Build founder ID mapping')
xl = pd.read_excel(ROOT / 'data/ASSEMBLIES_Best_version_of_dataset.xlsx',
                   usecols=['Assembly_ID','Accession_ID'])
xl['Assembly_ID'] = xl['Assembly_ID'].astype(str)
xl['Accession_ID'] = xl['Accession_ID'].astype(str)
asm_to_acc = dict(zip(xl['Assembly_ID'], xl['Accession_ID']))
print(f'    loaded {len(asm_to_acc):,} Assembly_ID → Accession_ID rows')

# ---------- 2) Load cn_var meta + raw VCF sample header ----------
print('[2] Load cn_var meta + raw VCF sample list')
m = np.load(ROOT / 'poolfreq/data/cn_var_231_v3qc_v3.meta.npz', allow_pickle=True)
founders = list(m['founders'])
chrom_arr = m['chrom']; pos_arr = m['pos']
ref_arr   = m['ref'];   alt_arr = m['alt']
ref_len_arr = m['ref_len']; alt_len_arr = m['alt_len']
print(f'    cn_var founders: {len(founders)} (e.g., {founders[:5]})')

# Raw VCF column header: get the Assembly_ID list (82 samples)
raw_samples = None
with open(ROOT / 'scratch/f2_honest/raw_chr1.vcf.gz', 'rb') as f:
    import gzip
    with gzip.open(f, 'rt') as gz:
        for line in gz:
            if line.startswith('#CHROM'):
                raw_samples = line.rstrip().split('\t')[9:]
                break
print(f'    raw VCF samples: {len(raw_samples)} (e.g., {raw_samples[:5]})')

# Map raw col idx -> cn_var founder row idx (or -1 if not present)
founder_idx = {f: i for i, f in enumerate(founders)}
raw_col_to_cnvar_row = []
n_mapped = 0
n_excluded = 0
for asm_id in raw_samples:
    acc_id = asm_to_acc.get(asm_id)
    if acc_id is None:
        raw_col_to_cnvar_row.append(-1); n_excluded += 1; continue
    cnvar_idx = founder_idx.get(acc_id, -1)
    raw_col_to_cnvar_row.append(cnvar_idx)
    if cnvar_idx >= 0:
        n_mapped += 1
    else:
        n_excluded += 1
raw_col_to_cnvar_row = np.array(raw_col_to_cnvar_row, dtype=np.int32)
print(f'    raw samples mapped to cn_var founders : {n_mapped}')
print(f'    raw samples NOT in cn_var (excluded)  : {n_excluded}  '
      f'(expected: 101003=5772 dup and 100852=Ped-0 — excluded from v3qc)')

# ---------- 3) Stream the GT TSV to build dict (chrom,pos,ref,alt) -> 82-int8 mask ----------
print('[3] Stream GT TSV, build carrier dict (Chr1 biallelic SNPs only)')
import time
t0 = time.time()
carrier_dict = {}     # (pos,ref,alt) -> np.int8[82]  (1 = ALT, 0 = REF, -1 = missing)
n_total = 0
n_snp   = 0
n_dup   = 0
with open(GT_TSV) as f:
    for line in f:
        parts = line.rstrip('\n').split('\t')
        chrom = parts[0]; pos = int(parts[1])
        ref = parts[2]; alt = parts[3]
        n_total += 1
        # SNP only for this F2-honest test
        if len(ref) != 1 or len(alt) != 1:
            continue
        if ref == '*' or alt == '*':
            continue
        gts = parts[4:]
        # Carrier vec: '1' -> 1, '0' -> 0, '.' -> -1
        # Use ord comparison for speed
        carrier = np.zeros(len(gts), dtype=np.int8)
        for i, g in enumerate(gts):
            if g == '1':
                carrier[i] = 1
            elif g == '.':
                carrier[i] = -1
            # '0' (or anything else) stays 0
        key = (pos, ref, alt)
        if key in carrier_dict:
            # multiple LV records map to same biallelic — OR the carriers (any LV=carrier wins)
            prev = carrier_dict[key]
            # Treat -1 as "no info" → only OR positive carriers; keep prev if both are 0
            merged = np.where(prev == 1, 1, np.where(carrier == 1, 1, np.where((prev == -1) & (carrier == -1), -1, 0)))
            carrier_dict[key] = merged.astype(np.int8)
            n_dup += 1
        else:
            carrier_dict[key] = carrier
        n_snp += 1
        if n_total % 500000 == 0:
            print(f'    streamed {n_total:,} records ({time.time()-t0:.1f}s), SNPs kept: {n_snp:,}')
print(f'    total records streamed : {n_total:,}')
print(f'    biallelic SNPs kept    : {n_snp:,}')
print(f'    duplicate SNP keys (multiple LVs collapsed): {n_dup:,}')
print(f'    unique (pos,ref,alt) keys                  : {len(carrier_dict):,}')
print(f'    elapsed: {time.time()-t0:.1f}s')

# ---------- 4) Load existing cn_var; identify Chr1 SNP records at mixed bubbles ----------
print('[4] Load cn_var; identify Chr1 SNP records')
cv = sp.load_npz(ROOT / 'poolfreq/data/cn_var_231_v3qc_v3.cn_var.npz')
cv_csc = cv.tocsc()
print(f'    cn_var: {cv.shape}, {cv.nnz:,} nnz')

chr1_mask = (chrom_arr == 'Chr1')
snp_mask  = chr1_mask & (ref_len_arr == 1) & (alt_len_arr == 1)
snp_cols  = np.where(snp_mask)[0]
print(f'    Chr1 biallelic SNP cn_var cols: {len(snp_cols):,}')

# Load F1 test joined table to know which records are mixed bubbles (n_sv_carriers_at_p > 0)
F1 = pd.read_csv(ROOT / 'scratch/hetmask_diag/f1_test_joined.tsv', sep='\t')
print(f'    F1 joined table: {len(F1):,} (Chr1 biallelic SNPs joined with hapFIRE+mech2)')

# ---------- 5) Build override map: for each cn_var SNP col, list of cactus founders to 0->1 ----------
print('[5] Build override list (cactus founders that should be 0->1 at each SNP)')
# raw_col_to_cnvar_row[col] = cn_var row index (or -1)
valid_raw = raw_col_to_cnvar_row >= 0     # 80 valid mapped cols
valid_raw_idx = np.where(valid_raw)[0]    # indices into raw 82-vec
valid_cnvar_rows = raw_col_to_cnvar_row[valid_raw]  # cn_var rows for these

n_overrides_per_col = np.zeros(len(snp_cols), dtype=np.int32)
overrides_applied = []   # list of (cn_var_col, cn_var_row) tuples
not_found_in_raw = 0
found_in_raw    = 0

for ix, c in enumerate(snp_cols):
    pos = int(pos_arr[c]); ref = str(ref_arr[c]); alt = str(alt_arr[c])
    key = (pos, ref, alt)
    raw_vec = carrier_dict.get(key)
    if raw_vec is None:
        not_found_in_raw += 1
        continue
    found_in_raw += 1
    # current cn_var carrier mask for valid cactus founders only
    col_data = np.asarray(cv_csc[:, c].todense()).ravel().astype(np.int8)
    current_at_cactus = col_data[valid_cnvar_rows]    # 0/1 for cactus founders only
    raw_at_cactus     = raw_vec[valid_raw_idx]        # 0/1/-1 for cactus founders only
    # Override mask: where raw says ALT and current says REF
    override = (raw_at_cactus == 1) & (current_at_cactus == 0)
    if override.any():
        n_overrides_per_col[ix] = override.sum()
        rows_to_flip = valid_cnvar_rows[override]
        for r in rows_to_flip:
            overrides_applied.append((c, r))
print(f'    cn_var SNP records found in raw biallelic dict : {found_in_raw:,}')
print(f'    cn_var SNP records NOT in raw dict            : {not_found_in_raw:,}')
print(f'    cn_var SNP cols with at least 1 override      : {(n_overrides_per_col>0).sum():,}')
print(f'    total founder x record overrides              : {len(overrides_applied):,}')
print(f'    n_overrides_per_col stats:')
arr = n_overrides_per_col[n_overrides_per_col > 0]
print(f'      mean   : {arr.mean():.2f}   median: {int(np.median(arr))}   max: {arr.max()}')

# ---------- 6) Compute F2-honest AFs (recipe-space and cactus_em-approx) ----------
print('[6] Re-project AFs after F2-honest override')
# Need AC_post and AN_post per cn_var Chr1 SNP record (already on F1 joined)
# And we need to know n_overrides at the cn_var record level (not the F1-joined level)
# Build a per-cnvar-col override count
ovr_per_col = pd.DataFrame({
    'col': snp_cols,
    'n_override': n_overrides_per_col,
})
# Build meta DF for Chr1 SNP records (cn_var-col indexed)
snp_meta = pd.DataFrame({
    'col':  snp_cols,
    'chrom': chrom_arr[snp_cols],
    'pos':   pos_arr[snp_cols],
    'ref':   ref_arr[snp_cols],
    'alt':   alt_arr[snp_cols],
    'n_override': n_overrides_per_col,
})
# Read SEEDMIX_S1 cactus_em alt_freq for these rows
cem_tsv = pd.read_csv(ROOT / 'scratch/v3qc_v3_mixedloose_chr1/SEEDMIX_S1.tsv', sep='\t')
# TSV row order == cn_var col order
snp_meta['alt_freq_orig'] = cem_tsv.iloc[snp_cols]['alt_freq'].values

# Merge with F1 joined table for AC/AN/F_MISSING/hapfire/is_mixed_bubble/n_sv_carriers_at_p
M = F1.merge(snp_meta, on=['chrom','pos','ref','alt'], how='inner')
print(f'    merged F1 joined with override table: {len(M):,} records')
# Reconciliation: alt_freq must match
assert np.allclose(M['alt_freq_orig'], M['alt_freq'], atol=1e-5), 'alt_freq mismatch'

# RECIPE-space exact F2-honest:
M['recipe_F2honest'] = (M['AC_post'] + M['n_override']) / M['AN_post']
M['recipe_F2honest'] = M['recipe_F2honest'].clip(0,1)

# cactus_em-space F2-honest under uniform-h approximation
# AF_new = AF_orig + n_override / AN_post   (each override adds 1/231 of mass; D = AN_post/231)
M['cactus_em_F2honest'] = (M['alt_freq'] + M['n_override'] / M['AN_post']).clip(0,1)

# Residuals
def res_stats(df, est_col):
    r = df['hapfire_af'] - df[est_col]
    return (r.mean(), r.median(), r.abs().mean(),
            (r.abs()>0.05).mean()*100, (r.abs()>0.10).mean()*100)

print('\n=== Per-cell agreement (hapFIRE − estimator) ===\n')
cells = {
    (False,0): '[a] low_FMISS  pure_SNP',
    (False,1): '[b] low_FMISS  mixed_bub',
    (True,0):  '[c] high_FMISS pure_SNP',
    (True,1):  '[d] high_FMISS mixed_bub',
}
rows = []
for (lo,mb), label in cells.items():
    sub = M[(M['high_fmiss']==lo) & (M['is_mixed_bubble']==mb)]
    if len(sub)==0: continue
    r_orig = sub['hapfire_af'] - sub['alt_freq']
    r_f2h  = sub['hapfire_af'] - sub['cactus_em_F2honest']
    r_rec  = sub['hapfire_af'] - sub['recipe_F2honest']
    rows.append({
        'cell': label, 'n': len(sub),
        'mean_orig':    r_orig.mean(),
        'mean_F2honest':r_f2h.mean(),
        'mean_recipeF2':r_rec.mean(),
        'MAE_orig':     r_orig.abs().mean(),
        'MAE_F2honest': r_f2h.abs().mean(),
        'MAE_recipeF2': r_rec.abs().mean(),
        'out05_orig':   (r_orig.abs()>0.05).mean()*100,
        'out05_F2honest':(r_f2h.abs()>0.05).mean()*100,
    })
tab = pd.DataFrame(rows)
print(tab.round(4).to_string(index=False))

print('\n=== Genome-wide (all 516k joined Chr1 SNPs) ===\n')
for est, name in [('alt_freq', 'cactus_em (orig)'),
                  ('cactus_em_F2honest','cactus_em + F2-honest (raw VCF)'),
                  ('recipe_af', 'recipe (orig)'),
                  ('recipe_F2honest', 'recipe + F2-honest')]:
    r = M['hapfire_af'] - M[est]
    print(f'  {name:36s}  mean={r.mean():+.5f}  MAE={r.abs().mean():.5f}  '
          f'%|d|>0.05={(r.abs()>0.05).mean()*100:.2f}%  %|d|>0.10={(r.abs()>0.10).mean()*100:.2f}%')

# Save patched table
M.to_csv(ROOT / 'scratch/f2_honest/f2_honest_test.tsv', sep='\t', index=False)
print(f'\nWrote {ROOT / "scratch/f2_honest/f2_honest_test.tsv"}')
