#!/bin/bash
#SBATCH --job-name=chr1_proj_h
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --output=logs/A6_proj_%j.out
#SBATCH --error=logs/A6_proj_%j.err
mkdir -p logs
set -euo pipefail

# Phase 2 Job A6: re-project existing SEEDMIX_S1 h vector through NEW var_pa → per-record AF.
# Why we can do this: kmer_pa is unchanged → existing h is still correct.
# Only var_pa (the projection matrix) changed in Arch 3.

cd /global/scratch/users/tbellg/kmate/panel/arch3/chr1
PY=/global/home/users/tbellg/miniforge3/envs/hapfm/bin/python

H_PATH=/global/scratch/users/tbellg/kmate/scratch/v3qc_v3_mixedloose_chr1/SEEDMIX_S1.h_per_chrom.npz
CN_VAR=var_pa_231_arch3_chr1.var_pa.npz
VAR_CALLED=var_pa_231_arch3_chr1.var_called.npz
CN_VAR_META=var_pa_231_arch3_chr1.meta.npz

[ -s "$H_PATH" ] || { echo "ERROR: missing h vector"; exit 1; }
[ -s "$CN_VAR" ] || { echo "ERROR: missing var_pa (A5 not done)"; exit 1; }

OUT=SEEDMIX_S1_arch3_chr1.tsv

$PY -u <<EOF
import numpy as np
from scipy.sparse import load_npz
import time

print('[load] h vector')
h_data = np.load('$H_PATH', allow_pickle=True)
print(f'  keys: {list(h_data.files)}')
assert 'Chr1' in h_data.files, f"expected 'Chr1' key in h_per_chrom.npz, got {list(h_data.files)}"
assert 'founders' in h_data.files, f"expected 'founders' key in h_per_chrom.npz"
h_raw = h_data['Chr1']
h_founders = h_data['founders']
print(f'  h shape: {h_raw.shape}, sum: {h_raw.sum():.4f}')
print(f'  h founders (first 3): {h_founders[:3]}')

print('[load] var_pa + var_called')
t0 = time.time()
var_pa = load_npz('$CN_VAR').tocsr()
var_called = load_npz('$VAR_CALLED').tocsr()
print(f'  var_pa: {var_pa.shape}, {var_pa.nnz:,} nnz')
print(f'  var_called: {var_called.shape}, {var_called.nnz:,} nnz')
print(f'  loaded in {time.time()-t0:.1f}s')

print('[load] meta')
meta = np.load('$CN_VAR_META', allow_pickle=True)
chrom = meta['chrom']; pos = meta['pos']; ref = meta['ref']; alt = meta['alt']
ref_len = meta['ref_len']; alt_len = meta['alt_len']
cn_founders = meta['founders']
print(f'  meta records: {len(chrom):,}')
print(f'  var_pa founders (first 3): {cn_founders[:3]}')

# CRITICAL: align h to var_pa founder ordering, since h was built against a different
# panel build (v3qc_v3 var_pa) whose sample ordering may differ.
print('[align] verify h founders == var_pa founders, reindex if needed')
h_founders_s = [str(x) for x in h_founders]
cn_founders_s = [str(x) for x in cn_founders]
if h_founders_s == cn_founders_s:
    print('  identical ordering — no reindex needed')
    h = h_raw.astype(float)
else:
    print('  ordering DIFFERS — reindexing h to var_pa founder order')
    h_lookup = {str(name): float(h_raw[i]) for i, name in enumerate(h_founders)}
    missing = [n for n in cn_founders_s if n not in h_lookup]
    extra = [n for n in h_founders_s if n not in set(cn_founders_s)]
    if missing:
        raise SystemExit(f'ERROR: {len(missing)} var_pa founders missing from h: {missing[:5]}...')
    if extra:
        print(f'  WARN: {len(extra)} h founders not in var_pa (dropped): {extra[:5]}')
    h = np.array([h_lookup[n] for n in cn_founders_s], dtype=float)
    print(f'  reindexed h shape: {h.shape}  sum: {h.sum():.6f} (should be ~1.0)')
    assert abs(h.sum() - h_raw.sum()) < 1e-9, 'reindex changed h sum — extra founders dropped'
F = h.shape[0]
assert F == var_pa.shape[0], f'h has {F} entries but var_pa has {var_pa.shape[0]} founder rows'

print('[project] AF = h @ var_pa / h @ var_called')
t0 = time.time()
numer = h @ var_pa          # shape (n_records,)
denom = h @ var_called
af = np.where(denom > 0, numer / denom, np.nan)
print(f'  projection in {time.time()-t0:.1f}s')
print(f'  AF range: [{np.nanmin(af):.4f}, {np.nanmax(af):.4f}], NaNs: {np.isnan(af).sum()}')

print('[write] TSV')
with open('$OUT', 'w') as f:
    f.write('chrom\tpos\tref_len\talt_len\talt_freq\n')
    for i in range(len(chrom)):
        f.write(f'{chrom[i]}\t{pos[i]}\t{ref_len[i]}\t{alt_len[i]}\t{af[i]:.5f}\n')
print(f'  wrote $OUT')

# Spot-checks at the 3 canonical positions
print()
print('=== Spot-check AF at 3 canonical positions ===')
for tp, tr, ta in [(5870018,'T','A'), (10421645,'T','C'),
                    (13843898,'C','T')]:
    matches = np.where((pos == tp) & (np.array([str(r) for r in ref]) == tr) & (np.array([str(a) for a in alt]) == ta))[0]
    if len(matches) == 0:
        print(f'  Chr1:{tp} {tr}>{ta}: NOT IN var_pa meta')
        continue
    for idx in matches:
        print(f'  Chr1:{tp} {tr}>{ta}: AF={af[idx]:.4f}  (numer={numer[idx]:.4f} denom={denom[idx]:.4f})')
EOF

echo
ls -lh $OUT
echo "[$(date)] DONE A6"
