#!/usr/bin/env python3
"""
COORDINATE-AWARE base-at-position test at Chr1:10421645.

Find ALL cn_var records that span the target coordinate (not just records
starting at it). For each, compute the base in the ALT sequence at the target
coordinate. Count carriers across all records where the ALT-base matches.

If this lifts AF from current 0.02 to ~0.99 (matching xwu/hapFIRE), we've
identified that the position-shifting from `bcftools norm` is the real issue
with my earlier base-at-pos fix.
"""
import numpy as np
import scipy.sparse as sp
from pathlib import Path

ROOT = Path('/carnegie/nobackup/scratch/tbellagio/kmate')

print('Loading cn_var and meta...')
cv = sp.load_npz(ROOT / 'data/cn_var_231_v3qc_v3.cn_var.npz').tocsc()
m = np.load(ROOT / 'data/cn_var_231_v3qc_v3.meta.npz', allow_pickle=True)
chrom_arr = m['chrom']; pos_arr = m['pos']
ref_arr = m['ref']; alt_arr = m['alt']
ref_len_arr = m['ref_len']; alt_len_arr = m['alt_len']

TARGET_POS = 10421645
TARGET_BASE = 'C'  # the ALT base we want (T → C SNP)
TARGET_CHROM = 'Chr1'

# Find ALL cn_var records on Chr1 that SPAN TARGET_POS
# i.e., pos <= TARGET_POS < pos + ref_len
chr1_mask = (chrom_arr == TARGET_CHROM)
chr1_cols = np.where(chr1_mask)[0]

# Filter: pos <= TARGET_POS AND pos + ref_len > TARGET_POS
span_mask = (pos_arr[chr1_cols] <= TARGET_POS) & \
            (pos_arr[chr1_cols] + ref_len_arr[chr1_cols] > TARGET_POS)
span_cols = chr1_cols[span_mask]
print(f'\ncn_var records on Chr1 that span pos {TARGET_POS}: {len(span_cols)}')

# For each, compute the ALT base at the target coordinate
# offset = TARGET_POS - record_pos
# alt_base_at_target = ALT[offset]  (if offset < len(ALT))
# Special handling: if ALT is shorter than offset, the founder's sequence at
# TARGET_POS is undefined (deletion). Treat as 'missing'.

print(f'\nLooking for records where ALT_carrier has base "{TARGET_BASE}" at coord {TARGET_POS}:')
matching_cols = []
matching_records = []
for c in span_cols:
    p = int(pos_arr[c]); ref = str(ref_arr[c]); alt = str(alt_arr[c])
    offset = TARGET_POS - p
    # ALT-carrier's sequence at coord TARGET_POS = ALT[offset], assuming alt-len > offset
    if offset >= len(alt):
        alt_base = 'DEL'  # carrier has deletion at this coord
    else:
        alt_base = alt[offset].upper()
    if alt_base == TARGET_BASE:
        ac = int(np.asarray(cv[:, c].todense()).ravel().sum())
        an = (cv[:, c].nnz / cv[:, c].shape[0]) if False else None  # placeholder
        matching_cols.append(c)
        matching_records.append({
            'pos': p, 'ref_len': len(ref), 'alt_len': len(alt),
            'ref_abbr': ref if len(ref)<=10 else ref[:5]+'...'+ref[-3:],
            'alt_abbr': alt if len(alt)<=10 else alt[:5]+'...'+alt[-3:],
            'offset_in_alt': offset,
            'alt_base_at_target': alt_base,
            'AC': ac,
        })

# Sort by position
matching_records.sort(key=lambda r: r['pos'])
print(f'\nMatching records (ALT has "{TARGET_BASE}" at coord {TARGET_POS}):')
print(f'  {"pos":>10} {"ref_len":>7} {"alt_len":>7} {"REF":>15} {"ALT":>15} {"offset":>6} {"base":>4} {"AC":>4}')
total_ac = 0
for r in matching_records:
    print(f'  {r["pos"]:>10} {r["ref_len"]:>7} {r["alt_len"]:>7} {r["ref_abbr"]:>15} '
          f'{r["alt_abbr"]:>15} {r["offset_in_alt"]:>6} {r["alt_base_at_target"]:>4} {r["AC"]:>4}')
    total_ac += r['AC']

print(f'\n  Total AC across all matching records (raw sum): {total_ac}')

# Compute UNIQUE founder mass — since founders are mutually exclusive across the same
# multi-allelic but DIFFERENT multi-allelics can overlap, we need to OR the columns
print('\nORing the carrier vectors across matching cn_var columns:')
combined = np.zeros(231, dtype=np.int8)
for c in matching_cols:
    col_data = np.asarray(cv[:, c].todense()).ravel().astype(np.int8)
    combined = np.maximum(combined, col_data)
total_carriers = int(combined.sum())
print(f'  Total UNIQUE founders carrying "{TARGET_BASE}" at coord {TARGET_POS}: {total_carriers}')

# Compare to current cn_var AC and to xwu/hapFIRE
print('\n=== Comparison ===')
print(f'  Current cn_var (T,C) record only:              AC=2   AN=182   AF=0.011')
print(f'  My naive base-at-pos (same-pos siblings only): AC=4   AN=182   AF=0.022 (only +2 extras)')
print(f'  Coordinate-aware (all spanning records):       AC={total_carriers}  AN=182   '
      f'AF={total_carriers/182:.3f}')
print(f'  xwu GrENE-Net (truth):                         AC=229 AN=231   AF=0.991')
print(f'  hapFIRE projection:                            AF=0.989')
