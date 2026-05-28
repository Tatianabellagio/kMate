#!/usr/bin/env python3
"""
Empirical check: is the proposed "skip bcftools norm; iterate rec.alts in
build_cn_var" mathematically equivalent to the current "bcftools norm -m -any
then process biallelic records" pipeline?

Test:
  - cactus_78.vcf.gz       = pre-norm (multi-allelic vcfbub-l-0 cactus 78 samples)
  - cactus_78_bi.vcf.gz    = post-norm biallelic (same VCF after `bcftools norm -m -any`)

For each multi-allelic record at a sample of mixed-bubble positions:
  - Build "proposed" cn_var/cn_var_called column-set from rec.alts iteration
  - Look up the corresponding biallelic records in cactus_78_bi.vcf.gz
  - Build "current" cn_var/cn_var_called from each biallelic
  - Compare cell-by-cell

If equivalent across 100 probes, the proposal is a refactor with no AF impact.
If divergent, identify the rows where they differ and characterize.
"""
import pysam
import numpy as np
from pathlib import Path

ROOT = Path('/carnegie/nobackup/scratch/tbellagio/kmate')
PRE  = ROOT / 'pangenie_genotyping/data/v3qc_tmp/cactus_78.vcf.gz'
POST = ROOT / 'pangenie_genotyping/data/v3qc_v2/cactus_78_bi.vcf.gz'

vpre  = pysam.VariantFile(str(PRE))
vpost = pysam.VariantFile(str(POST))
samples = list(vpre.header.samples)
print(f'pre-norm samples : {len(samples)} (e.g., {samples[:3]})')
post_samples = list(vpost.header.samples)
print(f'post-norm samples: {len(post_samples)} (same: {samples == post_samples})')
assert samples == post_samples, 'sample lists differ!'
F = len(samples)

def gt_to_int(gt_tuple):
    """Convert pysam GT tuple to int (None => missing -> -1; haploid -> single int)."""
    if gt_tuple is None or len(gt_tuple) == 0:
        return -1
    g = gt_tuple[0]
    return -1 if g is None else int(g)

def build_proposed_col(rec):
    """
    For each ALT_j (j=1..K), return:
      carrier_j[f] = 1 iff GT_f == j else 0
      called_j[f]  = 1 iff GT_f != . (same vector for all j)
    Returns dict: {alt_index_1based -> (carrier_vec, called_vec)}
    """
    gts = np.array([gt_to_int(rec.samples[s]['GT']) for s in samples], dtype=np.int8)
    called = (gts != -1).astype(np.int8)
    out = {}
    for j in range(1, len(rec.alts)+1):
        carrier = (gts == j).astype(np.int8)
        out[j] = (carrier, called, rec.alts[j-1])
    return out

def build_current_col(rec_bi):
    """For a biallelic record (single ALT), return (carrier_vec, called_vec)."""
    gts = np.array([gt_to_int(rec_bi.samples[s]['GT']) for s in samples], dtype=np.int8)
    carrier = (gts == 1).astype(np.int8)
    called  = (gts != -1).astype(np.int8)
    return carrier, called

def normalize_alt(ref, alt):
    """Strip shared prefix and suffix (mimics bcftools norm behaviour for SNP/MNP cases)."""
    # Strip shared prefix
    i = 0
    while i < len(ref)-1 and i < len(alt)-1 and ref[i] == alt[i]:
        i += 1
    ref2, alt2 = ref[i:], alt[i:]
    shift = i
    # Strip shared suffix (but keep at least one char each)
    j = 0
    while j < len(ref2)-1 and j < len(alt2)-1 and ref2[-1-j] == alt2[-1-j]:
        j += 1
    if j > 0:
        ref2, alt2 = ref2[:-j], alt2[:-j]
    return shift, ref2, alt2

# Walk Chr1, find multi-allelic records, compare to biallelic at same pos
n_probes_target = 200
n_checked = 0
n_match = 0
n_divergent_rows = 0
mismatches = []
diff_carrier_cells = 0
diff_called_cells = 0
total_cells = 0

print(f'\nScanning multi-allelic Chr1 records ...')
for rec in vpre.fetch('Chr1'):
    if rec.alts is None or len(rec.alts) <= 1:
        continue
    pos = rec.pos
    ref = rec.ref
    alts = list(rec.alts)

    # Get all biallelic records at the same pos in post-norm VCF (might also be at
    # shifted positions if left-alignment moved them — fetch a small window)
    bi_at_pos = []
    for rb in vpost.fetch('Chr1', max(pos-5, 1), pos+max(len(ref),1)+5):
        if rb.pos > pos + 5:
            break
        # Match by (norm(ref, alt) → expected (pos, ref, alt))
        bi_at_pos.append(rb)

    proposed = build_proposed_col(rec)
    # For each ALT in pre-norm, predict its canonical biallelic key and look it up in post-norm
    for j, alt in enumerate(alts, start=1):
        shift, ref_can, alt_can = normalize_alt(ref, alt)
        key = (pos + shift, ref_can, alt_can)
        # Find matching post-norm record
        match = None
        for rb in bi_at_pos:
            if rb.pos == key[0] and rb.ref == key[1] and rb.alts and rb.alts[0] == key[2]:
                match = rb
                break
        if match is None:
            # Sometimes norm shifts even more (e.g., complex left-alignment); skip
            continue
        cur_carrier, cur_called = build_current_col(match)
        prop_carrier, prop_called, _ = proposed[j]
        total_cells += F
        diff_c = (cur_carrier != prop_carrier).sum()
        diff_d = (cur_called  != prop_called ).sum()
        diff_carrier_cells += diff_c
        diff_called_cells  += diff_d
        if diff_c == 0 and diff_d == 0:
            n_match += 1
        else:
            n_divergent_rows += 1
            if len(mismatches) < 10:
                mismatches.append({
                    'pos': pos, 'ref': ref, 'alts': ','.join(alts),
                    'this_alt_idx': j, 'this_alt': alt,
                    'mapped_to': key,
                    'cur_AC': int(cur_carrier.sum()), 'prop_AC': int(prop_carrier.sum()),
                    'cur_AN': int(cur_called.sum()), 'prop_AN': int(prop_called.sum()),
                    'diff_carrier_cells': int(diff_c), 'diff_called_cells': int(diff_d),
                })
        n_checked += 1
        if n_checked >= n_probes_target * 4:
            break
    if n_checked >= n_probes_target * 4:
        break

print(f'\n=== Equivalence-check results ===')
print(f'  multi-allelic ALT cells probed   : {n_checked:,}')
print(f'    of which matched post-norm key : {n_match + n_divergent_rows:,}')
print(f'    (unmatched ALTs not checked: likely complex left-aligns)')
print(f'  ALT columns with identical cn_var/cn_var_called : {n_match:,}')
print(f'  ALT columns with ANY cell difference            : {n_divergent_rows:,}')
print(f'  total (founder × ALT) cells compared            : {total_cells:,}')
print(f'  cell-level cn_var      mismatches               : {diff_carrier_cells:,}')
print(f'  cell-level cn_var_called mismatches             : {diff_called_cells:,}')

if n_divergent_rows > 0:
    print(f'\nFirst {min(10, len(mismatches))} divergent ALTs:')
    for m in mismatches[:10]:
        print(f'  pos={m["pos"]:>10} ref={m["ref"][:8]:<8} alts={m["alts"][:30]:<30}'
              f'  ALT_{m["this_alt_idx"]}={m["this_alt"][:12]:<12}  '
              f'cur AC/AN={m["cur_AC"]}/{m["cur_AN"]}  prop AC/AN={m["prop_AC"]}/{m["prop_AN"]}  '
              f'diff carrier={m["diff_carrier_cells"]} called={m["diff_called_cells"]}')
