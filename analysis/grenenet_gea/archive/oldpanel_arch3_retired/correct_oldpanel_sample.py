#!/usr/bin/env python
"""Correct one OLD-panel (10.33M-record) sample's global-mode TSVs to the new
segregating-only panel (8.49M) by row-subsetting with old2new_mask.

Rigorous because the segregating filter only changed V_pa (variant projection),
NOT K_pa (h estimation) -> global h identical -> kept-row AFs byte-identical to a
new-panel rerun. We subset LINE-WISE so kept rows are preserved byte-for-byte; the
dropped rows are the monomorphic founder variants the new panel excludes by design.

Originals are moved to ARCHIVE_DIR first. Idempotent: skips if already archived.

Usage: correct_oldpanel_sample.py SAMPLE ARCH3_DIR MASK_NPY ARCHIVE_DIR
"""
import sys, os, shutil
import numpy as np

SAMPLE, DIR, MASK_PATH, ARCH = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]
# fixed panel constants (old per-chrom record counts, Chr1..5); sum = 10,325,364
OLD_COUNTS = {"Chr1": 2617370, "Chr2": 1793129, "Chr3": 1936098, "Chr4": 1655579, "Chr5": 2323188}
NEW_COUNTS = {"Chr1": 2154423, "Chr2": 1463385, "Chr3": 1592677, "Chr4": 1368120, "Chr5": 1911041}
N_OLD = sum(OLD_COUNTS.values())

mask = np.load(MASK_PATH)
assert mask.shape == (N_OLD,), f"mask len {mask.shape} != {N_OLD}"
# per-chrom slices (genome-wide TSV concatenation order = Chr1..Chr5)
slices, off = {}, 0
for c in ["Chr1", "Chr2", "Chr3", "Chr4", "Chr5"]:
    slices[c] = mask[off:off + OLD_COUNTS[c]]; off += OLD_COUNTS[c]
assert off == N_OLD

os.makedirs(ARCH, exist_ok=True)

def subset_file(path, m, expect_out):
    """Move path -> ARCH, then stream-subset back to path keeping header + True rows."""
    base = os.path.basename(path)
    apath = os.path.join(ARCH, base)
    if os.path.exists(apath):
        print(f"  skip (already archived): {base}"); return None
    shutil.move(path, apath)
    n_in = n_out = 0
    with open(apath) as fi, open(path, "w") as fo:
        fo.write(fi.readline())  # header
        for i, line in enumerate(fi):
            if m[i]:
                fo.write(line); n_out += 1
            n_in += 1
    assert n_in == m.size, f"{base}: read {n_in} data rows, mask {m.size}"
    assert n_out == expect_out, f"{base}: wrote {n_out}, expected {expect_out}"
    return n_out

# per-chrom TSVs
for c in ["Chr1", "Chr2", "Chr3", "Chr4", "Chr5"]:
    p = os.path.join(DIR, f"{SAMPLE}_{c}.tsv")
    r = subset_file(p, slices[c], NEW_COUNTS[c])
    if r is not None:
        print(f"  {c}: {OLD_COUNTS[c]:,} -> {r:,}")
# genome-wide concatenated TSV
gw = os.path.join(DIR, f"{SAMPLE}.tsv")
r = subset_file(gw, mask, sum(NEW_COUNTS.values()))
if r is not None:
    print(f"  genome-wide: {N_OLD:,} -> {r:,}")
print(f"{SAMPLE}: corrected.")
