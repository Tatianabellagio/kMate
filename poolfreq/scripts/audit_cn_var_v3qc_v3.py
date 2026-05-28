"""Audit cn_var / cn_var_called consistency vs the haploid VCF F_MISSING.

Hypothesis (Mechanism 3): cn_var_called over-counts non-missing cells, so
projection denominator (h @ cn_var_called) is inflated and AF cactus_em
pulled toward zero.

Diagnostics:
  1. Per-record nnz of cn_var_called (= called founders).
  2. Compare to VCF F_MISSING (231 - F_MISSING*231) for a Chr1 sample.
  3. Identity check: nnz(cn_var) <= nnz(cn_var_called) per record.
  4. Sanity: are any GT codings being silently mapped to ALT or to "called"?
"""
from __future__ import annotations
import numpy as np
from pathlib import Path
from scipy.sparse import load_npz
import pysam
import sys, time, os

OUT_DIR = "/tmp"
BASE = str(Path(__file__).resolve().parents[2])
CV   = f"{BASE}/poolfreq/data/cn_var_231_v3qc_v3.cn_var.npz"
CVC  = f"{BASE}/poolfreq/data/cn_var_231_v3qc_v3.cn_var_called.npz"
META = f"{BASE}/poolfreq/data/cn_var_231_v3qc_v3.meta.npz"
VCF  = f"{BASE}/pangenie_genotyping/data/v3qc_v3/founders_231_v3qc_v3.haploid.vcf.gz"

t0 = time.time()

# Step 1: per-column nnz via bincount on raw CSR indices (no tocsc, no densify).
def col_nnz_from_csr(path, n_cols):
    print(f"[{time.time()-t0:6.0f}s] loading {path}", flush=True)
    m = load_npz(path)
    print(f"            shape={m.shape} nnz={m.nnz:,}", flush=True)
    # csr indices = column indices of nz entries (one per nz)
    out = np.bincount(m.indices, minlength=n_cols).astype(np.int32)
    del m
    return out

# Get n_records cheaply: load only the indptr/shape from cv
print(f"[{time.time()-t0:6.0f}s] peeking shape", flush=True)
_z = np.load(CV)
shape = tuple(_z['shape'])
del _z
F, N = shape
print(f"            F={F}  N_records={N:,}", flush=True)

ac = col_nnz_from_csr(CV, N)
an = col_nnz_from_csr(CVC, N)
print(f"[{time.time()-t0:6.0f}s] AC/AN computed", flush=True)

print("\n--- AC (cn_var per-record nnz) ---")
print(f"  min {ac.min()}  max {ac.max()}  median {int(np.median(ac))}")
print(f"  AC==0    : {(ac==0).sum():,}   (after AC=0 cleanup; expect 0)")
print(f"  AC==F    : {(ac==F).sum():,}   (fixed-ALT records)")

print("\n--- AN (cn_var_called per-record nnz = 2*(1-F_MISSING)*F for haploid F-counts) ---")
print(f"  min {an.min()}  max {an.max()}  median {int(np.median(an))}")
print(f"  AN==F    : {(an==F).sum():,}  ({100*(an==F).sum()/N:.2f}%)  -- no missing")
print(f"  AN==0    : {(an==0).sum():,}")
print(f"  AC > AN  : {(ac>an).sum():,}  (must be 0; ALT-carrier subset of called)")
for p in [1,5,25,50,75,95,99]:
    print(f"  AN p{p:>2} : {np.percentile(an, p):.0f}")

fm = 1.0 - an/F
print("\n--- F_MISSING (= 1 - AN/F) implied by cn_var_called ---")
bins = [(-1e-9, 0.0001), (0.0001, 0.05), (0.05, 0.1), (0.1, 0.2),
        (0.2, 0.5), (0.5, 0.8), (0.8, 1.01)]
for lo, hi in bins:
    if lo < 0:
        m = fm < hi
        label = f'F_MISSING == 0'
    else:
        m = (fm > lo) & (fm <= hi)
        label = f'{lo:.4f} < F <= {hi:.4f}'
    print(f"  {label:32s}: {m.sum():>10,}  ({100*m.sum()/N:.2f}%)")

# Step 2: compare to VCF on a Chr1 sample
print(f"\n[{time.time()-t0:6.0f}s] loading meta", flush=True)
meta = np.load(META, allow_pickle=True)
chrom_arr = np.asarray(meta['chrom']).astype(str)
pos_arr   = np.asarray(meta['pos']).astype(np.int64)
print(f"  meta chrom unique: {sorted(set(chrom_arr.tolist()))[:10]}")

# Step 3: re-read VCF, compute per-record AN_vcf for first N_CHECK records and
#         compare to AN-from-cn_var_called.
N_CHECK = 100000
print(f"\n[{time.time()-t0:6.0f}s] re-reading VCF Chr1 (first {N_CHECK} records) to check AN", flush=True)
vcf = pysam.VariantFile(VCF)
samples = list(vcf.header.samples)
assert len(samples) == F, f'sample-count mismatch: {len(samples)} vs F={F}'
vcf_an = np.zeros(N_CHECK, dtype=np.int32)
vcf_ac = np.zeros(N_CHECK, dtype=np.int32)
vcf_chrom = []
vcf_pos = []
vcf_ref = []
vcf_alt = []
gt_value_seen = {}    # distinct GT raw values across all sample cells
i = 0
for rec in vcf.fetch('Chr1'):
    if i >= N_CHECK:
        break
    if not rec.alts:
        continue
    vcf_chrom.append(rec.chrom)
    vcf_pos.append(rec.pos)
    vcf_ref.append(rec.ref)
    vcf_alt.append(rec.alts[0])
    called = 0
    alt_carrier = 0
    for sname in samples:
        gt = rec.samples[sname]['GT']
        if gt is None or len(gt) == 0:
            key = 'empty'
        else:
            key = repr(gt)
        gt_value_seen[key] = gt_value_seen.get(key, 0) + 1
        if gt is None or len(gt) == 0:
            continue
        is_missing = all(a is None for a in gt)
        if is_missing:
            continue
        called += 1
        if any(a is not None and a > 0 for a in gt):
            alt_carrier += 1
    vcf_an[i] = called
    vcf_ac[i] = alt_carrier
    i += 1

n = i
vcf_an = vcf_an[:n]
vcf_ac = vcf_ac[:n]
print(f"  read {n} Chr1 records from VCF")

# Match by (chrom, pos, ref, alt) — multi-allelic records need this
chr1_mask = chrom_arr == 'Chr1'
chr1_idx_in_cnvar = np.where(chr1_mask)[0]
# meta retains record order; assume first N_CHECK Chr1 in meta == first N_CHECK Chr1 in VCF
matched = min(n, len(chr1_idx_in_cnvar))
ac_cnvar = ac[chr1_idx_in_cnvar[:matched]]
an_cnvar = an[chr1_idx_in_cnvar[:matched]]
vcf_an_m = vcf_an[:matched]
vcf_ac_m = vcf_ac[:matched]

an_diff = an_cnvar - vcf_an_m
ac_diff = ac_cnvar - vcf_ac_m
print(f"\n--- AN  (cn_var_called)  vs  AN_vcf  on Chr1 first {matched} records ---")
print(f"  AN     mean: cn_var={an_cnvar.mean():.3f}  vcf={vcf_an_m.mean():.3f}")
print(f"  AN diff: nonzero records {(an_diff!=0).sum()} / {matched}")
print(f"  AN diff: min {an_diff.min()}  max {an_diff.max()}  mean {an_diff.mean():.4f}")
print(f"\n--- AC  (cn_var)         vs  AC_vcf ---")
print(f"  AC     mean: cn_var={ac_cnvar.mean():.3f}  vcf={vcf_ac_m.mean():.3f}")
print(f"  AC diff: nonzero records {(ac_diff!=0).sum()} / {matched}")
print(f"  AC diff: min {ac_diff.min()}  max {ac_diff.max()}  mean {ac_diff.mean():.4f}")

# Show distinct GT codings encountered
print("\n--- distinct GT codings encountered (raw repr(gt) from pysam) ---")
total = sum(gt_value_seen.values())
for k, v in sorted(gt_value_seen.items(), key=lambda x: -x[1])[:20]:
    print(f"  {v:>12,} ({100*v/total:5.2f}%)  {k}")

# F_MISSING distribution from VCF (haploid → 1 - AN_vcf/F)
fm_vcf = 1 - vcf_an_m/F
print("\n--- F_MISSING from VCF (first 100k Chr1) ---")
for lo, hi in bins:
    if lo < 0:
        m_ = fm_vcf < hi
        label = 'F_MISSING == 0'
    else:
        m_ = (fm_vcf > lo) & (fm_vcf <= hi)
        label = f'{lo:.4f} < F <= {hi:.4f}'
    print(f"  {label:32s}: {m_.sum():>10,}  ({100*m_.sum()/matched:.2f}%)")

np.savez(f"{OUT_DIR}/cn_var_v3qc_v3_audit.npz",
         ac=ac, an=an,
         vcf_an_chr1_first=vcf_an_m, vcf_ac_chr1_first=vcf_ac_m,
         cnvar_an_chr1_first=an_cnvar, cnvar_ac_chr1_first=ac_cnvar)
print(f"\n[{time.time()-t0:6.0f}s] saved {OUT_DIR}/cn_var_v3qc_v3_audit.npz")
