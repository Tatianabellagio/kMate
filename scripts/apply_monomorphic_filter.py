#!/usr/bin/env python
"""Regenerate the production panel to SEGREGATING-ONLY (drop monomorphic records).

Monomorphic = AC=0 (ALT carried by no founder) or AC=AN (ALT carried by all called).
This matches the patched jobA4 Step-3 filter:  bcftools view -e 'AC=0 || AC=AN'.

For each chrom (overwrites in place, via temp+atomic-replace):
  1. final VCF        : bcftools view -e 'INFO/AC=0 || INFO/AC=INFO/AN'  (re-index)
  2. var_pa/var_called: column-subset existing CSR by mask (car>0)&(car<cal)
  3. meta             : subset chrom/pos/ref/alt/ref_len/alt_len by the same mask

Subsetting the existing matrices is provably identical to rebuilding build_var_pa on
the filtered VCF (V_pa is 1:1 with the VCF; we remove whole columns whose carrier/called
sets are unchanged). Verified by: per-chrom count match (matrix mask == filtered VCF
record count) + first/last position spot-check + a final invariant re-audit (separate).
"""
import os, gc, subprocess
import numpy as np
import scipy.sparse as sp

PROJ = "/global/scratch/users/tbellg/kmate"
BCF = "/global/home/users/tbellg/miniforge3/envs/kmate/bin/bcftools"
CHROMS = [1, 2, 3, 4, 5]
EXPECT = {1: 2154423, 2: 1463385, 3: 1592677, 4: 1368120, 5: 1911041}  # from invariant audit

grand_in = grand_out = 0
for c in CHROMS:
    cl = f"chr{c}"
    vcf = f"{PROJ}/panel/arch3/{cl}/merged_231_{cl}_final.vcf.gz"
    pre = f"{PROJ}/panel/arch3/{cl}/var_pa_231_arch3_{cl}"
    vp_f, vc_f, vm_f = pre + ".var_pa.npz", pre + ".var_called.npz", pre + ".meta.npz"

    print(f"\n=== Chr{c} ===", flush=True)
    vp = sp.load_npz(vp_f).tocsc()       # 231 x N ; CSC for fast column subsetting
    vc = sp.load_npz(vc_f).tocsc()
    N = vp.shape[1]
    car = np.asarray(vp.sum(axis=0)).ravel()
    cal = np.asarray(vc.sum(axis=0)).ravel()
    mask = (car > 0) & (car < cal)       # segregating: AC>0 and AC<AN
    n_keep = int(mask.sum())
    print(f"  records {N:,} -> segregating {n_keep:,} "
          f"(drop AC=0:{int((car==0).sum()):,}  AC=AN:{int((car==cal).sum()):,})", flush=True)
    assert n_keep == EXPECT[c], f"mask count {n_keep} != expected {EXPECT[c]}"
    grand_in += N; grand_out += n_keep

    # --- 1. filter VCF in place (atomic) ---
    tmp_vcf = vcf + ".seg.tmp.vcf.gz"
    subprocess.run([BCF, "view", "-e", "INFO/AC=0 || INFO/AC=INFO/AN",
                    vcf, "--threads", "4", "-Oz", "-o", tmp_vcf], check=True)
    subprocess.run([BCF, "index", "-t", "-f", tmp_vcf], check=True)
    n_vcf = int(subprocess.run([BCF, "index", "-n", tmp_vcf],
                               capture_output=True, text=True).stdout.strip())
    assert n_vcf == n_keep, f"filtered VCF {n_vcf} != matrix mask {n_keep}"
    print(f"  filtered VCF records: {n_vcf:,}  (== matrix mask) OK", flush=True)

    # --- 2. subset matrices ---
    vp2 = vp[:, mask].tocsr()
    vc2 = vc[:, mask].tocsr()
    sp.save_npz(vp_f + ".tmp.npz", vp2)   # save_npz appends .npz
    sp.save_npz(vc_f + ".tmp.npz", vc2)

    # --- 3. subset meta ---
    m = np.load(vm_f, allow_pickle=True)
    meta_out = {}
    for k in m.files:
        a = m[k]
        meta_out[k] = a[mask] if (a.ndim == 1 and a.shape[0] == N) else a
    np.savez(vm_f + ".tmp", **meta_out)

    # --- spot-check meta order vs filtered VCF (first & last 3 positions) ---
    vpos = subprocess.run([BCF, "query", "-f", "%POS\n", tmp_vcf],
                          capture_output=True, text=True).stdout.split()
    vfirst = [int(x) for x in vpos[:3]]; vlast = [int(x) for x in vpos[-3:]]
    mpos = meta_out["pos"]
    assert list(mpos[:3]) == vfirst and list(mpos[-3:]) == vlast, \
        f"order mismatch! VCF {vfirst}/{vlast} vs meta {list(mpos[:3])}/{list(mpos[-3:])}"
    print(f"  order spot-check OK (first {vfirst}, last {vlast})", flush=True)

    # --- atomic replace ---
    os.replace(tmp_vcf, vcf); os.replace(tmp_vcf + ".tbi", vcf + ".tbi")
    os.replace(vp_f + ".tmp.npz", vp_f); os.replace(vc_f + ".tmp.npz", vc_f)
    os.replace(vm_f + ".tmp.npz", vm_f)
    print(f"  overwrote VCF + var_pa + var_called + meta in place", flush=True)
    del vp, vc, vp2, vc2, car, cal, mask, m; gc.collect()

print(f"\nDONE. Genome-wide: {grand_in:,} -> {grand_out:,} segregating "
      f"({grand_in-grand_out:,} monomorphic dropped, {100*(grand_in-grand_out)/grand_in:.2f}%)")
