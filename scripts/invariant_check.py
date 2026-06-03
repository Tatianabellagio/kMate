#!/usr/bin/env python
"""Invariant-site audit of the production panel (VCF + V_pa + K_pa).

Invariant = non-segregating ALT allele:
  - AC=0          ALT carried by no founder  (monomorphic reference)
  - AC=AN (>0)    ALT carried by all called founders (monomorphic alt)
  - AN=0          no founder called at all (fully missing)

Matrix level (definitive — the actual genotype matrix):
  car = per-record carriers (var_pa col-sum); cal = per-record called (var_called col-sum)
VCF level (independent cross-check, from INFO/AC, INFO/AN).
K_pa is already invariant-filtered by filt2inv (verified separately by AC range 2..230).
"""
import os, gc, subprocess
import numpy as np

PROJ = "/global/scratch/users/tbellg/kmate"
BCF = "/global/home/users/tbellg/miniforge3/envs/kmate/bin/bcftools"
NF = 231
CHROMS = [1, 2, 3, 4, 5]


def colcounts(npz_path, ncols):
    z = np.load(npz_path, allow_pickle=False)
    c = np.bincount(z["indices"], minlength=ncols).astype(np.int64)
    del z; gc.collect()
    return c


print("chrom |        N |   matrix: AC=0   AC=AN(>0)   AN=0 |    VCF: AC=0  AC=AN(>0)  AN=0")
print("-" * 92)
tot = dict(N=0, m_ac0=0, m_fix=0, m_an0=0, v_ac0=0, v_fix=0, v_an0=0)

for c in CHROMS:
    cl = f"chr{c}"
    vp = f"{PROJ}/panel/arch3/{cl}/var_pa_231_arch3_{cl}.var_pa.npz"
    vc = f"{PROJ}/panel/arch3/{cl}/var_pa_231_arch3_{cl}.var_called.npz"
    vm = f"{PROJ}/panel/arch3/{cl}/var_pa_231_arch3_{cl}.meta.npz"
    vcf = f"{PROJ}/panel/arch3/{cl}/merged_231_{cl}_final.vcf.gz"

    N = int(np.load(vm, allow_pickle=True)["pos"].shape[0])
    car = colcounts(vp, N)
    cal = colcounts(vc, N)
    m_ac0 = int((car == 0).sum())
    m_fix = int(((car == cal) & (cal > 0)).sum())
    m_an0 = int((cal == 0).sum())
    del car, cal; gc.collect()

    # VCF cross-check: stream AC/AN
    p = subprocess.run([BCF, "query", "-f", "%INFO/AC\t%INFO/AN\n", vcf],
                       capture_output=True, text=True)
    v_ac0 = v_fix = v_an0 = 0
    for line in p.stdout.splitlines():
        ac_s, an_s = line.split("\t")
        try:
            ac = int(ac_s.split(",")[0]); an = int(an_s)
        except ValueError:
            continue
        if an == 0:
            v_an0 += 1
        elif ac == 0:
            v_ac0 += 1
        elif ac == an:
            v_fix += 1

    print(f"Chr{c}   | {N:>8,} |        {m_ac0:>7,} {m_fix:>9,} {m_an0:>6,} |     "
          f"{v_ac0:>7,} {v_fix:>9,} {v_an0:>6,}")
    tot["N"] += N
    tot["m_ac0"] += m_ac0; tot["m_fix"] += m_fix; tot["m_an0"] += m_an0
    tot["v_ac0"] += v_ac0; tot["v_fix"] += v_fix; tot["v_an0"] += v_an0

print("-" * 92)
print(f"TOTAL  | {tot['N']:>8,} |        {tot['m_ac0']:>7,} {tot['m_fix']:>9,} {tot['m_an0']:>6,} |     "
      f"{tot['v_ac0']:>7,} {tot['v_fix']:>9,} {tot['v_an0']:>6,}")
m_inv = tot["m_ac0"] + tot["m_fix"] + tot["m_an0"]
v_inv = tot["v_ac0"] + tot["v_fix"] + tot["v_an0"]
print(f"\nV_pa invariant records (any of the 3): {m_inv:,}  ({100*m_inv/tot['N']:.3f}% of {tot['N']:,})")
print(f"VCF  invariant records (any of the 3): {v_inv:,}  ({100*v_inv/tot['N']:.3f}% of {tot['N']:,})")
