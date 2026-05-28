"""Count what MAC<2 would drop on the now-biallelic-correctly-filtered pangenie_153_qc_v2 VCF.

Reports:
  - Records by AC class (0, 1, 2-10, >10, AN-AC<2)
  - Records by variant class (SNP / small_indel / SV / MNP)
  - "Hom-backed" drops — answers Tatiana's concern: of records that MAC would drop, how many
    have AC_Hom >= 1 (i.e., some real hom-alt evidence)?
"""
import subprocess, sys
from pathlib import Path
import numpy as np

BCF = '/global/home/users/tbellg/miniforge3/envs/sequencing_pipeline/bin/bcftools'
VCF = str(Path(__file__).resolve().parents[2] / 'panel/pangenie_genotyping/data/v3qc_v2/pangenie_153_qc_v2.vcf.gz')

CHROMS = sys.argv[1:] if len(sys.argv) > 1 else ['Chr1']

for chrom in CHROMS:
    print(f'\n========== {chrom} ==========')
    # Re-fill-tags to get AC_Hom (separate from AC_Het) — fill-tags only added AC_Het, we also need AC_Hom
    cmd = f"""{BCF} +fill-tags {VCF} -r {chrom} -- -t AC,AN,AC_Het,AC_Hom 2>/dev/null | \\
{BCF} query -f '%REF\\t%ALT\\t%INFO/AC\\t%INFO/AN\\t%INFO/AC_Het\\t%INFO/AC_Hom\\n'"""
    res = subprocess.run(cmd, shell=True, capture_output=True, text=True, check=True)
    lines = res.stdout.strip().split('\n')
    n = len(lines)
    print(f'records on {chrom}: {n:,}')

    ref_len = np.zeros(n, dtype=int)
    alt_len = np.zeros(n, dtype=int)
    ac = np.zeros(n, dtype=int)
    an = np.zeros(n, dtype=int)
    ac_het = np.zeros(n, dtype=int)
    ac_hom = np.zeros(n, dtype=int)
    for i, line in enumerate(lines):
        parts = line.split('\t')
        try:
            ref_len[i] = len(parts[0])
            alt_len[i] = len(parts[1])
            ac[i] = int(parts[2])
            an[i] = int(parts[3])
            ac_het[i] = int(parts[4])
            ac_hom[i] = int(parts[5])
        except Exception:
            pass

    # Variant class
    abs_d = np.abs(alt_len - ref_len)
    cls = np.where((ref_len == 1) & (alt_len == 1), 'SNP',
        np.where(abs_d == 0, 'MNP',
        np.where(abs_d < 50, 'small_indel', 'SV>=50bp')))

    # MAC drops
    mac_lower = ac < 2          # AC=0 or 1
    mac_upper = (an - ac) < 2   # mono-ALT or AC=AN-1
    mac_drop = mac_lower | mac_upper

    print(f'\n  MAC<2 would drop: {mac_drop.sum():,} ({100*mac_drop.sum()/n:.2f}%)')
    print(f'    AC<2 (lower-end): {mac_lower.sum():,}')
    print(f'    AN-AC<2 (upper-end): {mac_upper.sum():,}')
    print(f'    overlap: {(mac_lower & mac_upper).sum():,}')

    # By class
    print(f'\n  Drop by variant class:')
    for c in ['SNP', 'small_indel', 'SV>=50bp', 'MNP']:
        m = cls == c
        d = (mac_drop & m).sum()
        total = m.sum()
        print(f'    {c:>14}: {d:>10,} / {total:>10,} ({100*d/total if total else 0:.2f}%)')

    # ★ User's concern: do MAC drops have hom-alt backing?
    print(f'\n  ★ MAC drops by hom-alt backing (AC_Hom):')
    mac_drop_lower = mac_drop & mac_lower
    n_drop = mac_drop_lower.sum()
    ac_hom_dropped = ac_hom[mac_drop_lower]
    print(f'    AC<2 drops: {n_drop:,}')
    print(f'      AC_Hom = 0 (pure-het or mono-REF, no hom backing): {(ac_hom_dropped == 0).sum():,} ({100*(ac_hom_dropped==0).sum()/max(1,n_drop):.2f}%)')
    print(f'      AC_Hom >= 1 (HAS hom backing, user-concern case): {(ac_hom_dropped >= 1).sum():,} ({100*(ac_hom_dropped>=1).sum()/max(1,n_drop):.2f}%)')

    # AC distribution sample
    print(f'\n  AC=0,1,AN sample (first 10 each):')
    ac0 = (ac == 0).sum()
    ac1 = (ac == 1).sum()
    ac2 = (ac == 2).sum()
    acaN = (ac == an).sum()
    print(f'    AC=0: {ac0:,}')
    print(f'    AC=1: {ac1:,}')
    print(f'    AC=2: {ac2:,}')
    print(f'    AC=AN: {acaN:,}')
