#!/usr/bin/env python3
"""
Forensic trace of Chr1:13843898 (C→T in cn_var, AF 0.08 cactus vs 0.99 hapFIRE).

For each step in the pipeline, show what's there and how it's interpreted:

  Step 1 — RAW cactus deconstruct VCF at this position
  Step 2 — vcfbub -l 0 output at this position
  Step 3 — bcftools norm -m -any decomposition at this position
  Step 4 — merged 231-panel cn_var at this position
  Step 5 — base-at-position truth for each founder (from raw multi-allelic GT)
  Step 6 — what xwu's GrENE-Net SNP catalog says at this position
  Step 7 — what hapFIRE projects
"""
import subprocess
from pathlib import Path

ROOT = Path('/carnegie/nobackup/scratch/tbellagio/hapfire_sv')
BCF = '/home/tbellagio/miniforge3/envs/sequencing_pipeline/bin/bcftools'

POS = 13843898
print(f'=' * 80)
print(f'  FORENSIC TRACE  Chr1:{POS}  (REF=C, ALT=T in cn_var)')
print(f'=' * 80)

# ---------- Step 1: RAW cactus deconstruct (LV=2 nested) ----------
print(f'\n[Step 1] RAW cactus deconstruct VCF — pre-vcfbub')
print(f'  File: pang_1001gplus_82acc.raw.vcf')
print(f'  Search for any record overlapping pos {POS}:')
raw = '/home/tbellagio/scratch/pang/pang_1001gplus/pang/output/pang_1001gplus_82acc.raw.vcf'
res = subprocess.run(['awk', '-F\t', f'!/^#/ && $1=="Chr1" && $2 >= {POS-2} && $2 <= {POS+2}'
                      ' {print "  pos="$2," REF="$4," ALT="$5," QUAL="$6," INFO[0:200]="substr($8,1,200)}',
                      raw], capture_output=True, text=True)
print(res.stdout)

# ---------- Step 2: vcfbub-l-0 (drops LV>=1) ----------
print(f'\n[Step 2] vcfbub --max-level 0 output (production cactus VCF)')
print(f'  File: pang_1001gplus_82acc.vcf.gz')
vcfbub = '/home/tbellagio/scratch/pang/pang_1001gplus/pang/output/pang_1001gplus_82acc.vcf.gz'
res = subprocess.run([BCF, 'view', '-H', '-r', f'Chr1:{POS-2}-{POS+2}', vcfbub],
                     capture_output=True, text=True)
for line in res.stdout.strip().split('\n')[:5]:
    if line:
        f = line.split('\t')
        print(f'  pos={f[1]} REF={f[3]} ALT={f[4]} INFO[0:150]={f[7][:150]}')

# ---------- Step 3: bcftools norm -m -any decomposition ----------
print(f'\n[Step 3] After bcftools norm -m -any')
print(f'  Decomposes the multi-allelic into K biallelic records:')
print()
# Get the multi-allelic record from raw and decompose it manually
mult = subprocess.run([BCF, 'view', '-H', '-r', f'Chr1:{POS-2}-{POS+2}', vcfbub],
                      capture_output=True, text=True).stdout.strip().split('\n')
for line in mult:
    if not line: continue
    f = line.split('\t')
    pos = int(f[1]); ref = f[3]; alts = f[4].split(',')
    info = f[7]
    if pos != POS: continue
    # AC field
    ac = next((s for s in info.split(';') if s.startswith('AC=')), 'AC=?')
    an = next((s for s in info.split(';') if s.startswith('AN=')), 'AN=?')
    print(f'  Multi-allelic: pos={pos} REF={ref} ALT={alts}  {ac}  {an}')
    print(f'  Per-ALT decomposition (canonical biallelic form after norm):')
    for i, alt in enumerate(alts, 1):
        # Strip shared prefix
        ix = 0
        while ix < len(ref)-1 and ix < len(alt)-1 and ref[ix]==alt[ix]:
            ix += 1
        r2, a2 = ref[ix:], alt[ix:]
        # Strip shared suffix
        sx = 0
        while sx < len(r2)-1 and sx < len(a2)-1 and r2[-1-sx]==a2[-1-sx]:
            sx += 1
        if sx > 0:
            r2, a2 = r2[:-sx], a2[:-sx]
        is_snp = len(r2)==1 and len(a2)==1
        print(f'    ALT_{i}: original  REF={ref} ALT={alt}   →   biallelic REF={r2} ALT={a2}  '
              f'{"(SNP)" if is_snp else ""}')

# ---------- Step 4: cn_var records at this position ----------
print(f'\n[Step 4] cn_var records on the merged 231-panel haploid VCF')
print(f'  File: founders_231_v3qc_v3.haploid.vcf.gz')
merged = ROOT / 'pangenie_genotyping/data/v3qc_v3/founders_231_v3qc_v3.haploid.vcf.gz'
res = subprocess.run([BCF, 'query', '-r', f'Chr1:{POS}-{POS}',
                      '-f', '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AC\t%INFO/AN\t%INFO/F_MISSING\n',
                      str(merged)], capture_output=True, text=True)
print(f'  cn_var records at pos {POS}:')
for line in res.stdout.strip().split('\n'):
    if line:
        f = line.split('\t')
        print(f'    REF={f[2]} ALT={f[3]}  AC={f[4]}  AN={f[5]}  F_MISSING={f[6]}')

# ---------- Step 5: Base-at-position truth from raw multi-allelic GTs ----------
print(f'\n[Step 5] Base at pos {POS} per founder — derived from raw multi-allelic GTs')
print(f'  For each founder, their GT (0/1/2/3) maps to an ALT path, which has a specific base at pos {POS}:')
# Get GTs from the vcfbub-l-0 file (which has the multi-allelic record)
res = subprocess.run([BCF, 'query', '-r', f'Chr1:{POS}-{POS}',
                      '-f', '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n', vcfbub],
                     capture_output=True, text=True)
samples = subprocess.run([BCF, 'view', '-h', vcfbub], capture_output=True, text=True)
sample_line = [l for l in samples.stdout.split('\n') if l.startswith('#CHROM')][0]
sample_names = sample_line.rstrip().split('\t')[9:]
n_samples = len(sample_names)
print(f'  cactus samples in vcfbub-l-0 VCF: {n_samples}')

for line in res.stdout.strip().split('\n'):
    if not line: continue
    f = line.split('\t')
    pos = int(f[1]); ref = f[2]; alts = f[3].split(',')
    if pos != POS: continue
    gts = f[4:]
    # Map each GT index to base at pos POS within (ref/alt sequence)
    # base_at_pos_per_allele[0] = ref[0]; base_at_pos_per_allele[i] = alts[i-1][0] (since alts/ref start at POS)
    bases = [ref[0]] + [a[0] for a in alts]
    print(f'  Multi-allelic: REF={ref}({bases[0]}) ' + ', '.join(
        [f'ALT_{i+1}={alts[i]}({bases[i+1]})' for i in range(len(alts))]))
    # Tally per-founder bases at pos POS
    base_counts = {}
    n_missing = 0
    for g in gts:
        g = g.strip()
        if g in ('.', './.', '.|.'):
            n_missing += 1; continue
        try:
            idx = int(g.split('/')[0].split('|')[0])
        except:
            n_missing += 1; continue
        if idx < len(bases):
            b = bases[idx]
            base_counts[b] = base_counts.get(b, 0) + 1
        else:
            n_missing += 1
    print(f'  Base-at-pos tally across {n_samples} cactus founders:')
    for b, c in sorted(base_counts.items(), key=lambda x: -x[1]):
        print(f'    base "{b}": {c} founders  ({100*c/n_samples:.1f}%)')
    print(f'    missing:        {n_missing} founders')
    total_called = sum(base_counts.values())
    if 'T' in base_counts:
        af_T_at_pos = base_counts['T'] / total_called
        print(f'\n  ⇒ TRUE base-at-pos AF for T: {base_counts["T"]} / {total_called} = {af_T_at_pos:.3f}')
        print(f'    (this is what hapFIRE/xwu should be measuring)')
    # Count how many founders are coded GT==2 (TT path, the "pure SNP" decomposition)
    n_gt2 = sum(1 for g in gts if g.strip() == '2')
    print(f'\n  vs. cn_var path-frequency for (C, T) biallelic:')
    print(f'    only founders with GT==2 (TT path) get cn_var=1 → {n_gt2} cactus founders')
    print(f'    Founders with GT==1 (TC path) or GT==3 (TG path) are coded GT=0 at (C,T) biallelic,')
    print(f'    even though their actual base at pos {POS} is "T".')

# ---------- Step 6: xwu's GrENE-Net SNP catalog at this position ----------
print(f'\n[Step 6] xwu GrENE-Net SNP catalog at pos {POS}')
# TODO: generate chr1 split from /global/scratch/projects/fc_moilab/projects/grenenet-phase1/vcf/greneNet_final_v1.1.recode.vcf via "bcftools view -r 1"
XWU = '/global/scratch/users/tbellg/hapfire_sv/data/greneNet_final_v1.1_chr1.recode.vcf'
res = subprocess.run(['awk', '-F\t', f'!/^#/ && $1=="1" && $2=={POS}',
                      XWU], capture_output=True, text=True)
for line in res.stdout.strip().split('\n'):
    if not line: continue
    f = line.split('\t')
    pos = int(f[1]); ref = f[3]; alt = f[4]
    gts = f[9:]
    ac = an = 0
    for g in gts:
        if g in ('.', './.', '.|.', '.|', '|.'): continue
        an += 1
        # phased "a|b" or unphased "a/b"
        try:
            a,b = (int(x) for x in g.replace('|','/').split('/')[:2])
            if a>0 or b>0: ac += 1
        except:
            pass
    print(f'  Record: pos={pos} REF={ref} ALT={alt}  AC={ac}/{an}  →  AF={ac/max(an,1):.3f}')

# ---------- Step 7: hapFIRE projection ----------
HAPF = '/global/scratch/projects/fc_moilab/projects/grenenet-phase1/frequency/hapFIRE_frequencies/seed_mix/s1_snp_frequency.txt'
res = subprocess.run(['awk', '-F\t', f'$1=="1" && $2=={POS}', HAPF],
                     capture_output=True, text=True)
print(f'\n[Step 7] hapFIRE projection (s1_snp_frequency.txt):')
print(f'  {res.stdout.strip()}')

print(f'\n' + '=' * 80)
print(f'  CONCLUSION')
print(f'=' * 80)
print('''  The merge join (chrom, pos, ref, alt) is CORRECT — both cn_var and hapFIRE
  refer to TAIR10 pos 13843898 REF=C ALT=T. They are NOT comparing different SNPs.

  What differs is the *semantic* of AF:
    cn_var encodes "founder took the C→T pure-SNP allele branch in the cactus snarl"
    hapFIRE/xwu encodes  "founder genome has base T at this TAIR10 coordinate"

  At a multi-allelic snarl where multiple ALTs share the SNP base at pos p:
    - All ALT-carriers have T at pos p in their actual sequence
    - But only one ALT (the one that, after bcftools-norm trimming, becomes "C→T")
      gets path-coded as cn_var=1 at the SNP biallelic
    - Founders on the other ALTs are coded as path=0 at the SNP biallelic

  This is the bias source. The fix is to compute base-at-position frequency from
  the multi-allelic record directly (which is what Step 5 above shows we can do).''')
