#!/usr/bin/env python3
"""
Independent validation of Chr1:5870018 via minimap2 alignment of a 1kb TAIR10
window to founder genome FASTAs. For each founder, find the projected coord
and read the actual base.

Cactus pangenome (post-atomize): 70/78 founders carry T->A here
xwu/hapFIRE: 5/231 founders carry T->A  (~2/78 cactus)

Goal: which is right? minimap2 alignment to founder genomes is independent of both.
"""
import subprocess, re, os
from pathlib import Path

ROOT = Path('/carnegie/nobackup/scratch/tbellagio/hapfire_sv')
FA_DIR = Path('/home/tbellagio/scratch/pang/pang_1001gplus/20260209_Exposito-Alonso/chr_only')
TAIR = FA_DIR / 'TAIR10.chr.fa'
MINIMAP = '/home/tbellagio/miniforge3/envs/sequencing_pipeline/bin/minimap2'
SAMTOOLS = '/home/tbellagio/miniforge3/envs/sequencing_pipeline/bin/samtools'
OUT_DIR = ROOT / 'scratch/atomize_test/validation'
OUT_DIR.mkdir(parents=True, exist_ok=True)

POS = 5870018
WIN = 500
QUERY_FA = OUT_DIR / 'tair10_query.fa'

# Extract TAIR10 query window
subprocess.run([SAMTOOLS, 'faidx', str(TAIR),
                f'Chr1:{POS-WIN}-{POS+WIN}'],
                stdout=open(QUERY_FA, 'w'), check=True)
query_seq = ''.join(open(QUERY_FA).read().split('\n')[1:])
target_q_pos = WIN   # 0-indexed offset of POS in query
print(f'TAIR10 reference base at coord {POS}: {query_seq[target_q_pos]}')
print(f'(query 1kb window pos {POS-WIN}-{POS+WIN}; SNP at offset {target_q_pos})')

# Sample 12 founders (Assembly_IDs)
import random
random.seed(42)
all_fastas = sorted(FA_DIR.glob('1[0-9][0-9][0-9][0-9][0-9].chr.fa'))
sample = random.sample(all_fastas, min(12, len(all_fastas)))

def align_and_project(query_fa, founder_fa, target_q):
    """minimap2 align query to founder; project target_q onto founder coord."""
    proc = subprocess.run(
        [MINIMAP, '-a', '-x', 'asm5', '--secondary=no',
         str(founder_fa), str(query_fa)],
        capture_output=True, text=True
    )
    if proc.returncode != 0:
        return None
    for line in proc.stdout.split('\n'):
        if line.startswith('@') or not line.strip():
            continue
        fields = line.split('\t')
        flag = int(fields[1])
        if flag & 4:   # unmapped
            continue
        rname = fields[2]
        start = int(fields[3])    # 1-based
        cigar = fields[5]
        is_rev = bool(flag & 16)
        # Walk CIGAR to project target_q to founder coord
        q_pos = 0
        r_pos = start - 1   # 0-based founder coord
        for length, op in re.findall(r'(\d+)([MIDNSH=X])', cigar):
            length = int(length)
            if op in ('S', 'H'):
                continue
            if op in ('M', '=', 'X'):
                if q_pos + length > target_q:
                    delta = target_q - q_pos
                    return rname, r_pos + delta + 1, is_rev
                q_pos += length
                r_pos += length
            elif op == 'I':
                if q_pos + length > target_q:
                    return rname, None, is_rev   # in inserted region
                q_pos += length
            elif op in ('D', 'N'):
                r_pos += length
        return None
    return None

def get_base(fa_path, contig, coord, is_rev):
    res = subprocess.run([SAMTOOLS, 'faidx', str(fa_path), f'{contig}:{coord}-{coord}'],
                         capture_output=True, text=True)
    seq = res.stdout.split('\n')[1].strip()
    if not seq:
        return None
    base = seq[0].upper()
    if is_rev:
        rc = {'A':'T','T':'A','G':'C','C':'G','N':'N'}
        base = rc.get(base, 'N')
    return base

print(f'\nValidating Chr1:{POS} across 12 founders via minimap2 + samtools:')
print(f'{"FOUNDER":<10} {"contig":<14} {"founder_coord":<14} {"strand":<8} {"BASE":<6}')
print('-' * 60)
results = {'A': 0, 'T': 0, 'C': 0, 'G': 0, 'N': 0, 'INS': 0, 'NOHIT': 0}
for fa_path in sample:
    fid = fa_path.stem.replace('.chr', '')
    proj = align_and_project(QUERY_FA, fa_path, target_q_pos)
    if proj is None:
        print(f'{fid:<10} {"":<14} {"NO ALIGNMENT":<14}')
        results['NOHIT'] += 1
        continue
    contig, coord, is_rev = proj
    if coord is None:
        print(f'{fid:<10} {contig:<14} {"INSERTION":<14}')
        results['INS'] += 1
        continue
    base = get_base(fa_path, contig, coord, is_rev)
    strand = 'reverse' if is_rev else 'forward'
    print(f'{fid:<10} {contig:<14} {coord:<14} {strand:<8} {base}')
    results[base if base in 'ACGT' else 'N'] += 1

print()
print('=== Tally (per minimap2 founder genome alignment) ===')
for k, v in results.items():
    print(f'  base/state "{k}": {v}')
n_called = sum(results[b] for b in 'ACGT')
print()
print('=== Interpretation ===')
if n_called > 0:
    print(f'  Of {n_called} founders with successful base call:')
    print(f'    fraction T (TAIR10 REF):                 {results["T"]/n_called:.3f}')
    print(f'    fraction A (the alleged "T->A SNP" ALT): {results["A"]/n_called:.3f}')
    print(f'    fraction C:                              {results["C"]/n_called:.3f}')
    print(f'    fraction G:                              {results["G"]/n_called:.3f}')
print()
print('=== Pipeline claims (78 cactus founders) ===')
print('  Cactus pangenome (post-atomize):      70 carry T->A  → AF=0.90')
print('  Cactus pangenome (no atomize - current pipeline): 2/78 → AF=0.026')
print('  xwu GrENE-Net short-read SNP catalog:  5/231 → AF=0.022')
print('  hapFIRE projection:                                AF=0.022')
