#!/usr/bin/env python3
"""
Tie-breaker: at cell [d] positions where cactus and xwu disagree on per-founder
calls, look up the actual base in each founder's long-read assembly FASTA.

Method:
  1. Pick 5 high-leverage probes where cactus and xwu disagree strongly.
  2. For each probe, take a 200bp window from TAIR10 centered on the SNP.
  3. For each founder FASTA, find the unique flanking-context match
     (use 80bp left flank as the search key; should be unique enough at most positions).
  4. Read the base at position +80 (the SNP position relative to the match).
  5. Tabulate: how many founders show REF, SNP-ALT, or neither?

Independent of both cactus pangenome graph and xwu short-read VCF.
"""
import subprocess
import pandas as pd
import numpy as np
from pathlib import Path

ROOT  = Path('/carnegie/nobackup/scratch/tbellagio/hapfire_sv')
FA_DIR= Path('/home/tbellagio/scratch/pang/pang_1001gplus/20260209_Exposito-Alonso/chr_only')
TAIR  = FA_DIR / 'TAIR10.chr.fa'
SAMTOOLS = '/home/tbellagio/miniforge3/envs/sequencing_pipeline/bin/samtools'

V = pd.read_csv(ROOT / 'scratch/f2_honest/celld_3source_validation.tsv', sep='\t')
# Pick top-5 most-disagreed probes (cactus AC very low, xwu AC very high)
V['disagreement'] = V['AC_xwu_any']/V['AN_xwu'] - V['AC_cactus']/V['AN_cactus']
top5 = V.nlargest(5, 'disagreement')[['pos','ref','alt','AC_cactus','AN_cactus',
                                     'AC_xwu_any','AN_xwu','cactus_em','hapfire_af']]
print('5 top-disagreement cell [d] probes:')
print(top5.to_string(index=False))

# Load Assembly_ID -> Accession_ID
xl = pd.read_excel(ROOT / 'data/ASSEMBLIES_Best_version_of_dataset.xlsx',
                   usecols=['Assembly_ID','Accession_ID'])
xl['Assembly_ID'] = xl['Assembly_ID'].astype(str)
xl['Accession_ID'] = xl['Accession_ID'].astype(str)
asm_to_acc = dict(zip(xl['Assembly_ID'], xl['Accession_ID']))
acc_to_asm = {}
for asm,acc in asm_to_acc.items():
    acc_to_asm.setdefault(acc, []).append(asm)

# All available founder FASTAs
founder_fastas = sorted(FA_DIR.glob('1*.chr.fa'))
print(f'\nFounder FASTAs available: {len(founder_fastas)}')

def get_tair_flank(pos, ref, alt, flank=80):
    """Get [flank_left][SNP][flank_right] from TAIR10 chr1."""
    start = pos - flank
    end   = pos + flank
    out = subprocess.run([SAMTOOLS, 'faidx', str(TAIR), f'Chr1:{start}-{end}'],
                        capture_output=True, text=True)
    seq = ''.join(out.stdout.split('\n')[1:])
    assert seq[flank] == ref, f'TAIR10 base mismatch: expected {ref} got {seq[flank]} at pos {pos}'
    return seq[:flank], seq[flank+1:flank+1+flank]   # left, right (each 80bp)

def search_founder_fasta(fa_path, left_flank, right_flank):
    """
    Find the left_flank in the founder FASTA; return the base immediately after.
    Returns (base, found_count, location_or_None).
    """
    # samtools faidx if indexed; else simple grep approach
    # We grep for the left_flank as text, with case-insensitive matching
    # Plus also check reverse-complement for safety
    rc = {'A':'T','T':'A','G':'C','C':'G','N':'N'}
    left_rc = ''.join(rc.get(c.upper(),'N') for c in reversed(left_flank))
    right_rc = ''.join(rc.get(c.upper(),'N') for c in reversed(right_flank))

    # Read fasta as continuous Chr1 sequence
    seq_lines = []
    in_chr1 = False
    with open(fa_path) as f:
        for line in f:
            if line.startswith('>'):
                hdr = line[1:].strip().split()[0]
                in_chr1 = hdr.lower() in ('chr1','1','chromosome1') or hdr.lower().startswith('chr1')
                if seq_lines and not in_chr1:
                    break
                continue
            if in_chr1:
                seq_lines.append(line.strip().upper())
    if not seq_lines:
        return ('?', 0, None)
    seq = ''.join(seq_lines)

    # Forward search
    fwd_matches = []
    start = 0
    while True:
        idx = seq.find(left_flank.upper(), start)
        if idx < 0: break
        # Check that there's a base after the left flank
        if idx + len(left_flank) >= len(seq):
            break
        base = seq[idx + len(left_flank)]
        # Confirm with right flank
        next_chunk = seq[idx + len(left_flank) + 1: idx + len(left_flank) + 1 + len(right_flank)]
        score = sum(1 for a,b in zip(next_chunk, right_flank.upper()) if a==b)
        if score >= len(right_flank)*0.85:  # 85% right-flank similarity
            fwd_matches.append((idx, base, score))
        start = idx + 1

    # RC search
    rc_matches = []
    start = 0
    while True:
        idx = seq.find(right_rc.upper(), start)
        if idx < 0: break
        if idx + len(right_rc) >= len(seq):
            break
        base = seq[idx + len(right_rc)]
        next_chunk = seq[idx + len(right_rc) + 1: idx + len(right_rc) + 1 + len(left_rc)]
        score = sum(1 for a,b in zip(next_chunk, left_rc.upper()) if a==b)
        if score >= len(left_rc)*0.85:
            base_complement = rc.get(base, 'N')
            rc_matches.append((idx, base_complement, score))
        start = idx + 1

    all_matches = fwd_matches + rc_matches
    if len(all_matches) == 0:
        return ('NOT_FOUND', 0, None)
    if len(all_matches) > 1:
        # Multiple hits — pick best by right-flank score
        all_matches.sort(key=lambda t: -t[2])
    return (all_matches[0][1], len(all_matches), all_matches[0][0])

# For each probe, query each founder FASTA
print('\n=== FASTA tie-breaker results ===\n')
all_results = []
for _, p in top5.iterrows():
    pos = int(p['pos']); ref = str(p['ref']); alt = str(p['alt'])
    print(f'--- Chr1:{pos} REF={ref} ALT={alt}  (cactus AC={p["AC_cactus"]}/{p["AN_cactus"]}, '
          f'xwu AC={p["AC_xwu_any"]}/{p["AN_xwu"]}) ---')
    left, right = get_tair_flank(pos, ref, alt, flank=80)

    # Sample a subset of founders for tractability — pick 30 random ones
    rng = np.random.default_rng(seed=pos)
    sampled = rng.choice(founder_fastas, size=min(30, len(founder_fastas)), replace=False)

    counts = {'REF': 0, 'ALT': 0, 'OTHER': 0, 'NOT_FOUND': 0, 'MULTI': 0}
    for fa in sampled:
        asm_id = fa.stem.replace('.chr','')
        base, n_hits, loc = search_founder_fasta(fa, left, right)
        if n_hits == 0:
            counts['NOT_FOUND'] += 1
        elif n_hits > 1:
            counts['MULTI'] += 1   # ambiguous — multiple flank hits
            # but still record the best
            if base == ref: counts['REF'] += 1
            elif base == alt: counts['ALT'] += 1
            else: counts['OTHER'] += 1
        else:
            if base == ref: counts['REF'] += 1
            elif base == alt: counts['ALT'] += 1
            else: counts['OTHER'] += 1
    n_called = counts['REF'] + counts['ALT'] + counts['OTHER']
    print(f'   founders sampled: 30')
    print(f'     REF ({ref}):      {counts["REF"]:>3}  (cactus would say GT=0)')
    print(f'     ALT ({alt}):      {counts["ALT"]:>3}  (xwu/hapFIRE would say GT=1)')
    print(f'     OTHER base:       {counts["OTHER"]:>3}')
    print(f'     not_found:        {counts["NOT_FOUND"]:>3}')
    print(f'     multi-hit:        {counts["MULTI"]:>3}')
    if n_called > 0:
        af_fasta = counts['ALT'] / n_called
    else:
        af_fasta = float('nan')
    print(f'     ⇒ FASTA AF (ALT/called): {af_fasta:.3f}    '
          f'(vs cactus recipe {p["AC_cactus"]/p["AN_cactus"]:.3f}, '
          f'xwu {p["AC_xwu_any"]/p["AN_xwu"]:.3f})')
    all_results.append({
        'pos': pos, 'ref': ref, 'alt': alt,
        'AF_cactus': p['AC_cactus']/p['AN_cactus'],
        'AF_xwu':    p['AC_xwu_any']/p['AN_xwu'],
        'AF_hapfire':p['hapfire_af'],
        'AF_FASTA':  af_fasta,
        'n_REF':     counts['REF'],
        'n_ALT':     counts['ALT'],
        'n_OTHER':   counts['OTHER'],
        'n_not_found': counts['NOT_FOUND'],
    })
    print()

R = pd.DataFrame(all_results)
print('\n=== Summary table ===')
print(R.round(3).to_string(index=False))
print('\n=== Verdict ===')
print('  At each probe, who is closer to FASTA-derived AF?')
for _, r in R.iterrows():
    d_cac = abs(r['AF_cactus'] - r['AF_FASTA'])
    d_xwu = abs(r['AF_xwu']    - r['AF_FASTA'])
    winner = 'cactus' if d_cac < d_xwu else 'xwu/hapfire'
    print(f'  Chr1:{r["pos"]:>10}  AF_FASTA={r["AF_FASTA"]:.3f}  '
          f'AF_cactus={r["AF_cactus"]:.3f}  AF_xwu={r["AF_xwu"]:.3f}  → closer: {winner}')

R.to_csv(ROOT / 'scratch/f2_honest/celld_fasta_tiebreak.tsv', sep='\t', index=False)
