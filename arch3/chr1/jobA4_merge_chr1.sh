#!/bin/bash
#SBATCH --job-name=chr1_merge
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --output=logs/A4_merge_%j.out
#SBATCH --error=logs/A4_merge_%j.err
mkdir -p logs
set -euo pipefail

# Phase 2 Job A4: merge cactus_78 + PG_153 (both Chr1 haploid biallelic) + post-merge AN=0 filter.
# Outputs the final 231-panel for Chr1.

cd /global/scratch/users/tbellg/kmate/arch3/chr1

BCF=/global/home/users/tbellg/miniforge3/envs/gwas/bin/bcftools
TABIX=/global/home/users/tbellg/miniforge3/envs/gwas/bin/tabix
PY=/global/home/users/tbellg/miniforge3/envs/hapfm/bin/python

CACTUS_HAP=cactus_78_chr1_haploid.vcf.gz
PG_HAP=pg_153_chr1_haploid.vcf.gz

[ -s $CACTUS_HAP ] || { echo "ERROR: missing $CACTUS_HAP"; exit 1; }
[ -s $PG_HAP ]     || { echo "ERROR: missing $PG_HAP"; exit 1; }

echo "[$(date)] === Step 1: bcftools merge (both haploid biallelic) ==="
# --merge none: same-(CHROM,POS,REF,ALT) records combine sample sets; different ALTs become
# separate biallelic rows. Default would create multi-allelic records which build_cn_var
# misattributes (stores ALT[0] only but counts any non-zero allele as carrier of ALT[0]).
MERGED=merged_231_chr1.vcf.gz
$BCF merge --merge none $CACTUS_HAP $PG_HAP --threads 8 -Oz -o $MERGED
$TABIX -p vcf $MERGED
N_MERGED=$($BCF view -H $MERGED | wc -l)
N_MULTI=$($BCF view -H $MERGED | awk -F'\t' '$5 ~ /,/' | wc -l)
echo "  merged records: $N_MERGED"
echo "  multi-allelic rows (should be 0): $N_MULTI"
echo "  samples: $($BCF query -l $MERGED | wc -l) (expected 231)"

echo
echo "[$(date)] === Step 2: fill-tags ==="
MERGED_FILLED=merged_231_chr1_filled.vcf.gz
$BCF +fill-tags $MERGED --threads 8 -Oz -o $MERGED_FILLED -- -t AC,AN,F_MISSING
$TABIX -p vcf $MERGED_FILLED

echo
echo "[$(date)] === Step 3: post-merge AN=0 filter ==="
MERGED_FINAL=merged_231_chr1_final.vcf.gz
N_PRE=$($BCF view -H $MERGED_FILLED | wc -l)
$BCF view -e 'INFO/AN=0' $MERGED_FILLED --threads 8 -Oz -o $MERGED_FINAL
$TABIX -p vcf $MERGED_FINAL
N_POST=$($BCF view -H $MERGED_FINAL | wc -l)
echo "  records pre-AN=0:  $N_PRE"
echo "  records post-AN=0: $N_POST"
echo "  dropped: $((N_PRE - N_POST))"

echo
echo "[$(date)] === Step 4: panel-level AC at spot-checks + F_MISSING dist ==="
$PY <<'PYEOF'
import gzip
SPOT = [(5870018,'T','A'), (10421645,'T','C'),
        (13843898,'C','T'), (13843898,'CT','TC'), (13843898,'CT','TG')]
results = {}
samples = None
fm_buckets = [0]*6
fm_labels = ['F=0', '0<F<=0.05', '0.05<F<=0.1', '0.1<F<=0.3', '0.3<F<=0.5', 'F>0.5']
total = 0
with gzip.open('merged_231_chr1_final.vcf.gz','rt') as f:
    for line in f:
        if line.startswith('##'):
            continue
        if line.startswith('#CHROM'):
            samples = line.rstrip().split('\t')[9:]
            continue
        parts = line.rstrip().split('\t')
        # extract F_MISSING
        fm = None
        for kv in parts[7].split(';'):
            if kv.startswith('F_MISSING='):
                try: fm = float(kv[10:])
                except: fm = None
                break
        if fm is not None:
            total += 1
            if fm == 0: fm_buckets[0] += 1
            elif fm <= 0.05: fm_buckets[1] += 1
            elif fm <= 0.1: fm_buckets[2] += 1
            elif fm <= 0.3: fm_buckets[3] += 1
            elif fm <= 0.5: fm_buckets[4] += 1
            else: fm_buckets[5] += 1
        if (int(parts[1]), parts[3], parts[4]) in SPOT:
            results[(int(parts[1]), parts[3], parts[4])] = parts

print(f'\nSamples: {len(samples)}')
print(f'Total Chr1 records: {total:,}')
print()
print('=== F_MISSING distribution ===')
for lbl, n in zip(fm_labels, fm_buckets):
    print(f'  {lbl:>14}: {n:>10,} ({100*n/max(total,1):.2f}%)')

print()
print('=== Spot-check carrier counts ===')
for k in SPOT:
    if k not in results:
        print(f'  Chr1:{k[0]} {k[1]}>{k[2]}: not present')
        continue
    p = results[k]
    gts = p[9:]
    ac = sum(1 for g in gts if g == '1')
    an = sum(1 for g in gts if g in ('0','1'))
    print(f'  Chr1:{k[0]} {k[1]}>{k[2]}: AC={ac}/{an}  AF={ac/max(an,1):.4f}')
PYEOF

echo
echo "[$(date)] DONE A4 (merged Chr1 panel)"
ls -lh merged_231_chr1_final.vcf.gz
