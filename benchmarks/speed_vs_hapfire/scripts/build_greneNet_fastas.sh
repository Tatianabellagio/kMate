#!/bin/bash
#SBATCH --job-name=gren_fastas
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=2:00:00
#SBATCH --output=logs/gren_fastas_%j.out
#SBATCH --error=logs/gren_fastas_%j.out

# Build 231 per-founder haploid consensus FASTAs from the greneNet_final_v1.1
# SNP panel (Chr1), for the FAIR hapFIRE benchmark: hapFIRE's simulated reads
# come from these greneNet-derived genomes and hapFIRE is tested against the
# same greneNet VCF -- self-consistent panel (mirror of kMate on arch3).
# greneNet is SNP-only + fully-called (0% missing), so consensus is clean.
# Founder IDs are IDENTICAL to arch3, and we reuse the arch3 founders-meta +
# same seeds downstream so the simulated POOLS match the arch3 (kMate) sims.
set -euo pipefail
mkdir -p logs
ROOT=/global/scratch/users/tbellg/kmate
BASE=$ROOT/benchmarks/speed_vs_hapfire
BCF=/global/home/users/tbellg/miniforge3/bin/bcftools
SAMTOOLS=/global/home/users/tbellg/miniforge3/envs/kmate/bin/samtools
TABIX=/global/home/users/tbellg/miniforge3/bin/tabix
REF_XING=/global/scratch/users/tbellg/pang/ref_xing/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa
VCF_IN=$BASE/work/greneNet_chr1.vcf.gz
OUTDIR=$BASE/greneNet_fastas
WORK=$BASE/work
mkdir -p "$OUTDIR"

# 1) rename VCF contig 1 -> Chr1 (match the sim pipeline's --chroms Chr1)
echo "1 Chr1" > "$WORK/rename_1_to_Chr1.txt"
VCF=$WORK/greneNet_chr1.Chr1.vcf.gz
if [ ! -s "$VCF" ]; then
    $BCF annotate --rename-chrs "$WORK/rename_1_to_Chr1.txt" "$VCF_IN" -Oz -o "$VCF"
    $TABIX -p vcf "$VCF"
fi

# 2) Chr1 reference with contig renamed 1 -> Chr1
REF=$WORK/tair10_Chr1.fa
if [ ! -s "$REF" ]; then
    $SAMTOOLS faidx "$REF_XING" 1 | sed 's/^>1.*/>Chr1/' > "$REF"
    $SAMTOOLS faidx "$REF"
fi

# 3) per-founder haploid consensus (-H 1 = first haplotype; founders ~fully
#    homozygous so this is deterministic and matches the assembly allele)
FOUNDERS=$($BCF query -l "$VCF")
n=0
for F in $FOUNDERS; do
    OUT=$OUTDIR/${F}.chr.fa
    if [ -s "$OUT.fai" ]; then n=$((n+1)); continue; fi
    $BCF consensus -H 1 -f "$REF" -s "$F" "$VCF" 2>/dev/null | sed 's/^>Chr1.*/>Chr1/' > "$OUT"
    $SAMTOOLS faidx "$OUT"
    n=$((n+1))
    [ $((n % 25)) -eq 0 ] && echo "[$(date)] $n founders done"
done
echo "[$(date)] DONE: $n founder FASTAs in $OUTDIR"
ls "$OUTDIR"/*.chr.fa | wc -l
