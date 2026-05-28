#!/bin/bash
#SBATCH --job-name=fasta_v3
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=1:00:00
#SBATCH --output=logs/fasta_v3_%A_%a.out
#SBATCH --error=logs/fasta_v3_%A_%a.err

# =============================================================================
# build_consensus_fastas_one.sh
#
# For one PanGenie founder (selected by SLURM_ARRAY_TASK_ID), build their
# unimputed per-founder consensus FASTA by applying the haploid VCF onto
# TAIR10 with `bcftools consensus -H A`.
#
# Symlinks the 80 cactus-assembly FASTAs are expected to be pre-created by
# the orchestrator (launch_v3_rebuild.sh) — this array only handles the 151
# PanGenie-genotyped ecotypes.
#
# Output: sims/visor_freqk/founder_fastas_231_v3/<eco>.chr.fa  (+ .fai)
# =============================================================================
set -euo pipefail

BASE=/global/scratch/users/tbellg/hapfire_sv
BCF=/global/home/users/tbellg/miniforge3/envs/sequencing_pipeline/bin/bcftools
SAMTOOLS=/global/home/users/tbellg/miniforge3/envs/sequencing_pipeline/bin/samtools

VCF=$BASE/pangenie_genotyping/data/merged/founders_231_chr.haploid.vcf.gz
REF=/global/scratch/users/tbellg/pang/pang_1001gplus/20260209_Exposito-Alonso/chr_only/TAIR10.chr.iupacN.fa
ECO_LIST=$BASE/pangenie_genotyping/data/v3/pangenie_151.txt
OUT_DIR=$BASE/sims/visor_freqk/founder_fastas_231_v3
mkdir -p $OUT_DIR $BASE/pangenie_genotyping/logs

IDX=${SLURM_ARRAY_TASK_ID:-1}
ECO=$(sed -n "${IDX}p" "$ECO_LIST")
[ -n "$ECO" ] || { echo "ERROR: no ecotype at task index $IDX (list has $(wc -l < $ECO_LIST) lines)" >&2; exit 1; }

OUT=$OUT_DIR/${ECO}.chr.fa
if [ -s "${OUT}.fai" ] && [ -s "$OUT" ]; then
    echo "[$(date)] $ECO: already built ($OUT)"
    exit 0
fi

echo "[$(date)] $ECO: bcftools consensus -H A → $OUT"
# -H A: take the first allele (haploid: that's the only allele anyway)
# -s $ECO: restrict to this sample
# Missing GTs ('.') default to REF unless --missing is set; we keep default
# behavior — '.' cells in the haploid VCF treat the founder as having ref
# at that position, which matches the carrier-status convention.
$BCF consensus -f "$REF" -s "$ECO" -H A "$VCF" > "$OUT" 2> "${OUT}.consensus.log"

# Sanity: bcftools consensus emits "Applied <N> variants" on stderr; if the
# chrom names don't match, N=0 and the FASTA is just TAIR10 unchanged. Fail
# loud so a chrom-name regression doesn't pass silently again.
APPLIED=$(grep -oE "Applied [0-9]+ variants" "${OUT}.consensus.log" | head -1 | awk '{print $2}' || echo "")
if [ -z "$APPLIED" ] || [ "$APPLIED" = "0" ]; then
    echo "ERROR: bcftools consensus applied $APPLIED variants for $ECO — chrom-name mismatch?" >&2
    cat "${OUT}.consensus.log" >&2
    rm -f "$OUT"
    exit 1
fi
echo "  applied $APPLIED variants"

# Index the FASTA so downstream pysam.FastaFile.fetch works
$SAMTOOLS faidx "$OUT"

echo "[$(date)] $ECO: DONE  size=$(du -h "$OUT" | cut -f1)"
ls -lh "$OUT" "${OUT}.fai"
