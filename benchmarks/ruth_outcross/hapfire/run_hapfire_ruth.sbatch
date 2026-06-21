#!/bin/bash
#SBATCH --job-name=hapfire_ruth
#SBATCH --account=fc_moilab
#SBATCH --partition=savio3
#SBATCH --qos=savio_normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=06:00:00
#SBATCH --output=benchmarks/ruth_outcross/hapfire/logs/hf_%A_%a.out
#SBATCH --error=benchmarks/ruth_outcross/hapfire/logs/hf_%A_%a.err

# ============================================================================
# Run hapFIRE on a few of Ruth's individually-sequenced plants, so we can compare
# hapFIRE's per-founder _ecotype_frequency vector to kMate's h on the SAME panel.
#
# hapFIRE estimates ECOTYPE (founder) frequencies from SNP base-matches; for a
# single individual that vector is hapFIRE's analog of kMate's h:
#   inbred   -> one ecotype ~ 1.0
#   outcross -> spread over >=2 ecotypes
#
# Shared panel: the 231-founder Chr1 VCF kMate's panel was built from
# (greneNet_chr1.vcf.gz, chrom "1"); Ruth's dedup BAMs are aligned to the Ensembl
# TAIR10 ref with chrom names "1".."5" -> they match the VCF + ref with NO
# renaming (verified). Chr1 only (the hapFIRE-ready VCF is Chr1).
#
# Phases 1-2 only (-haplotype omitted on purpose; argparse type=bool makes
# "-haplotype False" == True -> would wrongly enable the BigLD/R phase).
#
# Usage: sbatch --array=1-4 benchmarks/ruth_outcross/hapfire/run_hapfire_ruth.sbatch
# ============================================================================
set -uo pipefail
ROOT=/global/scratch/users/tbellg/kmate
cd "$ROOT"
B=$ROOT/benchmarks/ruth_outcross/hapfire
mkdir -p "$B/logs" "$B/results"

HF=$ROOT/external/HapFIRE
PY=/global/home/users/tbellg/miniforge3/envs/hapfire/bin/python
export PATH="$HF/src:$PATH"          # harp binary

VCF=$ROOT/benchmarks/speed_vs_hapfire/work/greneNet_chr1.vcf.gz
REF=/global/scratch/users/tbellg/pang/ref_xing/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa
MANIFEST=$B/samples.tsv

LINE=$((SLURM_ARRAY_TASK_ID + 1))    # +1 skip header
ROW=$(awk -F'\t' -v n=$LINE 'NR==n' "$MANIFEST")
SAMPLE=$(echo "$ROW" | cut -f1)
BAM=$(echo "$ROW" | cut -f2)
LABEL=$(echo "$ROW" | cut -f3)
[ -n "$SAMPLE" ] || { echo "ERROR: no sample at line $LINE"; exit 1; }
[ -s "$BAM" ]    || { echo "ERROR: BAM missing $BAM"; exit 1; }

OUT=$B/results/${SAMPLE}
if [ -s "${OUT}_ecotype_frequency.txt" ]; then
    echo "[$(date)] $SAMPLE ecotype_frequency already present — skipping"; exit 0
fi

echo "[$(date)] hapFIRE $SAMPLE ($LABEL) on $(hostname), ${SLURM_CPUS_PER_TASK:-?} cpus"
echo "  VCF=$VCF"; echo "  BAM=$BAM"; echo "  harp=$(command -v harp)"
cd "$B/results"
/usr/bin/time -v -o "${OUT}_time.txt" $PY "$HF/hapFIRE.py" \
    -vcf "$VCF" -bam "$BAM" -fa "$REF" -output "$SAMPLE" \
  || { echo "[$(date)] ERROR: hapFIRE failed $SAMPLE"; exit 1; }
echo "[$(date)] DONE $SAMPLE"
ls -la "${OUT}"*
