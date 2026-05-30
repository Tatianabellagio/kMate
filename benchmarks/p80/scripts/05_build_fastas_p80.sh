#!/bin/bash
#SBATCH --job-name=p80_a5_fa
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=4:00:00
#SBATCH --array=1-80
#SBATCH --requeue
#SBATCH --output=logs/05_fa_%A_%a.out
#SBATCH --error=logs/05_fa_%A_%a.err

# =============================================================================
# Phase A5 -- Build per-founder consensus FASTAs from the canonical p80 VCF.
#
# Each task builds one founder's Chr1 consensus by bcftools consensus -H 1.
# Sample names in the VCF are Accession_IDs (after A1's reheader).
# Output: fastas_80/<Accession_ID>.chr.fa  (+ .fai)
#
# DO NOT replace these by symlinks to v3's unimputed_fastas_v3/ or to raw cactus
# assemblies -- variant-set mismatch would recreate the v3 simulation bug.
# =============================================================================
mkdir -p logs
set -uo pipefail

CTRL=/global/scratch/users/tbellg/kmate/benchmarks/p80
VCF=$CTRL/data/pangenome_p80_chr1.vcf.gz
REF=/global/scratch/users/tbellg/pang/pang_1001gplus/20260209_Exposito-Alonso/chr_only/TAIR10.chr.iupacN.fa
OUT_DIR=$CTRL/fastas_80
SAMPLE_LIST=$CTRL/data/samples_80_acc_order.txt

BCF=/global/home/users/tbellg/miniforge3/envs/sequencing_pipeline/bin/bcftools
SAMTOOLS=/global/home/users/tbellg/miniforge3/envs/sequencing_pipeline/bin/samtools

mkdir -p $OUT_DIR

# Cache the canonical sample order (Accession_IDs) once.
[ -s "$SAMPLE_LIST" ] || $BCF query -l $VCF > $SAMPLE_LIST

SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $SAMPLE_LIST)
[[ -z "$SAMPLE" ]] && { echo "no sample at index $SLURM_ARRAY_TASK_ID" >&2; exit 1; }

echo "[$(date)] task=$SLURM_ARRAY_TASK_ID  sample=$SAMPLE"

OUT_FA=$OUT_DIR/${SAMPLE}.chr.fa
if [ -s "${OUT_FA}.fai" ]; then
    echo "  already done -- skip ($OUT_FA)"
    exit 0
fi

# bcftools consensus on the haploid VCF.
# `.` cells: bcftools consensus skips the variant (treats as REF). Matches v3
# convention; downstream kmer_pa/var_pa builders use the same `.` -> REF rule.
$BCF consensus -f $REF -H 1 -s $SAMPLE $VCF 2> ${OUT_FA}.consensus.log > $OUT_FA
$SAMTOOLS faidx $OUT_FA

echo "[$(date)] DONE -- $OUT_FA ($(du -h $OUT_FA | cut -f1))"
head -1 $OUT_FA
ls -l ${OUT_FA}.fai
