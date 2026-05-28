#!/bin/bash
#SBATCH --job-name=p231_a5_fa
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=4:00:00
#SBATCH --array=1-231
#SBATCH --requeue
#SBATCH --output=logs/05_fa_%A_%a.out
#SBATCH --error=logs/05_fa_%A_%a.err

# =============================================================================
# control_p231 Phase A5 -- per-founder consensus FASTAs from the ARCH3 canonical
# 231 VCF (merged_231_chr1_final.vcf.gz). Mirrors control_p80/05_build_fastas_p80.sh.
#
# Each task builds one founder's Chr1 consensus by bcftools consensus -H 1
# (VCF is haploid: GT cells are '.'/'0'). Sample names are Accession_IDs, order
# IDENTICAL to the cn_var/cn_full founders axis (verified).
# Output: fastas_231/<Accession_ID>.chr.fa (+ .fai)
#
# DO NOT symlink to v3/v3qc unimputed_fastas_* -- variant-set mismatch would
# recreate the v3 simulation bug (memory/project_v3_singleton_kmer_bug).
# =============================================================================
mkdir -p logs
set -uo pipefail

CTRL=/global/scratch/users/tbellg/kmate/control_p231
VCF=/global/scratch/users/tbellg/kmate/panel/arch3/chr1/merged_231_chr1_final.vcf.gz
REF=/global/scratch/users/tbellg/pang/pang_1001gplus/20260209_Exposito-Alonso/chr_only/TAIR10.chr.iupacN.fa
OUT_DIR=$CTRL/fastas_231
SAMPLE_LIST=$CTRL/data/founders_231_order.txt

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

# bcftools consensus on the haploid VCF. '.' cells -> REF (same rule as the
# cn_full/cn_var builders), so reads carry exactly the arch3 variant set.
$BCF consensus -f $REF -H 1 -s $SAMPLE $VCF 2> ${OUT_FA}.consensus.log > $OUT_FA
$SAMTOOLS faidx $OUT_FA

echo "[$(date)] DONE -- $OUT_FA ($(du -h $OUT_FA | cut -f1))"
head -1 $OUT_FA
ls -l ${OUT_FA}.fai
