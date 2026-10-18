#!/bin/bash
#SBATCH --job-name=p231_chain
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=1:00:00
#SBATCH --array=1-231
#SBATCH --requeue
#SBATCH --output=logs/05c_chain_%A_%a.out
# =============================================================================
# TAIR10 -> consensus chains for the fastas_231n founder genomes (made before 05n wrote chains), with the exact
# 05n_build_fastas_p231n.sh command, and a check that the regenerated consensus is byte-identical to the FASTA the sims
# use. make_recomb_mosaics.py needs the chain to cut each segment at the founder's own coordinates
# (benchmarks/archive/scripts_before_simulations_2026-10-17/p80/audit_sims/). Output: fastas_231n/<Accession_ID>.chain (+ .chaincheck: SAME / DIFFERENT).
# =============================================================================
set -euo pipefail
mkdir -p logs
CTRL=/global/scratch/users/tbellg/kmate/benchmarks/p231
VCF=/global/scratch/users/tbellg/kmate/panel/arch3/chr1/merged_231_chr1_final.vcf.gz
REF=/global/scratch/users/tbellg/pang/pang_1001gplus/20260209_Exposito-Alonso/chr_only/TAIR10.chr.iupacN.fa
BCF=/global/home/users/tbellg/miniforge3/envs/kmate/bin/bcftools
SAMTOOLS=/global/home/users/tbellg/miniforge3/envs/kmate/bin/samtools
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $CTRL/data/founders_231_order.txt)
OUT=$CTRL/fastas_231n/$SAMPLE
TMP=$CTRL/fastas_231n/.tmp_chain_$SAMPLE.fa
$BCF consensus -f $REF -H 1 -M N -s $SAMPLE -c $OUT.chain $VCF 2>/dev/null > $TMP
a=$($SAMTOOLS faidx $TMP Chr1 | tail -n +2 | tr -d '\n' | md5sum | cut -c1-32)
b=$($SAMTOOLS faidx $OUT.chr.fa Chr1 | tail -n +2 | tr -d '\n' | md5sum | cut -c1-32)
echo -e "$SAMPLE\t$([ "$a" = "$b" ] && echo SAME || echo DIFFERENT)" > $OUT.chaincheck
rm -f $TMP $TMP.fai
