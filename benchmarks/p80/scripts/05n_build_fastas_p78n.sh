#!/bin/bash
#SBATCH --job-name=p78n_fa
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=2:00:00
#SBATCH --array=1-78
#SBATCH --requeue
#SBATCH --output=logs/05n_fa_%A_%a.out
#SBATCH --error=logs/05n_fa_%A_%a.out
# =============================================================================
# Per-founder consensus FASTAs for the long-read-only benchmark, with ./. written as N (bcftools consensus -M N) instead
# of REF. fastas_80/ (05_build_fastas_p80.sh) writes REF at ./., but the production k-mer panel (kmer_pa_p80_ours,
# --treat-missing-as-n) gives a founder no k-mer at its ./. records -- so reads from fastas_80 carry REF k-mers the panel
# cannot place on the true founder (~1% of k-mer counts per founder), which the EM pushes onto other founders (3.6% of
# frequency misassigned in a 1-founder pool). With -M N the simulated genome and the panel agree: no reads / no k-mers
# at a founder's ./. records, as for the panel. Output: fastas_78n/<Accession_ID>.chr.fa (+ .fai), the 78 founders.
# Also writes fastas_78n/<Accession_ID>.chain (TAIR10 -> consensus coordinates), which make_recomb_mosaics.py needs to cut
# recombinant segments at each founder's own coordinates.
# =============================================================================
set -euo pipefail
mkdir -p logs
CTRL=/global/scratch/users/tbellg/kmate/benchmarks/p80
VCF=$CTRL/data/pangenome_p80_chr1.vcf.gz
REF=/global/scratch/users/tbellg/pang/pang_1001gplus/20260209_Exposito-Alonso/chr_only/TAIR10.chr.iupacN.fa
BCF=/global/home/users/tbellg/miniforge3/envs/kmate/bin/bcftools
SAMTOOLS=/global/home/users/tbellg/miniforge3/envs/kmate/bin/samtools
OUT_DIR=$CTRL/fastas_78n; mkdir -p $OUT_DIR
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $CTRL/data/samples_78_acc.txt)
[[ -z "$SAMPLE" ]] && { echo "no sample at index $SLURM_ARRAY_TASK_ID" >&2; exit 1; }
OUT_FA=$OUT_DIR/${SAMPLE}.chr.fa
[ -s "${OUT_FA}.fai" ] && { echo "already done -- $OUT_FA"; exit 0; }
$BCF consensus -f $REF -H 1 -M N -s $SAMPLE -c $OUT_DIR/${SAMPLE}.chain $VCF 2> ${OUT_FA}.consensus.log > ${OUT_FA}.tmp
mv ${OUT_FA}.tmp $OUT_FA
$SAMTOOLS faidx $OUT_FA
echo "[$(date)] DONE -- $OUT_FA"; cat ${OUT_FA}.consensus.log
