#!/bin/bash
#SBATCH --job-name=p231n_fa
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=2:00:00
#SBATCH --array=1-231
#SBATCH --requeue
#SBATCH --output=logs/05n_fa_%A_%a.out
#SBATCH --error=logs/05n_fa_%A_%a.out
# =============================================================================
# Per-founder consensus FASTAs from the arch3 231 VCF with ./. written as N (bcftools consensus -M N), not REF.
# fastas_231/ (05_build_fastas_p231.sh) writes REF at ./., but every kMate panel these benchmarks score with
# (data/kmer_pa_231_arch3_filt2inv, data/kmer_pa_p231{,_filt2inv}) is built --treat-missing-as-n: a founder gets no
# k-mer at its ./. records. Reads from fastas_231 therefore carry REF k-mers the panel cannot place on the true
# founder, and the EM moves frequency onto other founders (see p80/scripts/05n_build_fastas_p78n.sh, the same fix
# for the long-read benchmark). With -M N the simulated genomes and the panels agree.
# Output: fastas_231n/<Accession_ID>.chr.fa (+ .fai, + .chain: TAIR10 -> consensus coordinates for make_recomb_mosaics.py). The *_n scripts (06n, 06bn, 07gn, 07hn, 08n, 09n) use them.
# =============================================================================
set -euo pipefail
mkdir -p logs
CTRL=/global/scratch/users/tbellg/kmate/benchmarks/p231
VCF=/global/scratch/users/tbellg/kmate/panel/arch3/chr1/merged_231_chr1_final.vcf.gz
REF=/global/scratch/users/tbellg/pang/pang_1001gplus/20260209_Exposito-Alonso/chr_only/TAIR10.chr.iupacN.fa
BCF=/global/home/users/tbellg/miniforge3/envs/kmate/bin/bcftools
SAMTOOLS=/global/home/users/tbellg/miniforge3/envs/kmate/bin/samtools
OUT_DIR=$CTRL/fastas_231n; mkdir -p $OUT_DIR
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $CTRL/data/founders_231_order.txt)
[[ -z "$SAMPLE" ]] && { echo "no sample at index $SLURM_ARRAY_TASK_ID" >&2; exit 1; }
OUT_FA=$OUT_DIR/${SAMPLE}.chr.fa
[ -s "${OUT_FA}.fai" ] && { echo "already done -- $OUT_FA"; exit 0; }
$BCF consensus -f $REF -H 1 -M N -s $SAMPLE -c $OUT_DIR/${SAMPLE}.chain $VCF 2> ${OUT_FA}.consensus.log > ${OUT_FA}.tmp
mv ${OUT_FA}.tmp $OUT_FA
$SAMTOOLS faidx $OUT_FA
echo "[$(date)] DONE -- $OUT_FA"; tail -1 ${OUT_FA}.consensus.log
