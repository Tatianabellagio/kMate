#!/bin/bash
#SBATCH --job-name=p231_hf_psd
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=2:00:00
#SBATCH --output=logs/hapfire_psd_%j.out
#SBATCH --error=logs/hapfire_psd_%j.out

# =============================================================================
# hapFIRE-vs-kMate poolsize x depth sweep, p231 arm ONLY (hapFIRE's panel is
# greneNet -- the 231-founder SNP set kMate's p231 panel is built from; there
# is no equivalent SNP panel for p80).
#
# kMate consumes raw fastq directly (alignment-free); hapFIRE requires a BAM,
# so unlike the original single-point speed_vs_hapfire benchmark (which reused
# a prebuilt BAM from a different sim pipeline), we align fresh here from the
# SAME reads kMate's p231 poolsize_depth sweep uses. Both align (bwa mem) and
# hapFIRE (Phase 1-2: independent-LD partition + HARP like/freq, Phase 4/BigLD
# skipped exactly as in the original benchmark) are timed separately with
# /usr/bin/time -v so alignment cost can be included or excluded from the
# "speed" comparison as needed.
#
# Usage: sbatch run_hapfire_poolsize_depth.sh N COV SEED
#   N    = 2 | 5 | 20 | 50 | 150
#   COV  = 1 | 10
#   SEED = 42 | 43 | 44 | 45 | 46
# =============================================================================
set -euo pipefail
N=${1:?Usage: N COV SEED}
COV=${2:?Usage: N COV SEED}
SEED=${3:?Usage: N COV SEED}

ROOT=/global/scratch/users/tbellg/kmate
P231=$ROOT/benchmarks/p231
WORK=$ROOT/benchmarks/speed_vs_hapfire/work
OUT=$ROOT/benchmarks/speed_vs_hapfire/results/poolsize_depth/n${N}_cov${COV}_s${SEED}
mkdir -p "$OUT" logs

BWA=/global/home/users/tbellg/miniforge3/envs/kmate/bin/bwa
SAMTOOLS=/global/home/users/tbellg/miniforge3/envs/kmate/bin/samtools
HFPY=/global/home/users/tbellg/miniforge3/envs/hapfire/bin/python
HF=$ROOT/external/HapFIRE
REF=/global/scratch/users/tbellg/pang/ref_xing/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa
VCF=$WORK/greneNet_chr1.vcf.gz

READS=$P231/sims/cov${COV}_n${N}_g0_s${SEED}_hotspots_p231_chr1/reads
[ -s "$READS/r1.fq" ] || { echo "ERROR: missing $READS/r1.fq" >&2; exit 1; }
[ -s "${REF}.bwt" ]    || { echo "ERROR: missing bwa index $REF.bwt" >&2; exit 1; }

SAMPLE=p231_hapfire_psd_n${N}_cov${COV}_s${SEED}
cd "$OUT"

echo "[$(date)] ALIGN N=$N cov=$COV s=$SEED  (host=$(hostname), cpus=${SLURM_CPUS_PER_TASK:-8})"
/usr/bin/time -v -o "$OUT/align_time.txt" bash -c "
    set -euo pipefail
    $BWA mem -t ${SLURM_CPUS_PER_TASK:-8} $REF $READS/r1.fq $READS/r2.fq 2> $OUT/bwa.log \
        | $SAMTOOLS sort -@ ${SLURM_CPUS_PER_TASK:-8} -o $OUT/sim.srt.bam -
    $SAMTOOLS index $OUT/sim.srt.bam
"

# reheader Chr1 -> 1 so contigs match hapFIRE's greneNet VCF (Ensembl-style)
$SAMTOOLS view -H "$OUT/sim.srt.bam" | sed 's/SN:Chr1/SN:1/' > "$OUT/newhdr.sam"
$SAMTOOLS reheader "$OUT/newhdr.sam" "$OUT/sim.srt.bam" > "$OUT/sim.chr1named.bam"
$SAMTOOLS index "$OUT/sim.chr1named.bam"
rm -f "$OUT/sim.srt.bam" "$OUT/sim.srt.bam.bai" "$OUT/newhdr.sam"

echo "[$(date)] HAPFIRE N=$N cov=$COV s=$SEED"
export PATH="$HF/src:$PATH"
# -haplotype deliberately omitted (see run_hapfire.sh note: argparse bool bug
# makes "-haplotype False" evaluate True) -> Phase 4/BigLD stays skipped.
/usr/bin/time -v -o "$OUT/hapfire_time.txt" $HFPY "$HF/hapFIRE.py" \
    -vcf "$VCF" \
    -bam "$OUT/sim.chr1named.bam" \
    -fa  "$REF" \
    -output $SAMPLE

echo "[$(date)] DONE N=$N cov=$COV s=$SEED"
ls -la "$OUT"
