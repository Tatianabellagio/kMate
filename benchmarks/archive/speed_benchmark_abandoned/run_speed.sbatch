#!/bin/bash
#SBATCH --job-name=speed
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --exclusive
#SBATCH --time=06:00:00
#SBATCH --output=benchmarks/speed_benchmark/logs/%x_%j.out
#SBATCH --error=benchmarks/speed_benchmark/logs/%x_%j.out
set -eo pipefail
# CLEAN speed/CPU benchmark: ONE tool, ALONE on an exclusive node, matched 8 threads,
# end-to-end fastq->AF, /usr/bin/time -v (wall + user/sys CPU + %CPU + peak RSS).
# Usage: sbatch run_speed.sbatch <TOOL> <POOL> <PANEL> [VG_PRESET]
#   TOOL = kmate | hapfire | vg
cd /global/scratch/users/tbellg/kmate
TOOL=${1:?kmate|hapfire|vg}; POOL=${2:?pool}; PANEL=${3:?p80|p231}; VGPRE=${4:-}
T=8                                   # matched thread count for all tools
export OMP_NUM_THREADS=$T OPENBLAS_NUM_THREADS=$T MKL_NUM_THREADS=$T
mkdir -p benchmarks/speed_benchmark/logs benchmarks/speed_benchmark/results
source ~/miniforge3/etc/profile.d/conda.sh
SC=benchmarks/accuracy_vs_competitors/scripts
W=benchmarks/accuracy_vs_competitors/work
RES=benchmarks/speed_benchmark/results
D=benchmarks/$PANEL/sims/$POOL
REF=/global/scratch/users/tbellg/pang/pang_1001gplus/20260209_Exposito-Alonso/chr_only/TAIR10.chr.iupacN.fa
TF=$RES/speed_${TOOL}_${POOL}${VGPRE:+_$VGPRE}_time.txt
echo "[$(date)] SPEED $TOOL pool=$POOL panel=$PANEL host=$(hostname) threads=$T (exclusive node)"

case "$PANEL" in
  p80)  KP=benchmarks/p80/data/kmer_pa_p80_filt2/kmer_pa
        VP=benchmarks/p80/data/var_pa_p80.var_pa.npz; VM=benchmarks/p80/data/var_pa_p80.meta.npz ;;
  p231) KP=data/kmer_pa_231_arch3_filt2inv/kmer_pa
        VP=panel/arch3/chr1/var_pa_231_arch3_chr1.var_pa.npz; VM=panel/arch3/chr1/var_pa_231_arch3_chr1.meta.npz ;;
esac

case "$TOOL" in
  kmate)   # end-to-end: jellyfish k-mer count + EM (global) -> AF
    conda activate kmate
    /usr/bin/time -v -o "$TF" python -m kmate.per_sample_per_chrom \
      --kmer-pa-prefix "$KP" --var-pa "$VP" --var-meta "$VM" \
      --reads "$D/reads/r1.fq" "$D/reads/r2.fq" \
      --sample speed_${POOL} --threads "$T" --chroms Chr1 \
      --block-mode global --kmer-weight uniform \
      --out $RES/kmate_${POOL}.tsv ;;
  hapfire) # end-to-end: align reads to TAIR10 + HARP/hapFIRE -> per-SNP AF
    conda activate kmate    # for minimap2/samtools
    ROOT=$(pwd)             # absolute repo root (avoid $PWD-in-subshell expansion bugs)
    HF=$ROOT/external/HapFIRE; export PATH="$HF/src:$PATH"
    PYHF=/global/home/users/tbellg/miniforge3/envs/hapfire/bin/python
    BAM=$ROOT/$RES/${POOL}.hf.bam
    /usr/bin/time -v -o "$TF" bash -c "
      minimap2 -ax sr --MD --cs -Y --sam-hit-only -t $T -R '@RG\tID:illumina\tSM:bulk' \
        $REF $ROOT/$D/reads/r1.fq $ROOT/$D/reads/r2.fq | samtools sort -@ $T -o $BAM -
      samtools index $BAM
      cd $ROOT/$RES && $PYHF $HF/hapFIRE.py \
        -vcf $ROOT/$W/shared_snps_${PANEL}_Chr1.vcf.gz \
        -bam $BAM -fa $REF -output speed_hapfire_${POOL}" ;;
  vg)      # end-to-end: giraffe + surject + AF
    conda activate vg
    /usr/bin/time -v -o "$TF" bash $SC/run_vg.sh "$POOL" "$PANEL" "$T" "$VGPRE" ;;
esac
echo "[$(date)] DONE. timing -> $TF"
grep -iE "wall clock|User time|System time|Percent of CPU|Maximum resident" "$TF" 2>/dev/null
