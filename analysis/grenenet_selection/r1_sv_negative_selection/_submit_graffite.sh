#!/bin/bash
# GraffiTE, annotation + TSD only. We already have the pangenome graph and the genotypes, so the two
# halves of the pipeline that build and genotype a graph are skipped (--genotype false); --vcf enters the
# workflow directly at the RepeatMasker step.
#
# What this adds over the RepeatMasker run we already did:
#   OneCodeToFindThemAll   stitches fragmented hits into single copies (the 5' LTR / internal / 3' LTR
#                          problem), so family identity is per COPY rather than per fragment
#   -s sensitivity         GraffiTE runs RepeatMasker in slow mode, which recovers older, more diverged
#                          copies that the default setting misses
#   deletions              the REF sequence of a deletion is annotated exactly like the ALT of an
#                          insertion, so cargo exists for all 227,851 SV records, not only insertions
#   TSDs                   target site duplications: direct, library-independent evidence that an
#                          insertion arose by transposition
#
# NOTE on what the span means. GraffiTE leaves RepeatMasker's low-complexity and small-RNA masking ON, so
# its `total_match_span` is "covered by any repeat", not "covered by a TE" -- a microsatellite expansion
# can pass its 0.80 filter. Every hit carries its class, so both numbers are recoverable at parsing: the
# standard span, and a TE-only span excluding Simple_repeat / Low_complexity / Satellite / *RNA.
#SBATCH --job-name=graffite
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=32 --mem=128G
#SBATCH --time=72:00:00
#SBATCH --output=logs/graffite_%A.out
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
module load nextflow/24.10.4
T=/global/scratch/projects/fc_moilab/tbellg/kmate/telib
export NXF_SINGULARITY_CACHEDIR=$T/nxf_cache
export APPTAINER_CACHEDIR=$T/.apptainer_cache
export NXF_OPTS="-Xms2g -Xmx8g"
OUT=$PWD/results/sv_adaptive/graffite/run
mkdir -p "$OUT" "$T/nxf_work"
nextflow run $T/graffite_src/main.nf \
  --vcf  "$PWD/results/sv_adaptive/graffite/sv_sites.vcf" \
  --reference /global/scratch/users/tbellg/pang/pang_1001gplus/20260209_Exposito-Alonso/chr_only/TAIR10.chr.iupacN.fa \
  --TE_library "$T/combined_tair10_repbase.rm.lib" \
  --genotype false \
  --out "$OUT" \
  --repeatmasker_memory 100G --repeatmasker_time '48h' \
  -profile standard -work-dir "$T/nxf_work" -resume
echo "[$(date)] graffite done"
