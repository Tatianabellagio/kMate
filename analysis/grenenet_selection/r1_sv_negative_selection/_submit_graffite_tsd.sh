#!/bin/bash
# GraffiTE TSD search ONLY, re-run with the flank-naming fix (graffite_prepTSD_flank_names.patch).
#
# Upstream commit 3f062c1 made prepTSD.sh name the reference flanks by region, while TSD_Match_v2.sh looks
# them up by SV ID, so the original run (_submit_graffite.sh, job 39593195) searched every TSD without its
# flanking sequence (395 PASS of 42,837). The patch restores the pre-3f062c1 bedtools route. On 200 Chr2
# records it took PASS from 3 to 114, with superfamily-typical lengths (9 bp for MuDR, longer and variable
# for LINEs).
#
# Enters at GraffiTE's documented --RM_dir point: the published 2_Repeat_Filtering/{1..5} of the original
# run (genotypes_repmasked_filtered.vcf + repeatmasker_dir), so RepeatMasker is not repeated. Same
# GraffiTE checkout, container and reference. Fresh work and output directories; the original run is left
# untouched for comparison.
#
# Scope is GraffiTE's own: records in pangenome.vcf's input (OneCode span > 0.80) with a single TE copy
# (n_hits = 1). That excludes the replacement records (never masked by GraffiTE) and Helitrons (dropped
# by OneCode; they make no TSD anyway).
#SBATCH --job-name=gt_tsd
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=32 --mem=64G
#SBATCH --time=24:00:00
#SBATCH --exclude=n0041.savio4
#SBATCH --output=logs/graffite_tsd_%A.out
set -euo pipefail
# absolute, not $SLURM_SUBMIT_DIR: from inside an salloc shell that points at wherever salloc was run
cd /global/scratch/projects/fc_moilab/tbellg/kmate/analysis/grenenet_selection/r1_sv_negative_selection
module load nextflow/24.10.4
T=/global/scratch/projects/fc_moilab/tbellg/kmate/telib
grep -q -- '-nameOnly' "$T/graffite_src/bin/prepTSD.sh" || {
  echo "prepTSD.sh is not patched: git -C $T/graffite_src apply $PWD/graffite_prepTSD_flank_names.patch"; exit 1; }
export NXF_SINGULARITY_CACHEDIR=$T/nxf_cache
export APPTAINER_CACHEDIR=$T/.apptainer_cache
export NXF_OPTS="-Xms2g -Xmx8g"
OUT=$PWD/results/sv_adaptive/graffite/run_tsd
mkdir -p "$OUT" "$T/nxf_work_tsd"
nextflow run $T/graffite_src/main.nf \
  --RM_dir "$PWD/results/sv_adaptive/graffite/run/2_Repeat_Filtering" \
  --reference /global/scratch/users/tbellg/pang/pang_1001gplus/20260209_Exposito-Alonso/chr_only/TAIR10.chr.iupacN.fa \
  --genotype false \
  --out "$OUT" \
  -profile standard -work-dir "$T/nxf_work_tsd"
echo "[$(date)] graffite TSD rerun done: $(awk '$NF=="PASS"' "$OUT/3_TSD_search/TSD_summary.txt" | wc -l) PASS of $(wc -l < "$OUT/3_TSD_search/TSD_summary.txt")"
