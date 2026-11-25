#!/bin/bash
# GraffiTE's masking step on the 71,865 SV records GraffiTE skipped (COMPLEX / base-swap records that
# `bcftools view --types indels` rejects; see _graffite_prep_skipped.py). One array task per chromosome.
#
# Mirrors GraffiTE bin/repmask_vcf.sh exactly -- same container (graffite.sif is byte-identical to the
# image the Nextflow run used), same library, RepeatMasker -s, build_dictionary.pl, OneCode, same parse of
# the OneCode output -- so the skipped records are annotated on the same footing as the masked ones.
# What is deliberately NOT run: annotate_vcf.R / the span filter / the TSD search. The parse in
# _graffite_classify.py computes spans from the raw .out for every record, and a TSD is undefined for a
# replacement (no single insertion point).
#
# --home on scratch: RepeatMasker writes a blast cache to $HOME, and the home filesystem's locking makes
# makeblastdb fail inside apptainer (_submit_rm_cargo_libs.sh).
# TEST=N masks only the first N sequences of the task's FASTA (or of TEST_FA), into a separate test/
# directory.
#SBATCH --job-name=gt_skip
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=32 --mem=64G
#SBATCH --time=72:00:00
#SBATCH --array=1-5
#SBATCH --output=logs/gt_skip_%A_%a.out
set -euo pipefail
# absolute, not $SLURM_SUBMIT_DIR: from inside an salloc shell that points at wherever salloc was run
cd /global/scratch/projects/fc_moilab/tbellg/kmate/analysis/grenenet_selection/r1_sv_negative_selection
T=/global/scratch/projects/fc_moilab/tbellg/kmate/telib
BIN=$T/graffite_src/bin
LIB=$T/combined_tair10_repbase.rm.lib
CH=Chr${SLURM_ARRAY_TASK_ID:?}
SRC=$PWD/results/sv_adaptive/graffite/skipped/$CH
D=$SRC
if [ -n "${TEST:-}" ]; then
  D=$SRC/test; rm -rf "$D"; mkdir -p "$D"
  head -n $((2 * TEST)) "${TEST_FA:-$SRC/indels.fa}" > "$D/indels.fa"   # TEST_FA: e.g. GraffiTE's own query, to check parity
fi
mkdir -p "$T/rmhome_skip$CH"
cd "$D"
rm -rf repeatmasker_dir
mkdir repeatmasker_dir
NP=${SLURM_CPUS_PER_TASK:-$(nproc)}

apptainer exec -B /global/scratch:/global/scratch --home "$T/rmhome_skip$CH:/home/rmuser" "$T/graffite.sif" \
  bash -euo pipefail -c "
    cd '$D'
    REPMASK_DIR=repeatmasker_dir
    REPMASK_OUT=\$REPMASK_DIR/indels.fa.out
    REPMASK_ONECODE_OUT=\$REPMASK_DIR/indels.fa.onecode.out
    /opt/RepeatMasker/RepeatMasker -lib '$LIB' -s -dir \$REPMASK_DIR -pa $NP indels.fa
    echo 'Building onecode LTR dictionary...'
    perl '$BIN/build_dictionary.pl' --rm \$REPMASK_OUT --unknown > \$REPMASK_DIR/OneCode_LTR.dic
    echo 'Running onecode...'
    perl '$BIN/one_code_to_find_them_all.pl' --rm \$REPMASK_OUT --unknown --ltr \$REPMASK_DIR/OneCode_LTR.dic > \$REPMASK_DIR/onecode.log 2>&1
    rm -f \$REPMASK_DIR/ALL.elem_sorted.csv
    for i in \$REPMASK_DIR/*.elem_sorted.csv; do cat \"\$i\" >> \$REPMASK_DIR/ALL.elem_sorted.csv; done
    cat <(head -n 3 \$REPMASK_OUT | sed 's/ID/ID\tfragmts/g') <(grep '###' \$REPMASK_DIR/ALL.elem_sorted.csv | cut -f 1-16 | sed 's/###//g') > \$REPMASK_ONECODE_OUT
    mv \$REPMASK_DIR/ALL.elem_sorted.csv \$REPMASK_DIR/ALL.onecode.elem_sorted.bak
    find \$REPMASK_DIR -name '*.copynumber.csv' -delete
    find \$REPMASK_DIR -name '*.elem_sorted.csv' -delete
    find \$REPMASK_DIR -name '*.ltr.csv' -delete
    find \$REPMASK_DIR -name '*.transposons.csv' -delete
    samtools faidx indels.fa
    awk '{print \$1\"\t\"\$2}' indels.fa.fai > indels.length
  "
echo "[$(date)] $CH done: $(grep -c '^>' indels.fa) sequences, $(awk 'NR>3' repeatmasker_dir/indels.fa.out | wc -l) RM hits, $(awk 'NR>3' repeatmasker_dir/indels.fa.onecode.out | wc -l) OneCode copies"
