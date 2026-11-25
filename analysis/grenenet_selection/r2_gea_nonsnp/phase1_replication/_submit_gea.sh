#!/bin/bash
# The production GEA end to end, as one chain of Slurm jobs (README.md in this folder: what each step does and writes).
# Input made beforehand: the generation-9 pool matrices (assemble_lastgen.py --kind snp / --kind nonsnp).
#   1. class matrices, call-rate floor 0.50             run_class_matrices_cr50.sbatch
#   2. record-arm LFMM, snp + smallindel x 22 axes      multiaxis/ma_lfmm_cr50.sbatch            (after 1)
#   3. event-arm inputs (Truvari events)                genes/convergence/event_gea/run_event_matrix_cr50.sbatch
#   4. event-arm LFMM, sv + pooled non-SNP x 22 axes    genes/convergence/event_gea/run_event_lfmm_cr50.sbatch (after 3)
#   5. one table per class x axis, tiling blocks        multiaxis/run_reblock_cr50.sbatch        (after 2 and 4)
#   6. significant blocks per class                     genes/convergence/class_signal_counts.py
#      region / hit / content tables                    annotate_gea_universe.py, gea_hit_landscape_data.py, gea_sv_cargo_data.py
#      every GEA Results number                         genes/convergence/gea_results_numbers.py
#      BH / GIF / genome-fraction diagnostic            multiaxis/raw_block_significance.py   (all of step 6 after 5)
# Then, by hand on a compute node: notebooks/gea_draft.ipynb and genes/convergence/_build_gea_hit_landscape_nb.py.
# Usage (any node): [KMATE_PANEL=arch3v2] [FROM=n] [TO=n] bash _submit_gea.sh   -> job ids in results/_submit_gea[_arch3v2]_jobids.txt
# Steps 1-5 rewrite the production inputs (hours, ~180 GB per LFMM task); `FROM=6 bash _submit_gea.sh` reruns only step 6.
set -euo pipefail
GS=/global/scratch/projects/fc_moilab/tbellg/kmate/analysis/grenenet_selection
P=$GS/r2_gea_nonsnp/phase1_replication; EV=$GS/genes/convergence/event_gea; CONV=$GS/genes/convergence
K=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python
FROM=${FROM:-1}; TO=${TO:-6}
[[ "$FROM" =~ ^[1-6]$ && "$TO" =~ ^[1-6]$ ]] || { echo "FROM / TO must be 1-6, got '$FROM' / '$TO'" >&2; exit 1; }
# panel version (lib.PANEL, env KMATE_PANEL, exported to every job): every output directory gets the suffix $S, and the
# logs go to logs/gea$S/ (the #SBATCH --output lines name the arch3 folders, so they are overridden here)
V=${KMATE_PANEL:-arch3v2}; S=""; [ "$V" = arch3 ] || S="_$V"; export KMATE_PANEL=$V
LG=/global/scratch/projects/fc_moilab/tbellg/kmate/logs/gea$S; mkdir -p $LG
IDS=$P/results/_submit_gea${S}_jobids.txt; echo "# $(date -Is) KMATE_PANEL=$V FROM=$FROM TO=$TO" >> $IDS; log() { echo "$1 $2" | tee -a $IDS; }
SB="sbatch --parsable --export=ALL"
dep() { local d=""; for j in "$@"; do [ -n "$j" ] && d="$d:$j"; done; if [ -n "$d" ]; then echo "--dependency=afterok$d"; fi; }

C=""; L=""; M=""; E=""; R=""
run() { [ "$FROM" -le "$1" ] && [ "$TO" -ge "$1" ]; }
# step 3 checks its short indels against step 1's class matrices, so it follows step 1 when both run
if run 1; then C=$($SB -o $LG/class_matrices_%j.out -e $LG/class_matrices_%j.out $P/run_class_matrices_cr50.sbatch); log class_matrices $C; fi
if run 2; then L=$($SB $(dep $C) -o $LG/lfmm_records_%A_%a.out -e $LG/lfmm_records_%A_%a.out $P/multiaxis/ma_lfmm_cr50.sbatch); log lfmm_records $L; fi
if run 3; then M=$($SB $(dep $C) -o $LG/event_matrix_%j.out -e $LG/event_matrix_%j.out $EV/run_event_matrix_cr50.sbatch); log event_matrix $M; fi
if run 4; then E=$($SB $(dep $M) -o $LG/lfmm_events_%A_%a.out -e $LG/lfmm_events_%A_%a.out $EV/run_event_lfmm_cr50.sbatch); log lfmm_events $E; fi
if run 5; then R=$($SB $(dep $L $E) -o $LG/reblock_%A_%a.out -e $LG/reblock_%A_%a.out $P/multiaxis/run_reblock_cr50.sbatch); log reblock $R; fi
run 6 || exit 0
N=$($SB $(dep $R) --account=co_moilab --partition=savio4_htc --qos=savio_lowprio --cpus-per-task=16 --mem=96G \
    --time=6:00:00 --job-name=gea_numbers -o $LG/gea_numbers_%j.out -e $LG/gea_numbers_%j.out --wrap="set -e; cd $CONV; \
    $K class_signal_counts.py && $K annotate_gea_universe.py 16 && $K gea_hit_landscape_data.py && \
    $K gea_sv_cargo_data.py && $K gea_results_numbers.py && $K $P/multiaxis/raw_block_significance.py --model lfmm"); log numbers $N
