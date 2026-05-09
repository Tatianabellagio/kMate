#!/bin/bash
# =============================================================================
# launch_v3_rebuild.sh
#
# Submit the v3 panel rebuild chain. Reads founders_231_chr.haploid.vcf.gz
# (produced by haploidize_merged_vcf.sh) and produces:
#
#   poolfreq/data/cn_var_231_v3.{cn_var,meta}.npz             (founder × variant)
#   sims/visor_freqk/founder_fastas_231_v3/<eco>.chr.fa × 231 (per-founder seq)
#   poolfreq/data/cn_full_231_v3/cn_Chr{1..5}.{cn,meta}.npz   (founder × k-mer)
#
# Dependency graph:
#
#   [haploid VCF, already on disk]
#     │
#     ├─► cn_var_v3            (single SLURM job, ~30 min)
#     │
#     └─► FASTAs_v3            (SLURM array 1..151, ~20-40 min wall, 16 concurrent)
#              │
#              └─► cn_full_v3  (SLURM array 1..5 chroms, ~2-4h per chrom)
#
# cn_var_v3 fires in parallel with FASTAs (independent paths). cn_full_v3
# waits for FASTAs to finish. The cactus 80 founders' FASTAs are pre-symlinked
# from the existing direct-assembly source (no rebuild needed for those).
#
# Pass --dependency=<jobid> to chain the whole pipeline behind another job
# (e.g. --dependency=61000 after the haploidize job, but 61000 already
# completed and the haploid VCF is on disk).
# =============================================================================
set -euo pipefail

BASE=/carnegie/nobackup/scratch/tbellagio/hapfire_sv
BCF=/home/tbellagio/miniforge3/envs/sequencing_pipeline/bin/bcftools

# ---- arg parsing -----------------------------------------------------------
DEP_JID=""
DRY_RUN=0
while [[ $# -gt 0 ]]; do
    case "$1" in
        --dependency)  DEP_JID="$2"; shift 2 ;;
        --dry-run)     DRY_RUN=1;     shift   ;;
        *)             echo "unknown arg: $1" >&2; exit 1 ;;
    esac
done
DEP_FLAG=""
[ -n "$DEP_JID" ] && DEP_FLAG="--dependency=afterok:$DEP_JID"

# ---- prerequisites --------------------------------------------------------
HAPLOID_VCF=$BASE/pangenie_genotyping/data/merged/founders_231_chr.haploid.vcf.gz
[ -s "$HAPLOID_VCF" ] || { echo "ERROR: missing $HAPLOID_VCF" >&2; exit 1; }

V3_AUX=$BASE/pangenie_genotyping/data/v3
mkdir -p $V3_AUX $BASE/poolfreq/logs

# 1. Lock down the canonical 231 founder order (used by every downstream tool).
FOUNDERS_FILE=$V3_AUX/founders_231_order.txt
$BCF query -l "$HAPLOID_VCF" > $FOUNDERS_FILE
N_FOUNDERS=$(wc -l < $FOUNDERS_FILE)
[ "$N_FOUNDERS" = "231" ] || { echo "ERROR: VCF has $N_FOUNDERS samples (expected 231)" >&2; exit 1; }

# 2. Make the PG-151 list explicit (= all 231 minus the cactus 80).
CACTUS_80=$BASE/pangenie_genotyping/data/merged/cactus_overlap_80.txt
PG_151=$V3_AUX/pangenie_151.txt
sort -u $CACTUS_80 > /tmp/cactus_80.tmp
sort -u $FOUNDERS_FILE > /tmp/all_231.tmp
comm -23 /tmp/all_231.tmp /tmp/cactus_80.tmp > $PG_151
N_PG=$(wc -l < $PG_151)
[ "$N_PG" = "151" ] || { echo "ERROR: derived $N_PG PanGenie founders (expected 151)" >&2; exit 1; }

# 3. Pre-symlink the 80 cactus-assembly FASTAs into founder_fastas_231_v3/.
# (PG-151 .chr.fa files will land here from the FASTAs SLURM array.)
FASTAS_V3_DIR=$BASE/sims/visor_freqk/founder_fastas_231_v3
mkdir -p $FASTAS_V3_DIR
N_LINKED=0
N_MISSING=0
while read eco; do
    src=/home/tbellagio/scratch/pang/pang_1001gplus/20260209_Exposito-Alonso/chr_only/${eco}.chr.fa
    # The existing v2 founder_fastas_231 dir already symlinks these — fall back
    # to chasing that symlink if the direct path doesn't exist (covers id-mapping
    # quirks in the cactus assembly catalog).
    if [ ! -s "$src" ]; then
        existing=$BASE/sims/visor_freqk/founder_fastas_231/${eco}.chr.fa
        if [ -L "$existing" ]; then
            src=$(readlink -f "$existing")
        fi
    fi
    if [ -s "$src" ]; then
        ln -sf "$src" "$FASTAS_V3_DIR/${eco}.chr.fa"
        # also link the .fai if present
        [ -s "${src}.fai" ] && ln -sf "${src}.fai" "$FASTAS_V3_DIR/${eco}.chr.fa.fai"
        N_LINKED=$((N_LINKED + 1))
    else
        echo "  WARN: no source FASTA for cactus founder $eco" >&2
        N_MISSING=$((N_MISSING + 1))
    fi
done < $CACTUS_80
echo "  symlinked $N_LINKED of 80 cactus founders into $FASTAS_V3_DIR (missing: $N_MISSING)"

# ---- summary --------------------------------------------------------------
echo ""
echo "== v3 rebuild plan =="
echo "  haploid VCF:    $HAPLOID_VCF"
echo "  founders order: $FOUNDERS_FILE  (n=$N_FOUNDERS)"
echo "  PG-151 list:    $PG_151"
echo "  FASTAs out:     $FASTAS_V3_DIR (cactus 80 already symlinked, PG-151 will populate)"
echo "  cn_var out:     $BASE/poolfreq/data/cn_var_231_v3.{cn_var,meta}.npz"
echo "  cn_full out:    $BASE/poolfreq/data/cn_full_231_v3/cn_Chr{1..5}.{cn,meta}.npz"
[ -n "$DEP_JID" ] && echo "  dep on:         job $DEP_JID (afterok)"
[ $DRY_RUN -eq 1 ] && { echo ""; echo "(dry-run — no jobs submitted)"; exit 0; }

# ---- submit ---------------------------------------------------------------
echo ""
echo "== submitting =="

# (a) cn_var_v3 — depends only on haploid VCF (already on disk; chain to user-
#     provided dep job if any).
JID_CNVAR=$(sbatch --parsable $DEP_FLAG \
    $BASE/poolfreq/scripts/build_cn_var_v3.sh)
echo "  cn_var_v3:    job $JID_CNVAR"

# (b) FASTAs array — 1..151 PG founders, 16 concurrent.
JID_FASTAS=$(sbatch --parsable $DEP_FLAG \
    --array=1-${N_PG}%16 \
    $BASE/pangenie_genotyping/scripts/build_consensus_fastas_one.sh)
echo "  FASTAs array: job $JID_FASTAS  (${N_PG} tasks, 16 concurrent)"

# (c) cn_full_v3 array — 1..5 chroms, all concurrent. Depends on FASTAs array.
JID_CNFULL=$(sbatch --parsable \
    --dependency=afterok:$JID_FASTAS \
    --array=1-5 \
    $BASE/poolfreq/scripts/build_cn_full_v3_one.sh)
echo "  cn_full_v3:   job $JID_CNFULL  (5 chroms, depends on $JID_FASTAS)"

echo ""
echo "== submitted =="
echo "  squeue -u tbellagio --format='%.12i %.20j %.10T %.6M %R'"
echo ""
echo "== expected wall ==
  cn_var_v3:    ~30 min   (after haploid VCF; in parallel with FASTAs)
  FASTAs:       ~20-40 min (151 / 16 concurrent ≈ 10 batches × ~5-10 min each)
  cn_full_v3:   ~2-4 h per chrom, 5 chroms in parallel → ~2-4 h
  total:        ~3-5 h after FASTAs land
"
