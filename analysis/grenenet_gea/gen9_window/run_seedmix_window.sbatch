#!/bin/bash
#SBATCH --job-name=sm_window
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_lowprio
#SBATCH --cpus-per-task=8
#SBATCH --mem=72G
#SBATCH --time=06:00:00
#SBATCH --array=1-8
#SBATCH --requeue
#SBATCH --output=logs/sm_window_%A_%a.out
#SBATCH --error=logs/sm_window_%A_%a.err
# -----------------------------------------------------------------------------
# Window-mode kMate on the 8 SEEDMIX (gen-0 founding) reps, on the SAME dynld-K500
# block map + params as the cohort window run (logs/kmate_grenenet_*; min-kmers 50,
# inv_mb, var_called passed). Produces SEEDMIX_S{N}_Chr{M}.h_blocks_per_chrom.npz so
# the founding haplotype-frequency p0 projects onto the SAME units as the cohort
# hapfreq matrix (results/grenenet_gea/hapfreq/).
#
# SEEDMIX is PCR-free (Lucigen) -> use trim-only reads (_P.fq.gz), NOT dedup.
# Count-once: build the Jellyfish DB once per sample in node-local /dev/shm, then
# query each chrom (the cohort's scale-out I/O fix).
#
#   sbatch run_seedmix_window.sbatch        # all 8 reps (array task = sample idx)
# -----------------------------------------------------------------------------
set -euo pipefail
hostname
cd /global/scratch/users/tbellg/kmate
mkdir -p logs results/grenenet_kmate_window_seedmix

source ~/miniforge3/etc/profile.d/conda.sh && conda activate kmate

S=S${SLURM_ARRAY_TASK_ID}                                  # S1..S8
READS_DIR=/global/scratch/users/tbellg/pang/grenenet_reads/seed_mix
R1=$READS_DIR/${S}-1.1_P.fq.gz
R2=$READS_DIR/${S}-1.2_P.fq.gz
KP=data/kmer_pa_231_arch3_filt2inv/kmer_pa
OUTDIR=results/grenenet_kmate_window_seedmix
T=${SLURM_CPUS_PER_TASK:-8}

for f in "$R1" "$R2"; do [ -s "$f" ] || { echo "ERROR missing $f" >&2; exit 1; }; done

# --- count once: build k-mer DB in node-local RAM -----------------------------
DB=/dev/shm/SEEDMIX_${S}.${SLURM_JOB_ID:-x}.jf
trap 'rm -f "$DB"' EXIT
echo "[$(date)] build-kmer-db $S -> $DB"
kmate build-kmer-db --reads "$R1" "$R2" --out "$DB" --threads "$T" --hash-size 3G

# --- query each chrom in window mode on its dynld-K500 unit map ----------------
for N in 1 2 3 4 5; do
  VP=panel/arch3/chr${N}/var_pa_231_arch3_chr${N}.var_pa.npz
  VC=panel/arch3/chr${N}/var_pa_231_arch3_chr${N}.var_called.npz
  VM=panel/arch3/chr${N}/var_pa_231_arch3_chr${N}.meta.npz
  BLK=results/grenenet_gea/blocks_mcf90/chr${N}_units_dynld_K500.tsv
  OUT=$OUTDIR/SEEDMIX_${S}_Chr${N}.tsv
  echo "[$(date)] $S Chr${N} -> $OUT"
  kmate run \
    --kmer-pa-prefix "$KP" --var-pa "$VP" --var-called "$VC" --var-meta "$VM" \
    --reads "$R1" "$R2" --kmer-db "$DB" \
    --sample "SEEDMIX_${S}" --out "$OUT" \
    --threads "$T" --chroms Chr${N} \
    --block-mode window --kmer-weight inv_mb \
    --blocks-tsv "$BLK" --min-kmers-per-block 50
done
echo "[$(date)] DONE $S"
