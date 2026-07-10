#!/bin/bash
#SBATCH --job-name=ldblk_win
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=72G
#SBATCH --time=06:00:00
#SBATCH --output=benchmarks/ldblock_window_test/logs/%x_%j.out
#SBATCH --error=benchmarks/ldblock_window_test/logs/%x_%j.out
set -eo pipefail
# window-mode kMate using LD blocks as windows.
# Usage: sbatch run_ldblock_window.sbatch <POOL> <PANEL> <BLOCKS_TSV> [MIN_KMERS=50] [WEIGHT=inv_mb]
cd /global/scratch/users/tbellg/kmate
mkdir -p benchmarks/ldblock_window_test/logs
POOL=${1:?pool}; PANEL=${2:?p80|p231}; BLOCKS=${3:?blocks tsv}
MKB=${4:-50}; WEIGHT=${5:-inv_mb}
T=${SLURM_CPUS_PER_TASK:-8}
source ~/miniforge3/etc/profile.d/conda.sh && conda activate kmate
case "$PANEL" in
  p231) KP=data/kmer_pa_231_arch3_filt2inv/kmer_pa
        VP=panel/arch3/chr1/var_pa_231_arch3_chr1.var_pa.npz
        VM=panel/arch3/chr1/var_pa_231_arch3_chr1.meta.npz ;;
  p80)  KP=benchmarks/p80/data/kmer_pa_p80_filt2/kmer_pa
        VP=benchmarks/p80/data/var_pa_p80.var_pa.npz
        VM=benchmarks/p80/data/var_pa_p80.meta.npz ;;
  *) echo "bad panel $PANEL"; exit 1 ;;
esac
BTAG=$(basename "$BLOCKS" .tsv)          # encode the block map in the output name
OUT=benchmarks/ldblock_window_test/${POOL}_${PANEL}_${BTAG}_mkb${MKB}.tsv
echo "[$(date)] LD-block window kMate: pool=$POOL panel=$PANEL mkb=$MKB on $(hostname)"
echo "  kmer_pa=$KP  var_pa=$VP  blocks=$BLOCKS"
python -m kmate.per_sample_per_chrom \
  --kmer-pa-prefix "$KP" --var-pa "$VP" --var-meta "$VM" \
  --reads benchmarks/$PANEL/sims/$POOL/reads/r1.fq benchmarks/$PANEL/sims/$POOL/reads/r2.fq \
  --sample ldblk_${POOL} --threads "$T" --chroms Chr1 \
  --block-mode window --kmer-weight "$WEIGHT" \
  --blocks-tsv "$BLOCKS" --min-kmers-per-block "$MKB" \
  --out "$OUT"
echo "[$(date)] DONE -> $OUT"
