#!/bin/bash
#SBATCH --job-name=ldblk_coarse
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=72G
#SBATCH --time=06:00:00
#SBATCH --array=0-29
#SBATCH --output=benchmarks/ldblock_window_test/logs/%x_%A_%a.out
#SBATCH --error=benchmarks/ldblock_window_test/logs/%x_%A_%a.out
# Coarseness benchmark: window-mode kMate on base + floor{8,15,25,40} block maps,
# across 6 replicate seeds (read/sim randomness). p231, Chr1, recombinant g3 self97 pool.
# task = seed_idx*5 + level_idx ; seeds s42..s47 ; levels {base,8,15,25,40}.
set -eo pipefail
cd /global/scratch/users/tbellg/kmate
mkdir -p benchmarks/ldblock_window_test/logs
SEEDS=(42 43 44 45 46 47); LEVELS=(base 8 15 25 40)
SI=$(( SLURM_ARRAY_TASK_ID / 5 )); LI=$(( SLURM_ARRAY_TASK_ID % 5 ))
SEED=${SEEDS[$SI]}; LEVEL=${LEVELS[$LI]}
POOL=cov10_n50_g3_s${SEED}_self97_hotspots_dom500_p231_chr1
MKB=50; WEIGHT=inv_mb; T=${SLURM_CPUS_PER_TASK:-8}
source ~/miniforge3/etc/profile.d/conda.sh && conda activate kmate
KP=data/kmer_pa_231_arch3_filt2inv/kmer_pa
VP=panel/arch3/chr1/var_pa_231_arch3_chr1.var_pa.npz
VM=panel/arch3/chr1/var_pa_231_arch3_chr1.meta.npz
if [ "$LEVEL" = "base" ]; then
  BLOCKS=analysis/grenenet_gea/blocks_mcf90/chr1_clq0.9_blocks_clq0.9.tsv
else
  BLOCKS=analysis/grenenet_gea/blocks_mcf90/coarse/chr1_floor${LEVEL}.tsv
fi
OUT=benchmarks/ldblock_window_test/sweep_s${SEED}_${LEVEL}.tsv
if [ ! -d "benchmarks/p231/sims/$POOL/reads" ]; then echo "MISSING pool $POOL"; exit 0; fi
echo "[$(date)] seed=$SEED level=$LEVEL pool=$POOL blocks=$BLOCKS on $(hostname)"
python -m kmate.per_sample_per_chrom \
  --kmer-pa-prefix "$KP" --var-pa "$VP" --var-meta "$VM" \
  --reads benchmarks/p231/sims/$POOL/reads/r1.fq benchmarks/p231/sims/$POOL/reads/r2.fq \
  --sample sweep_s${SEED}_${LEVEL} --threads "$T" --chroms Chr1 \
  --block-mode window --kmer-weight "$WEIGHT" \
  --blocks-tsv "$BLOCKS" --min-kmers-per-block "$MKB" \
  --out "$OUT"
echo "[$(date)] DONE seed=$SEED level=$LEVEL -> $OUT"
