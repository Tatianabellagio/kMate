#!/bin/bash
#SBATCH --job-name=kal_other_regimes
#SBATCH --partition=bse
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=4:00:00
#SBATCH --output=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/sims/visor_freqk/logs/kal_other_%j.out
#SBATCH --error=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/sims/visor_freqk/logs/kal_other_%j.err
# Run kallisto-EM (Route 3) on the 3 regimes where it's currently absent.
# Used to fill the plot matrix in RECOMB_SWEEP_RESULTS_slim — kallisto_em
# already has n50_g3; need n200_g1, n50_g1, n50_g3_skewed.
set -euo pipefail
PYBIN=/home/tbellagio/miniforge3/envs/hapfm/bin/python
ROOT=/carnegie/nobackup/scratch/tbellagio/hapfire_sv

cd $ROOT

declare -a SIMS=(
  "pool_sweep_82_recomb        cov50_n200_g1_s42_hotspots_p231_chr1"
  "pool_sweep_82_recomb        cov50_n50_g1_s42_hotspots_p231_chr1"
  "pool_sweep_82_recomb_skewed cov50_n50_g3_s42_hotspots_dom500_p231_chr1"
)

for entry in "${SIMS[@]}"; do
    SUBDIR=$(echo "$entry" | awk '{print $1}')
    SIM=$(echo "$entry" | awk '{print $2}')
    WORK=$ROOT/sims/visor_freqk/$SUBDIR/$SIM
    OUT=$WORK/kallisto_em.tsv
    if [[ -s $OUT ]]; then
        echo "[skip] $OUT exists"; continue
    fi
    echo "[$(date)] kallisto-EM on $SIM ..."
    $PYBIN -u poolfreq/src/per_sample_kallisto_em.py \
        --cn-kmer-prefix poolfreq/data/cn_full_231_v2/cn \
        --cn-var       poolfreq/data/cn_var_231_v2.cn_var.npz \
        --cn-var-meta  poolfreq/data/cn_var_231_v2.meta.npz \
        --reads $WORK/reads/r1.fq $WORK/reads/r2.fq \
        --sample ${SIM}_kallisto \
        --out $OUT \
        --chroms Chr1 --window-bp 10000 \
        --save-h-blocks $WORK/kallisto_em.h_blocks.npz \
        --global-anchor-weight 0.0
    echo "[$(date)] kallisto $SIM done"
done
echo "[$(date)] cross-regime kallisto sweep complete."
