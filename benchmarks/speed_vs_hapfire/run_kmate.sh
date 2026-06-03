#!/bin/bash
#SBATCH --job-name=bench_kmate
#SBATCH --account=fc_moilab
#SBATCH --partition=savio3
#SBATCH --qos=savio_normal
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --output=benchmarks/speed_vs_hapfire/logs/kmate_%j.out
#SBATCH --error=benchmarks/speed_vs_hapfire/logs/kmate_%j.err
set -eo pipefail

# Speed benchmark: kMate on the SAME 231-founder Chr1 pool hapFIRE is timed on
# (visor_freqk del/rep1/cov50/var_del_1kb_n231_f30). Per-sample driver:
# k-mer count (jellyfish) + Poisson EM on the 231-founder simplex + var_pa
# projection -> per-record alt_freq. Production recipe (arch3 matrices, inv_mb,
# global mode). 8 threads, matched to hapFIRE's allocation.

ROOT=/global/scratch/users/tbellg/kmate
RES=$ROOT/benchmarks/speed_vs_hapfire/results
PY=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python
SIMDIR=/global/scratch/users/tbellg/visor_freqk/data/reads_var/del/rep1/cov50/var_del_1kb_n231_f30_err001

mkdir -p "$RES"; cd "$RES"

# Pin all BLAS/threadpool backends to 8 so kMate is genuinely capped at 8 cores
# (savio3 nodes have 32; without this numpy/OpenBLAS grabs more during EM).
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8 \
       NUMEXPR_NUM_THREADS=8 VECLIB_MAXIMUM_THREADS=8

echo "[$(date)] kMate benchmark start on $(hostname), ${SLURM_CPUS_PER_TASK:-?} cpus (threads pinned to 8)"

/usr/bin/time -v -o "$RES/kmate_time_pinned8.txt" $PY "$ROOT/src/per_sample_per_chrom.py" \
    --kmer-pa-prefix "$ROOT/data/kmer_pa_231_arch3_filt2inv/kmer_pa" \
    --var-pa     "$ROOT/panel/arch3/chr1/var_pa_231_arch3_chr1.var_pa.npz" \
    --var-called "$ROOT/panel/arch3/chr1/var_pa_231_arch3_chr1.var_called.npz" \
    --var-meta   "$ROOT/panel/arch3/chr1/var_pa_231_arch3_chr1.meta.npz" \
    --reads "$SIMDIR/r1.fq" "$SIMDIR/r2.fq" \
    --sample kmate_n231_chr1_cov50_f30_pin8 \
    --out kmate_n231_chr1_cov50_f30_pin8.tsv \
    --threads 8 --chroms Chr1 --kmer-weight inv_mb --block-mode global

echo "[$(date)] kMate benchmark done"
ls -la "$RES"
