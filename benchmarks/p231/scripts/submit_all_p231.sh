#!/bin/bash
# =============================================================================
# benchmarks/p231 -- full 231-panel benchmark DAG (paper headline).
# Mirrors benchmarks/p80 but: (a) RANDOM crossovers @ 4 cM/Mb (no BigLD forcing),
# (b) truth + projection on BOTH arch3 var_pas (atomized + raw), (c) kmer_pa
# rebuilt from the arch3 canonical VCF for single-source provenance.
#
# Reused (verified current v3qc+arch, see README): arch3 var_pas, merged_231 VCF.
# Rebuilt here: fastas_231 (A5), kmer_pa_p231[_filt2|_filt2inv] (A3/A3b/A3c), sims (B), EM (C).
#
# Usage: bash scripts/submit_all_p231.sh    (idempotent: stages skip if present)
# =============================================================================
set -euo pipefail
cd /global/scratch/users/tbellg/kmate/benchmarks/p231
S=scripts

A5=$(sbatch --parsable $S/05_build_fastas_p231.sh)                          # 231 FASTAs from arch3 VCF
A3=$(sbatch --parsable $S/03_build_kmer_pa_p231.sh)                         # kmer_pa from arch3 VCF
A3B=$(sbatch --dependency=afterok:$A3 --parsable $S/03b_build_kmer_pa_filt2_p231.sh)
A3C=$(sbatch --dependency=afterok:$A3 --parsable $S/03c_build_kmer_pa_filt2inv_p231.sh)   # 07f's kmer_pa

# Phase B: sims (need FASTAs). Random crossovers; truth on both var_pas.
B_n50g0=$(sbatch  --dependency=afterok:$A5 --parsable $S/06_run_sim_p231.sh 50  0)
B_n231g0=$(sbatch --dependency=afterok:$A5 --parsable $S/06_run_sim_p231.sh 231 0)
B_n50g1=$(sbatch  --dependency=afterok:$A5 --parsable $S/06_run_sim_p231.sh 50  1)
B_n231g1=$(sbatch --dependency=afterok:$A5 --parsable $S/06_run_sim_p231.sh 231 1)
B_n50g3=$(sbatch  --dependency=afterok:$A5 --parsable $S/06_run_sim_p231.sh 50  3)
B_dom=$(sbatch    --dependency=afterok:$A5 --parsable $S/06b_run_sim_p231_skewed.sh 50 3 42 50.0)

# Phase C: front-runner EM (07f: filt2inv kmer_pa + per_founder normalization,
# GLOBAL mode), both var_pa arms, per regime.
declare -A BJOB=( [n50_g0]=$B_n50g0 [n231_g0]=$B_n231g0 [n50_g1]=$B_n50g1
                  [n231_g1]=$B_n231g1 [n50_g3]=$B_n50g3 [n50_g3_dom500]=$B_dom )
CJOBS=""
for reg in n50_g0 n231_g0 n50_g1 n231_g1 n50_g3 n50_g3_dom500; do
  for cnvar in atomized raw; do
    cj=$(sbatch --dependency=afterok:${BJOB[$reg]}:$A3C --parsable \
         $S/07f_run_kmate_filt2inv_p231.sh $reg $cnvar)
    CJOBS="$CJOBS:$cj"
  done
done

# Phase D: score (after all EM).
sbatch --dependency=afterok${CJOBS} --wrap \
  "/global/home/users/tbellg/miniforge3/envs/hapfm/bin/python $PWD/$S/score_p231.py" \
  --job-name=p231_score --account=co_moilab --partition=savio4_htc --qos=moilab_htc4_normal --mem=32G --time=1:00:00 \
  -o logs/score_%j.out -e logs/score_%j.err

echo "submitted: A5=$A5 A3=$A3 A3B=$A3B A3C=$A3C sims=[$B_n50g0 $B_n231g0 $B_n50g1 $B_n231g1 $B_n50g3 $B_dom]"
echo "EM jobs:$CJOBS"
