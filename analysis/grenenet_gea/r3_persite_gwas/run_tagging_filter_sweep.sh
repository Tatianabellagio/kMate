#!/bin/bash
set -uo pipefail
cd /global/scratch/users/tbellg/kmate/analysis/grenenet_gea
PY=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python3
mkdir -p sv_snp_ld_v2_nofilter sv_snp_ld_v2_maconly logs

export OMP_NUM_THREADS=2
export OPENBLAS_NUM_THREADS=2
export MKL_NUM_THREADS=2

run_one() {
  local mac="$1" minpair="$2" outdir="$3" cls="$4" mode="$5" chrom="$6"
  local outfile="${outdir}/tagging_${cls}_${mode}_${chrom}.npz"
  if [[ -s "$outfile" ]]; then
    echo "[skip] $outfile already exists"
    return 0
  fi
  echo "[start] $outfile $(date)"
  "$PY" build_tagging_masked.py --chrom "$chrom" --anchor-class "$cls" --panel-mode "$mode" \
      --mac-floor "$mac" --min-pair-frac "$minpair" --out-dir "$outdir" \
      > "logs/tagging_$(basename "$outdir")_${cls}_${mode}_${chrom}.log" 2>&1
  echo "[done]  $outfile $(date)"
}
export -f run_one

MAXJOBS=2
njobs=0
for cls in sv indel; do
  for mode in panel shortread; do
    for chrom in Chr1 Chr2 Chr3 Chr4 Chr5; do
      run_one 1 0.0 sv_snp_ld_v2_nofilter "$cls" "$mode" "$chrom" &
      njobs=$((njobs+1))
      if (( njobs % MAXJOBS == 0 )); then wait -n; fi
      run_one 2 0.0 sv_snp_ld_v2_maconly "$cls" "$mode" "$chrom" &
      njobs=$((njobs+1))
      if (( njobs % MAXJOBS == 0 )); then wait -n; fi
    done
  done
done
wait
echo "ALL_TAGGING_SWEEP_JOBS_DONE"
