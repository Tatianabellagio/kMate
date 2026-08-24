#!/bin/bash
# Pull SLURM accounting data for all hapfire_1141 jobs and write a stats TSV.
#
# Usage:
#   bash collect_stats.sh [STATS_FILE]
#     STATS_FILE defaults to results/job_stats.tsv
#
# Run anytime - it overwrites the stats file with the current sacct snapshot.
# Re-run after fan-out completes to get final MaxRSS for every job.

set -eo pipefail

ROOT=/global/scratch/users/tbellg/hapfire_sv/contamination_test
OUT="${1:-${ROOT}/results/job_stats.tsv}"
mkdir -p "$(dirname "${OUT}")"

# Find all SLURM jobs whose name starts with "hapfire_" launched by us.
# MaxRSS lives on the .batch step (not the job-level row), so we pull it
# separately and join.
SACCT_RAW=$(sacct --user="$USER" --partition=savio3_xlmem \
                  --starttime=2026-05-05 \
                  --format=JobID,JobName%50,State,Elapsed,MaxRSS,NodeList%30,Start,End,AllocCPUs,ReqMem \
                  -P -n 2>/dev/null \
            | awk -F'|' '$2 ~ /^(hapfire|batch|extern)/')

# build jobid -> MaxRSS map (from .batch step rows)
declare -A RSS_BY_JOB
while IFS='|' read jid _ _ _ rss _ _ _ _ _; do
  if [[ "$jid" == *.batch ]]; then
    RSS_BY_JOB[${jid%.batch}]="$rss"
  fi
done <<<"$SACCT_RAW"

{
  echo -e "jobid\tjob_name\tsample\tchrom\tstate\telapsed\telapsed_sec\tmax_rss\tmax_rss_gb\tnode\tstart\tend\tcpus\tmem_alloc"
  echo "$SACCT_RAW" \
  | awk -F'|' '$1 !~ /\.(batch|extern|[0-9]+)$/ {print}' \
  | while IFS='|' read jobid jobname state elapsed _ node start end cpus mem; do
      maxrss="${RSS_BY_JOB[$jobid]:-}"
      [[ -z "$jobid" ]] && continue
      # parse hapfire_<sample>_chr<N> -> sample, chrom
      if [[ "$jobname" =~ ^hapfire_(SEEDMIX_S[0-9]+)_chr([0-9]+)$ ]]; then
        sample="${BASH_REMATCH[1]}"
        chrom="${BASH_REMATCH[2]}"
      else
        # fall back: try to recover from the job's stdout log
        sample=""; chrom=""
        for ext in "_4294967294.out" ".out"; do
          for log in "${ROOT}/logs/hapfire_${jobname}_${jobid}${ext}" "${ROOT}/logs/hapfire_hapfire_1141_${jobid}${ext}"; do
            [[ -f "$log" ]] || continue
            line=$(grep -m1 "hapFIRE start: sample=" "$log" 2>/dev/null)
            if [[ -n "$line" ]]; then
              sample=$(echo "$line" | sed -nE 's/.*sample=([A-Z_0-9]+) .*/\1/p')
              chrom=$(echo "$line" | sed -nE 's/.* chr=([0-9]+).*/\1/p')
              break 2
            fi
          done
        done
      fi
      # convert elapsed HH:MM:SS or D-HH:MM:SS -> seconds
      esec=$(echo "$elapsed" | awk -F'[-:]' '{
        if (NF==4) print $1*86400 + $2*3600 + $3*60 + $4
        else if (NF==3) print $1*3600 + $2*60 + $3
        else print 0
      }')
      # convert MaxRSS K/M/G to GB float
      rss_gb=$(echo "$maxrss" | awk '{
        if ($0 ~ /K$/) printf "%.2f", substr($0,1,length($0)-1)/1024/1024
        else if ($0 ~ /M$/) printf "%.2f", substr($0,1,length($0)-1)/1024
        else if ($0 ~ /G$/) printf "%.2f", substr($0,1,length($0)-1)
        else if ($0 ~ /T$/) printf "%.2f", substr($0,1,length($0)-1)*1024
        else if ($0+0 > 0) printf "%.2f", $0/1024/1024  # bare number assumed K
        else print ""
      }')
      printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "$jobid" "$jobname" "$sample" "$chrom" "$state" "$elapsed" "$esec" "$maxrss" "$rss_gb" "$node" "$start" "$end" "$cpus" "$mem"
    done | sort -k4,4n -k3,3
} > "${OUT}"

echo "wrote ${OUT}"
n=$(($(wc -l < "${OUT}") - 1))
echo "${n} jobs recorded"
echo
echo "--- summary by chrom (mean elapsed_min / max_rss_gb / n_complete) ---"
awk -F'\t' 'NR>1 && $5=="COMPLETED" {n[$4]++; tsec[$4]+=$7; tmem[$4]+=$9}
            END{
              printf "chrom\tn\tmean_min\tmean_gb\n"
              for(c in n) printf "%s\t%d\t%.1f\t%.1f\n", c, n[c], tsec[c]/n[c]/60, tmem[c]/n[c]
            }' "${OUT}" | sort -k1,1n
