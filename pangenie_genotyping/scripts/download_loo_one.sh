#!/bin/bash
#SBATCH --job-name=loo_dl
#SBATCH --partition=bse
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH --time=4:00:00
#SBATCH --output=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/pangenie_genotyping/logs/loo_dl_%A_%a.out
#SBATCH --error=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/pangenie_genotyping/logs/loo_dl_%A_%a.err

# =============================================================================
# download_loo_one.sh
# Same logic as download_one.sh but reads from loo_ena_manifest.tsv and writes
# to data/loo_raw_fastqs/<ecotype>/. Indexed 1..N over manifest data rows.
# =============================================================================
set -euo pipefail

BASE=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/pangenie_genotyping
MANIFEST=$BASE/data/loo_ena_manifest.tsv
OUT_DIR=$BASE/data/loo_raw_fastqs
mkdir -p $OUT_DIR

IDX=${SLURM_ARRAY_TASK_ID:-1}
LINE=$(sed -n "$((IDX+1))p" $MANIFEST)
[ -n "$LINE" ] || { echo "ERROR: no row at $IDX" >&2; exit 1; }

ECOTYPE=$(echo "$LINE" | cut -f1)
RUN=$(echo "$LINE" | cut -f2)
FASTQ_FTP=$(echo "$LINE" | cut -f4)
FASTQ_MD5=$(echo "$LINE" | cut -f5)

echo "[$(date)] LOO ecotype=$ECOTYPE run=$RUN files=$FASTQ_FTP"

ECO_DIR=$OUT_DIR/$ECOTYPE
mkdir -p $ECO_DIR

IFS=';' read -ra URLS <<< "$FASTQ_FTP"
IFS=';' read -ra MD5S <<< "$FASTQ_MD5"

for i in "${!URLS[@]}"; do
    URL="${URLS[$i]}"
    EXP_MD5="${MD5S[$i]:-}"
    [[ "$URL" =~ ^(ftp|http)s?:// ]] || URL="https://$URL"
    OUT="$ECO_DIR/$(basename "$URL")"

    if [ -f "$OUT" ] && [ -n "$EXP_MD5" ]; then
        ACT_MD5=$(md5sum "$OUT" | cut -d' ' -f1)
        if [ "$ACT_MD5" = "$EXP_MD5" ]; then
            echo "  [skip] $OUT (md5 ok)"; continue
        else
            echo "  [redo] $OUT (md5 mismatch)"; rm -f "$OUT"
        fi
    elif [ -f "$OUT" ]; then
        echo "  [skip] $OUT (exists, no md5 check)"; continue
    fi

    echo "  [dl] $URL"
    wget --quiet --continue --tries=5 --timeout=120 -O "$OUT" "$URL"

    if [ -n "$EXP_MD5" ]; then
        ACT_MD5=$(md5sum "$OUT" | cut -d' ' -f1)
        [ "$ACT_MD5" = "$EXP_MD5" ] || { echo "ERROR: md5 mismatch on $OUT"; exit 1; }
        echo "  [ok] md5 verified"
    fi
done

du -sh "$ECO_DIR"
echo "[$(date)] DONE LOO ecotype=$ECOTYPE"
