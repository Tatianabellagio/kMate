#!/bin/bash
# Independently validate the base at Chr1:5870018 in 8 cactus founders.
# Cactus pangenome says 70/78 have T->A at this coord.
# xwu/hapFIRE says 5/231 have T->A (~2/78 cactus).
# Use minimap2 to align a 1kb window from TAIR10 to each founder genome
# and read the projected base.
set -e
MINIMAP=/global/home/users/tbellg/miniforge3/envs/sequencing_pipeline/bin/minimap2
SAMTOOLS=/global/home/users/tbellg/miniforge3/envs/sequencing_pipeline/bin/samtools
FA_DIR=/global/scratch/users/tbellg/pang/pang_1001gplus/20260209_Exposito-Alonso/chr_only
TAIR=$FA_DIR/TAIR10.chr.fa
OUT=/global/scratch/users/tbellg/kmate/scratch/atomize_test
mkdir -p $OUT/validation

POS=5870018
# 1kb window centered on the SNP from TAIR10
WIN_START=$((POS-500))
WIN_END=$((POS+500))

# Extract TAIR10 query window
$SAMTOOLS faidx $TAIR Chr1:${WIN_START}-${WIN_END} > $OUT/validation/tair10_query.fa
echo "TAIR10 base at coord $POS:"
$SAMTOOLS faidx $TAIR Chr1:${POS}-${POS} | tail -1

echo ""
echo "TAIR10 query window (1kb around SNP, offset of SNP in window = 500):"
echo -n "  Confirm: query[500] = "
tail -1 $OUT/validation/tair10_query.fa | cut -c501

# Pick 8 founders to test
# Pick from the cactus pangenome's T->A carriers (per atomize output) and T (REF) holders.
# We don't know exactly which founders are which without examining the per-sample GTs,
# but the cactus VCF said ~70/78 have T->A. So most random samples should match cactus's claim.
FOUNDERS=(100042 100043 100053 100108 100130 100208 100226 100302)

echo ""
echo "Aligning TAIR10 1kb window to each of 8 founder genomes via minimap2..."
echo "(Goal: find the founder's actual base at the coord that aligns to TAIR10 pos $POS)"
echo ""
printf "  %-8s  %-12s  %s\n" "FOUNDER" "FOUNDER_COORD" "FOUNDER_BASE"
echo "  ------------------------------"
for FID in "${FOUNDERS[@]}"; do
    FA=$FA_DIR/${FID}.chr.fa
    if [ ! -f "$FA" ]; then
        printf "  %-8s  %-12s  %s\n" "$FID" "(no FASTA)" "-"
        continue
    fi
    # Align query window to founder genome
    $MINIMAP -a -x asm5 --secondary=no $FA $OUT/validation/tair10_query.fa 2>/dev/null \
        | $SAMTOOLS view -F 0x900 2>/dev/null \
        | awk -F'\t' -v fid=$FID -v offset_in_query=500 '$2 != 4 {
            # Parse CIGAR to map query position 500 to reference position
            cigar = $6
            # rname (founder contig), pos (1-based start on founder)
            rname = $3
            pos = $2 == 16 ? $4 : $4  # start pos
            seq = $10
            flag = $2
            # For simplicity, use samtools tview approach: just print founder coord
            # walk through CIGAR
            q_pos = 0
            r_pos = pos
            # If reverse aligned (FLAG 16), we need to handle differently
            n_ops = 0
            while (match(cigar, /[0-9]+[MIDNSH=X]/)) {
                op = substr(cigar, RSTART+RLENGTH-1, 1)
                len = substr(cigar, RSTART, RLENGTH-1) + 0
                # query consume: M I S = X
                # ref consume:   M D N = X
                if (op == "S" || op == "H") {
                    cigar = substr(cigar, RSTART+RLENGTH)
                    continue
                }
                if (q_pos + len > offset_in_query) {
                    # Within this CIGAR op, find where in it
                    delta = offset_in_query - q_pos
                    if (op == "M" || op == "=" || op == "X") {
                        # query coord 500 maps to founder coord r_pos + delta
                        founder_coord = r_pos + delta
                        # Get founder base at that coord
                        char_idx = founder_coord - pos + 1  # in seq, accounting for cigar
                        # Easier: print founder coord; we'll extract base via samtools
                        printf "  %-8s  Chr1:%-10s  reverse_flag=%s\n", fid, founder_coord, (and(flag, 16) ? "REV" : "FWD")
                        exit
                    } else if (op == "I") {
                        # query has insertion, no ref coord
                        printf "  %-8s  INSERTION    -\n", fid
                        exit
                    } else if (op == "D" || op == "N") {
                        # query already consumed; this op consumes ref
                        # query position not advanced; continue
                    }
                }
                if (op == "M" || op == "=" || op == "X") { q_pos += len; r_pos += len }
                else if (op == "I") { q_pos += len }
                else if (op == "D" || op == "N") { r_pos += len }
                cigar = substr(cigar, RSTART+RLENGTH)
            }
        }'
done

echo ""
echo "Now check the base at each founder coord via samtools faidx..."
echo "  (Done separately because CIGAR parsing in awk was approximate)"
echo ""
echo "Per-founder true-base check via minimap2 PAF + samtools faidx:"
printf "  %-8s  %-15s  %-10s  %-10s\n" "FOUNDER" "FND_COORD" "FND_BASE" "TAIR10_BASE_EXPECTED"
echo "  --------------------------------------------------"
for FID in "${FOUNDERS[@]}"; do
    FA=$FA_DIR/${FID}.chr.fa
    [ ! -f "$FA" ] && continue
    # Use minimap2 PAF output for cleaner parsing
    # paftools.js mapped + tview would work too
    # Simpler: align and use samtools tview
    BAM=$OUT/validation/${FID}.bam
    if [ ! -f "$BAM" ]; then
      $MINIMAP -a -x asm5 --cs --secondary=no $FA $OUT/validation/tair10_query.fa 2>/dev/null \
        | $SAMTOOLS sort -o $BAM 2>/dev/null
      $SAMTOOLS index $BAM 2>/dev/null
    fi
    # Use samtools tview-like extract: get founder contig name from BAM
    HIT=$($SAMTOOLS view $BAM | head -1)
    if [ -z "$HIT" ]; then
        printf "  %-8s  no_hit\n" $FID
        continue
    fi
    FND_CONTIG=$(echo "$HIT" | cut -f3)
    FND_START=$(echo "$HIT" | cut -f4)
    FND_FLAG=$(echo "$HIT" | cut -f2)
    CIGAR=$(echo "$HIT" | cut -f6)

    # Use cs tag to get the alignment details
    CS_TAG=$(echo "$HIT" | grep -oE 'cs:Z:[^ 	]+' | head -1)

    # Use Python to parse CIGAR and project query pos 500 (= TAIR10 coord 5870018) to founder coord
    PROJ=$(python3 -c "
import sys, re
flag=int('$FND_FLAG')
start=int('$FND_START')
cigar='$CIGAR'
target_q=500   # 0-indexed query position
q_pos = 0
r_pos = start - 1  # 0-indexed founder pos
for length, op in re.findall(r'(\d+)([MIDNSH=X])', cigar):
    length = int(length)
    if op in ('S','H'):
        continue
    if op in ('M','=','X'):
        if q_pos + length > target_q:
            delta = target_q - q_pos
            print(f'{r_pos + delta + 1}')   # 1-based founder coord
            sys.exit(0)
        q_pos += length
        r_pos += length
    elif op == 'I':
        if q_pos + length > target_q:
            print('INSERTION')
            sys.exit(0)
        q_pos += length
    elif op in ('D','N'):
        r_pos += length
print('?')
")
    if [ -z "$PROJ" ] || [ "$PROJ" = "?" ] || [ "$PROJ" = "INSERTION" ]; then
        printf "  %-8s  %-15s  %-10s\n" "$FID" "$PROJ" "-"
        continue
    fi
    # Extract base at founder coord
    if [ $((FND_FLAG & 16)) -eq 16 ]; then
        # Reverse-strand alignment: founder base is reverse complement
        BASE_RAW=$($SAMTOOLS faidx $FA "$FND_CONTIG:$PROJ-$PROJ" 2>/dev/null | tail -1)
        BASE=$(echo "$BASE_RAW" | tr 'ACGTacgt' 'TGCAtgca')
    else
        BASE=$($SAMTOOLS faidx $FA "$FND_CONTIG:$PROJ-$PROJ" 2>/dev/null | tail -1)
    fi
    printf "  %-8s  %s:%-7s  %-10s  %s\n" "$FID" "$FND_CONTIG" "$PROJ" "$BASE" "(TAIR10 has T at $POS)"
done
