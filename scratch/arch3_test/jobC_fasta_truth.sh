#!/bin/bash
#SBATCH --job-name=archC_fasta
#SBATCH --partition=bse
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/scratch/arch3_test/jobC_fasta_%j.out
#SBATCH --error=/carnegie/nobackup/scratch/tbellagio/hapfire_sv/scratch/arch3_test/jobC_fasta_%j.err
set -euo pipefail

# Job C: open-loop FASTA-truth validation at 50 random positions in test region.
# For each position, read each founder's consensus FASTA at coord, compute "true" AC.
# Compare to merged biallelic AC from Job B.

cd /carnegie/nobackup/scratch/tbellagio/hapfire_sv/scratch/arch3_test

PY=/home/tbellagio/miniforge3/envs/hapfm/bin/python
MERGED=merged_231_test_final.vcf.gz

[ -s $MERGED ] || { echo "ERROR: missing $MERGED from Job B"; exit 1; }

# Try to locate per-founder consensus FASTAs
FASTA_DIR=""
for candidate in \
  /carnegie/nobackup/scratch/tbellagio/hapfire_sv/sims/visor_freqk/unimputed_fastas_v3 \
  /carnegie/nobackup/scratch/tbellagio/hapfire_sv/unimputed_fastas_v3 \
  /home/tbellagio/scratch/hapfire_sv/unimputed_fastas_v3; do
  if [ -d "$candidate" ]; then
    FASTA_DIR=$candidate
    break
  fi
done
[ -n "$FASTA_DIR" ] || { echo "ERROR: cannot locate unimputed_fastas_v3 directory"; exit 1; }
echo "Using FASTA dir: $FASTA_DIR"
echo "Sample FASTAs: $(ls $FASTA_DIR/*.fa* 2>/dev/null | head -3)"

$PY <<EOF
import gzip, random
import subprocess
from pathlib import Path

# Step 1: load 50 random SNP biallelics from merged VCF in test region
target_recs = []
all_snps = []
with gzip.open("$MERGED", "rt") as f:
    samples = None
    for line in f:
        if line.startswith("##"): continue
        if line.startswith("#CHROM"):
            samples = line.rstrip().split("\t")[9:]
            continue
        parts = line.rstrip().split("\t")
        if parts[0] != "Chr1": continue
        ref = parts[3]; alt = parts[4]
        if len(ref) == 1 and len(alt) == 1:
            all_snps.append(parts)

random.seed(42)
sample_idx = random.sample(range(len(all_snps)), min(50, len(all_snps)))
print(f"Loaded {len(all_snps):,} biallelic SNPs in test region; sampling {len(sample_idx)} for FASTA truth.")

# Step 2: for each sampled record, get truth AC from FASTAs
fasta_dir = Path("$FASTA_DIR")
fasta_paths = {p.stem: p for p in fasta_dir.glob("*.fa")}
print(f"Found {len(fasta_paths)} per-founder FASTAs")

# Index 5 sample FASTAs to validate access works
sample_fastas = list(fasta_paths.values())[:5]
print("Sample FASTA paths:")
for p in sample_fastas:
    print(f"  {p}")

# Skip if no FASTAs found
if len(fasta_paths) == 0:
    print("WARN: no founder FASTAs found — cannot do FASTA-truth validation")
    print("Skipping Job C.")
    exit(0)

# pyfaidx for fast random access
try:
    from pyfaidx import Fasta
except ImportError:
    print("ERROR: pyfaidx not installed; trying samtools faidx fallback")
    import subprocess
    SAMTOOLS = "/home/tbellagio/miniforge3/envs/BIOS424/bin/samtools"
    Fasta = None

# Open all founder FASTAs (lazy)
print("Building per-position truth at 50 sampled positions...")
results = []
for i, idx in enumerate(sample_idx):
    rec = all_snps[idx]
    pos = int(rec[1]); ref = rec[3]; alt = rec[4]
    pipeline_gts = rec[9:]
    pipeline_carriers = sum(1 for g in pipeline_gts if g == "1")
    pipeline_called = sum(1 for g in pipeline_gts if g in ("0", "1"))

    # Get FASTA-truth: for each sample, look up base at coord
    truth_carriers = 0; truth_called = 0
    for s in samples:
        # try matching FASTA by sample name
        candidates = [p for k, p in fasta_paths.items() if k == s or k.endswith("_"+s) or k.startswith(s+"_") or s in k]
        if not candidates:
            continue
        fa_path = candidates[0]
        # Use samtools faidx to extract base
        try:
            out = subprocess.check_output(["${BIOS424_SAMTOOLS:-/home/tbellagio/miniforge3/envs/BIOS424/bin/samtools}",
                                            "faidx", str(fa_path), f"Chr1:{pos}-{pos}"],
                                            stderr=subprocess.DEVNULL).decode().strip().split("\n")
            base = out[1].upper() if len(out) > 1 else "N"
        except Exception as e:
            base = "N"
        if base == "N":
            continue
        truth_called += 1
        if base == alt.upper():
            truth_carriers += 1

    pipeline_af = pipeline_carriers / max(pipeline_called, 1)
    truth_af = truth_carriers / max(truth_called, 1)
    diff = pipeline_af - truth_af
    results.append((pos, ref, alt, pipeline_carriers, pipeline_called, truth_carriers, truth_called, pipeline_af, truth_af, diff))
    if (i+1) % 10 == 0:
        print(f"  done {i+1}/{len(sample_idx)}")

# Summary
import statistics
diffs = [r[9] for r in results]
abs_diffs = [abs(d) for d in diffs]
print()
print(f"=== FASTA-truth vs Arch 3 pipeline AC comparison ({len(results)} SNPs) ===")
print(f"  mean |Δ AF|:  {statistics.mean(abs_diffs):.5f}")
print(f"  median |Δ AF|: {statistics.median(abs_diffs):.5f}")
print(f"  max |Δ AF|:   {max(abs_diffs):.5f}")
print(f"  # records with |Δ AF| > 0.05: {sum(1 for d in abs_diffs if d > 0.05)}")
print(f"  # records with |Δ AF| > 0.10: {sum(1 for d in abs_diffs if d > 0.10)}")
print()
print("Sample of 10 records (sorted by |Δ AF|):")
results.sort(key=lambda r: abs(r[9]), reverse=True)
print(f"  {'pos':>10} {'REF':>4} {'ALT':>4} {'pipe_AC':>8} {'pipe_AN':>8} {'truth_AC':>9} {'truth_AN':>9} {'pipe_AF':>8} {'truth_AF':>9} {'diff':>8}")
for r in results[:10]:
    print(f"  {r[0]:>10} {r[1]:>4} {r[2]:>4} {r[3]:>8} {r[4]:>8} {r[5]:>9} {r[6]:>9} {r[7]:>8.4f} {r[8]:>9.4f} {r[9]:>+8.4f}")
EOF

echo
echo "[$(date)] DONE Job C"
