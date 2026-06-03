#!/bin/bash
#SBATCH --job-name=p231_skew
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=3:00:00
#SBATCH --output=/global/scratch/users/tbellg/kmate/benchmarks/p231/logs/08_skew_%j.out
#SBATCH --error=/global/scratch/users/tbellg/kmate/benchmarks/p231/logs/08_skew_%j.err
#
# 08 -- class-skewed n50 g0 pools for the h-imbalance grid (NO VISOR needed).
# For a g0 pool, each founder = one full haplotype; VISOR SHORtS just wgsim's each
# clone and pools by fraction. We reuse the per-founder Chr1 haplotypes already in
# the n231_g0 sim (haps/s_ind*/h1.fa) and wgsim each selected founder at cov10/50,
# pooling to a 50-founder pool with a chosen long-read(cactus)/short-read(PG) skew.
#
#   cact : 40 long-read + 10 short-read  -> long-read truth mass 0.80
#   pg   :  5 long-read + 45 short-read  -> long-read truth mass 0.10
#
# Then runs the 3 grid arms (raw / filt2 uniform / filt2 inv_mb) count-once.
#
# Usage: sbatch 08_build_run_skewed_p231.sh {cact|pg}
set -euo pipefail
SKEW=${1:?Usage: cact|pg}
[[ "$SKEW" == "cact" || "$SKEW" == "pg" ]] || { echo "SKEW must be cact|pg" >&2; exit 1; }

ROOT=/global/scratch/users/tbellg/kmate
CTRL=$ROOT/benchmarks/p231
PY=/global/home/users/tbellg/miniforge3/envs/kmate/bin/python
WGSIM=/global/home/users/tbellg/miniforge3/envs/kmate/bin/wgsim
DRIVER=$ROOT/src/per_sample_per_chrom.py
SRC_SIM=$CTRL/sims/cov10_n231_g0_s42_hotspots_p231_chr1   # supplies per-founder haplotypes

REGIME=n50_${SKEW}heavy
WORK=$CTRL/sims/cov10_${REGIME}_s42_hotspots_p231_chr1
READS=$WORK/reads
mkdir -p $READS $CTRL/logs

# ---- STAGE 1: select founders + wgsim per-founder reads, pool them ----
echo "[$(date)] STAGE 1: build $SKEW-heavy n50 pool with wgsim"
$PY - "$SKEW" "$SRC_SIM" "$WORK" "$WGSIM" <<'PY'
import sys, os, csv, json, random, subprocess
skew, src_sim, work, wgsim = sys.argv[1:5]
ROOT = "/global/scratch/users/tbellg/kmate"
sp  = json.load(open(f"{ROOT}/data/founder_split_cactus_pg.json"))
CAC = set(map(str, sp["cactus"]))
# founder -> haplotype fasta (via ancestry: ind_id -> founder ; haps/s_ind###/h1.fa)
f2hap = {}
for r in csv.DictReader(open(f"{src_sim}/ancestry.tsv"), delimiter="\t"):
    ind = r["ind_id"]                       # ind001 ..
    num = ind.replace("ind", "")
    hap = f"{src_sim}/haps/s_ind{num}/h1.fa"
    f2hap[str(r["founder"])] = hap
cac = sorted(f for f in f2hap if f in CAC)
pg  = sorted(f for f in f2hap if f not in CAC)
rng = random.Random(42)
n_c, n_p = (40, 10) if skew == "cact" else (5, 45)
sel = sorted(rng.sample(cac, n_c)) + sorted(rng.sample(pg, n_p))
assert len(sel) == 50

L, RL, COV, N = 30427671, 150, 10, 50            # Chr1 len, read len, total cov, n founders
npairs = int(round((COV / N) * L / (2 * RL)))    # per-founder read pairs (~20k)
reads = f"{work}/reads"
r1o, r2o = f"{reads}/r1.fq", f"{reads}/r2.fq"
open(r1o, "w").close(); open(r2o, "w").close()
for i, f in enumerate(sel):
    hap = f2hap[f]
    t1, t2 = f"{reads}/_t1.fq", f"{reads}/_t2.fq"
    cmd = [wgsim, "-N", str(npairs), "-1", str(RL), "-2", str(RL),
           "-e", "0.001", "-r", "0", "-R", "0", "-S", str(1000 + i), hap, t1, t2]
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    for t, o in ((t1, r1o), (t2, r2o)):
        with open(o, "a") as oh, open(t) as th:
            oh.write(th.read())
        os.remove(t)
# truth: uniform 1/50 over the 50 selected founders
with open(f"{work}/pool_weights.tsv", "w") as fh:
    fh.write("founder\tcount\tweight\n")
    for f in sel:
        fh.write(f"{f}\t1\t{1.0/N}\n")
with open(f"{work}/founders_selected.txt", "w") as fh:
    fh.write("\n".join(sel) + "\n")
print(f"  pool {skew}-heavy: {n_c} long-read + {n_p} short-read; {npairs} pairs/founder; "
      f"long-read truth mass {n_c/N:.2f}")
PY
echo "[$(date)] reads: $(du -h $READS/r1.fq | cut -f1) + $(du -h $READS/r2.fq | cut -f1)"

# ---- STAGE 2: 3 grid arms, count-once over the pool ----
RAW_PREFIX=$CTRL/data/kmer_pa_p231/kmer_pa
F2_PREFIX=$CTRL/data/kmer_pa_p231_filt2/kmer_pa
CN_VAR=$ROOT/panel/arch3/chr1/var_pa_231_arch3_chr1.var_pa.npz
CN_VAR_META=$ROOT/panel/arch3/chr1/var_pa_231_arch3_chr1.meta.npz

DB=${TMPDIR:-/tmp}/skew_${SKEW}.jf
echo "[$(date)] build_kmer_db -> $DB"
$PY - "$READS/r1.fq" "$READS/r2.fq" "$DB" <<'PY'
import sys; sys.path.insert(0, "/global/scratch/users/tbellg/kmate/src")
from kmer_count import build_kmer_db
build_kmer_db([sys.argv[1], sys.argv[2]], sys.argv[3], k=31, threads=4)
PY

run_arm () {  # arm_tag  prefix  weight
    local TAG=$1 PREFIX=$2 WEIGHT=$3
    local OUT_DIR=$CTRL/results/kmate_global_${TAG}/${REGIME}; mkdir -p "$OUT_DIR"
    local SAMPLE=p231_${TAG}_${REGIME}_cov10_s42
    echo "[$(date)] arm=$TAG weight=$WEIGHT"
    $PY -u $DRIVER --kmer-pa-prefix "$PREFIX" \
        --var-pa $CN_VAR --var-meta $CN_VAR_META \
        --reads "$READS/r1.fq" "$READS/r2.fq" --kmer-db "$DB" \
        --sample "$SAMPLE" --out "$OUT_DIR/${SAMPLE}.tsv" \
        --threads 4 --chroms Chr1 --block-mode global --kmer-weight $WEIGHT
}
run_arm raw_raw     "$RAW_PREFIX" uniform
run_arm filt2u_raw  "$F2_PREFIX"  uniform
run_arm filt2mb_raw "$F2_PREFIX"  inv_mb
rm -f "$DB"
echo "[$(date)] DONE skewed $SKEW -> results/kmate_global_*/$REGIME/"
