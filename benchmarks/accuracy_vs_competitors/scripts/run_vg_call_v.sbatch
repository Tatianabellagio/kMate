#!/bin/bash
#SBATCH --job-name=vg_call_v
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=04:00:00
#SBATCH --output=benchmarks/accuracy_vs_competitors/logs/vgcallv_%x_%j.out
#SBATCH --error=benchmarks/accuracy_vs_competitors/logs/vgcallv_%x_%j.out
set -eo pipefail
# Fair vg-SV comparator (BENCHMARK_DESIGN.md §3): genotype the KNOWN panel SVs
# against the pool's graph read-support (vg call -v), the SV analog of how the SNP
# benchmark scores vg at the shared panel sites. Emits AD at every panel SV ->
# AF = AD_alt/(AD_ref+AD_alt), placed at the panel svidx by sequence identity.
cd /global/scratch/users/tbellg/kmate
POOL=${1:?pool}
SVVCF=${2:-benchmarks/accuracy_vs_competitors/work/panel_sv_p80_Chr1.vcf.gz}
source ~/miniforge3/etc/profile.d/conda.sh && conda activate vg
G=/global/scratch/users/tbellg/pang/pang_1001gplus/pang/output/pang_1001gplus_82acc.d2.gbz
PACK=benchmarks/accuracy_vs_competitors/work/${POOL}.pack
OUT=benchmarks/accuracy_vs_competitors/work/${POOL}_vgcallv_Chr1.vcf
test -f "$PACK" || { echo "missing pack: $PACK"; exit 1; }
# NB: minigraph-cactus gbz has NO _alt_ paths, so the default VCFTraversalFinder
# finds 0 alt traversals. -z makes vg call find alleles via GBZ HAPLOTYPE traversals
# instead -- and our panel SV VCF was deconstructed from THIS gbz, so its alt
# traversals ARE gbz haplotype paths. (Without -z the run emits 0 variant lines.)
echo "[$(date)] vg call -v -z (genotype panel SVs via gbz haplotypes, TAIR10#0#Chr1) on $(hostname): $POOL"
vg call "$G" -k "$PACK" -v "$SVVCF" -z -p TAIR10#0#Chr1 -t "${SLURM_CPUS_PER_TASK:-8}" > "$OUT"
echo "[$(date)] done: $(grep -vc '^#' "$OUT") variant lines -> $OUT"
