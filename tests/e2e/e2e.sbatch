#!/bin/bash
#SBATCH --job-name=kmate_e2e
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --qos=moilab_htc4_normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=12:00:00
#SBATCH --output=e2e_%j.log
# End-to-end check of a kMate release as a new user gets it: fresh env, PyPI install,
# Chr1 panel built from the 135-assembly Cactus graph, one pool run, then compared
# against the production (arch3) matrices by tests/e2e/compare.py.
#
# Usage:  mkdir -p $W && cd $W && cp <repo>/tests/e2e/compare.py . && \
#         sbatch --export=ALL,KMATE_VERSION=0.1.4,W=$PWD <repo>/tests/e2e/e2e.sbatch
#
# 0.1.3 result (2026-10-05, job 39654836, ~4 h): build-index 1h18 / 24 GB, decompose
# Chr1 42 min / 8 GB, build-kmer-pa 1h32 / 22 GB, run 2 min / 16 GB. Index identical to
# pang_135_haploid in all bubbles (cols 1-4; production also filled the overhang
# column); for the 73 founders mapped to production, var_pa / var_called (157M cells)
# and kmer_pa (664M cells) 0 mismatches; 2,154,416 of 2,154,423 production records
# present. SEEDMIX S1 alt_freq vs the 231 panel r=0.963 -- different founder sets.
set -euo pipefail
W=${W:?set W to the work dir}; cd $W
V=${KMATE_VERSION:?set KMATE_VERSION}
PANG=/global/scratch/users/tbellg/pang/pang_1001gplus/pang_all/output
REF=/global/scratch/users/tbellg/pang/pang_1001gplus/20260209_Exposito-Alonso/chr_only/TAIR10.chr.fa
PROD=/global/home/users/tbellg/kmate
R1=/global/scratch/users/tbellg/pang/grenenet_reads/seed_mix/S1-1.1_P.fq.gz
R2=/global/scratch/users/tbellg/pang/grenenet_reads/seed_mix/S1-1.2_P.fq.gz
T=${SLURM_CPUS_PER_TASK:-8}
stamp() { echo "[$(date +%T)] $*"; }

stamp "0. fresh env (documented install: conda deps + PyPI kmate)"
export CONDA_PKGS_DIRS=/local/job${SLURM_JOB_ID}/pkgs CONDA_SOLVER=libmamba
[ -x $W/env/bin/kmate ] || /global/home/users/tbellg/miniforge3/bin/mamba create -y -q -p $W/env -c conda-forge -c bioconda \
    "python>=3.9,<3.13" numpy scipy pysam "kmer-jellyfish >=2.3.1 py*" samtools bcftools pip pandas > env.log 2>&1
export PATH=$W/env/bin:/usr/bin:/bin
[ -x $W/env/bin/kmate ] || pip install -q "kmate==$V"
kmate --version; python --version; python -c "import dna_jellyfish; print('dna_jellyfish OK')"
kmate selftest 2>&1 | tail -1

stamp "1. build-index, genome-wide graph VCF (as production), on node-local disk"
# build-index keeps its jellyfish DB next to --out and memory-maps it; on Lustre the
# random lookups crawl (job 39652179 stalled). Build on /local, then copy back.
L=/local/job${SLURM_JOB_ID}/kmate_e2e; mkdir -p $L/index index
/usr/bin/time -v kmate build-index --vcf $PANG/pang_1001gplus_all.vcf.gz --ref $REF --out $L/index/ours \
    -k 31 --haploid --jellyfish-threads $T --jellyfish-hash 3000000000 > index.log 2> index.err
grep -E "Elapsed|Maximum resident" index.err
cp $L/index/ours_*_kmers.tsv.gz index/ && rm -rf $L

stamp "2. decompose Chr1 (--gfa)"
bcftools view -r Chr1 $PANG/pang_1001gplus_all.vcf.gz -Oz -o graph_chr1.vcf.gz; bcftools index -t -f graph_chr1.vcf.gz
/usr/bin/time -v kmate decompose --gfa $PANG/pang_1001gplus_all.gfa.gz --genotyped-vcf graph_chr1.vcf.gz \
    --out panel_chr1.vcf.gz --threads $T > decompose.log 2> decompose.err
grep -E "Elapsed|Maximum resident" decompose.err

stamp "3. build-kmer-pa + build-var-pa (Chr1)"
mkdir -p kmer_pa var_pa
/usr/bin/time -v kmate build-kmer-pa --kmers index/ours_Chr1_kmers.tsv.gz --vcf panel_chr1.vcf.gz --ref $REF \
    --chrom Chr1 --out kmer_pa/kmer_pa_Chr1 --treat-missing-as-n --filter-production > kpa.log 2> kpa.err
grep -E "Elapsed|Maximum resident" kpa.err
kmate build-var-pa --vcf panel_chr1.vcf.gz --chrom Chr1 --out var_pa/var_pa_Chr1 > vpa.log 2>&1

stamp "4. run SEEDMIX S1 Chr1 against the new panel and the production panel"
kmate build-kmer-db --reads $R1 $R2 --out /local/job${SLURM_JOB_ID}/s1.jf --threads $T > db.log 2>&1
/usr/bin/time -v kmate run --kmer-pa-prefix kmer_pa/kmer_pa --var-pa-prefix var_pa/var_pa \
    --reads $R1 $R2 --kmer-db /local/job${SLURM_JOB_ID}/s1.jf --sample S1_new --out S1_new_Chr1.tsv \
    --chroms Chr1 --threads $T > run_new.log 2> run_new.err
grep -E "Elapsed|Maximum resident" run_new.err
kmate run --kmer-pa-prefix $PROD/data/kmer_pa_231_arch3_filt2inv/kmer_pa \
    --var-pa $PROD/panel/arch3/chr1/var_pa_231_arch3_chr1.var_pa.npz \
    --var-called $PROD/panel/arch3/chr1/var_pa_231_arch3_chr1.var_called.npz \
    --var-meta $PROD/panel/arch3/chr1/var_pa_231_arch3_chr1.meta.npz \
    --reads $R1 $R2 --kmer-db /local/job${SLURM_JOB_ID}/s1.jf --sample S1_prod --out S1_prod_Chr1.tsv \
    --chroms Chr1 --threads $T > run_prod.log 2>&1
rm -f /local/job${SLURM_JOB_ID}/s1.jf

stamp "5. compare against production"
python $W/compare.py
stamp "DONE"
