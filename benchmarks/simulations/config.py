"""Paths and fixed settings shared by simulate.py, fit.py, submit.py and score.py."""

K = "/global/scratch/projects/fc_moilab/tbellg/kmate"
B = f"{K}/benchmarks"
HERE = f"{B}/simulations"
POOLS = f"{HERE}/pools"            # one directory per simulated pool
FITS = f"{HERE}/fits"              # one directory per kMate fit
RESULTS = f"{HERE}/results"        # scored tables read by the benchmark notebook
LOGS = f"{HERE}/logs"

ENV = "/global/home/users/tbellg/miniforge3/envs/kmate/bin"
PYTHON, WGSIM, BWA, SAMTOOLS = f"{ENV}/python", f"{ENV}/wgsim", f"{ENV}/bwa", f"{ENV}/samtools"
HAPFIRE_PYTHON = "/global/home/users/tbellg/miniforge3/envs/hapfire/bin/python"
HAPFIRE = f"{K}/external/HapFIRE"
HAPFIRE_REF = "/global/scratch/users/tbellg/pang/ref_xing/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa"   # bwa-indexed, contig "1"
MOSAICS = f"{K}/sims/scripts/make_recomb_mosaics.py"
TRUTH = f"{K}/sims/scripts/compute_recomb_truth.py"
KMATE = f"{K}/src/per_sample_per_chrom.py"
CHROM = "Chr1"

# Founder sets a pool is drawn from. genomes: <founder>.chr.fa (./. written as N) and, with chains, <founder>.chain
# (TAIR10 -> the founder's own coordinates); order: the founder list whose order, with the pool seed, fixes the draw;
# truth: the variant panels the per-record truth is computed on.
# kMate's pools come from the panel's genomes (indels and SVs included); hapFIRE's from its own GrENE-Net SNP-only genomes
# (TAIR10 coordinates, no chains), scored against its panel (hapfire_vcf). Same order and seed: both tools get the same
# founders.
P231_ORDER = f"{K}/data/kmer_pa_231_arch3_filt2inv/kmer_pa_Chr1.meta.npz"
P78_ORDER = f"{B}/p80/data/kmer_pa_p78_raw/kmer_pa_Chr1.meta.npz"
FOUNDER_SETS = {
    "p231": dict(genomes=f"{B}/p231/fastas_231n", chains=True, order=P231_ORDER, truth=["p231"]),
    "p78": dict(genomes=f"{B}/p80/fastas_78n", chains=True, order=P78_ORDER, truth=["p78", "p231"]),
    "greneNet231": dict(genomes=f"{B}/vs_hapfire/greneNet_fastas", chains=False, order=P231_ORDER, truth=[],
                        hapfire_vcf=f"{B}/vs_hapfire/work/greneNet_chr1.vcf.gz"),
    "greneNet78": dict(genomes=f"{B}/vs_hapfire/greneNet_fastas", chains=False, order=P78_ORDER, truth=[],
                       hapfire_vcf=f"{B}/vs_hapfire/work/greneNet_chr1_p78.vcf.gz"),      # build_greneNet_p78_vcf.sh
}

# Variant panels: <prefix>.{var_pa,var_called,meta}.npz, and the r2 >= 0.1 LD blocks used by `--unit ld`.
VARIANT_PANELS = {
    "p231": dict(prefix=f"{K}/panel/arch3/chr1/var_pa_231_arch3_chr1", ld_blocks=f"{B}/p231/data/ld_blocks_r2_0.10_Chr1.tsv"),
    "p78": dict(prefix=f"{B}/p80/data/var_pa_p78", ld_blocks=f"{B}/p80/data/var_pa_p78.ld_blocks_r2_0.10_Chr1.tsv"),
}

# k-mer panels (<prefix>_Chr1.kmer_pa.npz). "filt" = the private k-mer filter (k-mers carried by 2 to F-1 founders).
KMER_PANELS = {
    "p231": f"{K}/data/kmer_pa_231_arch3_filt2inv/kmer_pa",              # production
    "p231_unfilt": f"{B}/p231/data/kmer_pa_p231/kmer_pa",                # founder-type figure: every k-mer ...
    "p231_filt": f"{B}/p231/data/kmer_pa_p231_filt2inv/kmer_pa",         # ... and the same panel filtered
    "p78_unfilt": f"{B}/p80/data/kmer_pa_p78_raw/kmer_pa",
    "p78_filt": f"{B}/p80/data/kmer_pa_p78_filt2inv/kmer_pa",
}

# Long-read founders: the 78 samples of the cactus VCF; the other 153 of the 231 were genotyped from short reads.
LONG_READ_VCF = f"{K}/panel/arch3/chr1/cactus_78_chr1_haploid.vcf.gz"

# Paired-end reads (wgsim): VISOR SHORtS's settings, used for every earlier simulation.
READ_LEN, INSERT, INSERT_SD = 150, 500, 50
ERROR, MUTATION, INDEL_FRAC, INDEL_EXT, MAX_N = 0.001, 0.001, 0.15, 0.25, 0.05

# Slurm
SBATCH = ["--account=co_moilab", "--partition=savio4_htc", "--qos=savio_lowprio"]
SIM_RES = ["--cpus-per-task=4", "--mem=48G", "--time=4:00:00"]
FIT_RES = ["--cpus-per-task=8", "--mem=96G", "--time=6:00:00"]
SCORE_RES = ["--cpus-per-task=4", "--mem=96G", "--time=6:00:00"]
# timed fits: a whole savio3 node with 40 cores and 96 GB (one CPU type), nothing else on it; the fit still uses FIT_THREADS
TIMED_SBATCH = ["--account=co_moilab", "--partition=savio3", "--qos=savio_lowprio", "--constraint=savio3_c40", "--exclusive"]
TIMED_RES = ["--time=6:00:00"]
FIT_THREADS = 8
