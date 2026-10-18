#!/usr/bin/env python
"""Fit one simulated pool with kMate or hapFIRE.

Usage: python fit.py FIT_ID      (FIT_ID as listed by runs.py)

kMate: the production driver (src/per_sample_per_chrom.py: uniform k-mer weights, repeat guard and haploblock collapse
at their defaults) on the pool's reads, with the fit's k-mer panel, estimation unit and normalization, projecting
allele frequencies on its variant panel (or founder frequencies only, --h-only). Writes fit.tsv, fit.h_per_chrom.npz
(or fit.h_blocks_per_chrom.npz for window units), time.txt.
hapFIRE: reads aligned to TAIR10 (bwa mem | samtools sort), then hapFIRE on the founder set's GrENE-Net VCF. Writes
fit_ecotype_frequency.txt, fit_snp_frequency.txt, align_time.txt, hapfire_time.txt (the BAM is deleted).
Every step timed with /usr/bin/time -v, with at most FIT_THREADS threads (timed fits pinned to that many cores); DONE
written last. Output: config.FITS/FIT_ID/.
"""
import os, sys, subprocess
import config as C
from runs import all_fits


# Numerical libraries start one thread per core they can see; on a whole node (timed fits) that is 40 idle-spinning threads
# whose CPU time would count against the tool. Cap them at the fit's threads, and pin timed fits to that many cores.
THREAD_ENV = {v: str(C.FIT_THREADS) for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS")}
PIN = []


def run(cmd, env=None, **kw):
    cmd = PIN + cmd
    print("+", " ".join(cmd), flush=True)
    subprocess.run(cmd, check=True, env={**(env or os.environ), **THREAD_ENV}, **kw)


def unit_args(fit):
    if fit.unit == "chrom":
        return ["--unit", "chrom"]
    if fit.unit == "ld":
        return ["--unit", "ld", "--ld-blocks", C.VARIANT_PANELS[fit.variants]["ld_blocks"]]
    return ["--unit", "bp", "--window-bp", str(int(fit.unit.removeprefix("bp").removesuffix("k")) * 1000)]   # bp<N>k


def kmate(fit, reads, out):
    if fit.variants:
        pre = C.VARIANT_PANELS[fit.variants]["prefix"]
        project = ["--var-pa", f"{pre}.var_pa.npz", "--var-called", f"{pre}.var_called.npz", "--var-meta", f"{pre}.meta.npz"]
    else:
        project = ["--h-only"]
    run(["/usr/bin/time", "-v", "-o", f"{out}/time.txt", C.PYTHON, "-u", C.KMATE,
         "--kmer-pa-prefix", C.KMER_PANELS[fit.kmers], *project,
         "--reads", f"{reads}/r1.fq", f"{reads}/r2.fq", "--sample", fit.id, "--out", f"{out}/fit.tsv",
         "--threads", str(C.FIT_THREADS), "--chroms", C.CHROM, *unit_args(fit),
         "--kmer-weight", "uniform", "--normalize", fit.normalize])


def hapfire(fit, reads, out):
    t, bam = C.FIT_THREADS, f"{out}/sim.srt.bam"
    run(["/usr/bin/time", "-v", "-o", f"{out}/align_time.txt", "bash", "-c",
         f"set -euo pipefail; {C.BWA} mem -t {t} {C.HAPFIRE_REF} {reads}/r1.fq {reads}/r2.fq 2> {out}/bwa.log"
         f" | {C.SAMTOOLS} sort -@ {t} -o {bam} -; {C.SAMTOOLS} index {bam}"])
    env = dict(os.environ, PATH=f"{C.HAPFIRE}/src:{os.environ['PATH']}")
    run(["/usr/bin/time", "-v", "-o", f"{out}/hapfire_time.txt", C.HAPFIRE_PYTHON, f"{C.HAPFIRE}/hapFIRE.py",
         "-vcf", C.FOUNDER_SETS[fit.pool.founders]["hapfire_vcf"], "-bam", bam, "-fa", C.HAPFIRE_REF, "-output", "fit"],
        cwd=out, env=env)
    os.remove(bam); os.remove(f"{bam}.bai")


def main(fit_id):
    fit = all_fits()[fit_id]
    if fit.timed:
        PIN[:] = ["taskset", "-c", f"0-{C.FIT_THREADS - 1}"]
    out = f"{C.FITS}/{fit_id}"
    os.makedirs(out, exist_ok=True)
    {"kmate": kmate, "hapfire": hapfire}[fit.tool](fit, f"{C.POOLS}/{fit.pool.id}/reads", out)
    open(f"{out}/DONE", "w").close()


if __name__ == "__main__":
    main(sys.argv[1])
