#!/usr/bin/env python3
"""Extract the greneNet_chr1.vcf.gz SNP genotype matrix (231 founders x N_records)
as a founder-dosage array, cached to work/greneNet_gt_matrix.npz. Feeds
score_af_hapfire_vs_kmate.py's truth-AF computation (pool_weights @ dosage).

This VCF has NO missing genotypes (fully Beagle-imputed, confirmed by a full-file
grep for './.'/'.|.') so dosage = allele-count/2 with no MAR denominator needed
(unlike sims/scripts/compute_recomb_truth.py's arch3 MAR truth, which does need
one for the raw/atomized var_pa panel where PG founders can be `.` at SV records).

Run with the `basic` env (bcftools resolved via absolute path below):
    /global/home/users/tbellg/miniforge3/envs/basic/bin/python \
        benchmarks/speed_vs_hapfire/scripts/extract_greneNet_gt_matrix.py
"""
import subprocess
import numpy as np
import pandas as pd

ROOT = "/global/scratch/users/tbellg/kmate"
VCF = f"{ROOT}/benchmarks/speed_vs_hapfire/work/greneNet_chr1.vcf.gz"
OUT = f"{ROOT}/benchmarks/speed_vs_hapfire/work/greneNet_gt_matrix.npz"
BCFTOOLS = "/global/home/users/tbellg/miniforge3/bin/bcftools"

GT_DOSAGE = {"0|0": 0.0, "0/0": 0.0, "1|1": 1.0, "1/1": 1.0,
             "0|1": 0.5, "1|0": 0.5, "0/1": 0.5, "1/0": 0.5}

print("querying sample names...", flush=True)
samples = subprocess.run([BCFTOOLS, "query", "-l", VCF], capture_output=True,
                         text=True, check=True).stdout.split()
F = len(samples)
print(f"  {F} founders", flush=True)

print("querying CHROM/POS/GT (streaming through bcftools)...", flush=True)
proc = subprocess.Popen([BCFTOOLS, "query", "-f", "%CHROM\t%POS[\t%GT]\n", VCF],
                        stdout=subprocess.PIPE, text=True)
df = pd.read_csv(proc.stdout, sep="\t", header=None, names=["chrom", "pos"] + samples)
ret = proc.wait()
if ret != 0:
    raise RuntimeError(f"bcftools query failed: {ret}")
print(f"  read {len(df):,} records x {F} founders", flush=True)

dosage = df[samples].replace(GT_DOSAGE).to_numpy(dtype=np.float32)
n_missing = np.isnan(dosage).sum()
print(f"  dosage matrix built, missing/unrecognized GT calls: {n_missing}", flush=True)
if n_missing:
    print("  WARNING: missing GTs present -- truth_af will be biased for those "
          "founders/records unless handled downstream", flush=True)

np.savez_compressed(OUT, chrom=df["chrom"].to_numpy(dtype=str),
                    pos=df["pos"].to_numpy(dtype=np.int64),
                    founders=np.array(samples), dosage=dosage)
print(f"wrote {OUT}", flush=True)
