#!/usr/bin/env python
"""Extract per-chrom founder SNP genotypes from the GrENE-net short-read VCF.

Produces, per chromosome, an npz the SV-SNP LD script consumes:
  snp_pos  [nSNP] int64        SNP positions (TAIR10)
  geno     [nSNP x 231] f32    ALT-allele dosage (0/1/2), missing mean-imputed
  founders [231]               sample (ecotype) IDs, in VCF column order

The short-read VCF chroms are named 1..5; we tag the output by Chr1..Chr5.
"""
from __future__ import annotations
import argparse, os, subprocess, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lib

BCF = "/global/home/users/tbellg/miniforge3/bin/bcftools"
GT2D = {"0/0": 0, "0|0": 0, "0/1": 1, "1/0": 1, "0|1": 1, "1|0": 1,
        "1/1": 2, "1|1": 2}   # everything else (./., .|., .) -> NaN


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--vcf", required=True)
    ap.add_argument("--chrom-vcf", required=True, help="chrom name in the VCF (e.g. 1)")
    ap.add_argument("--chrom-out", required=True, help="output tag (e.g. Chr1)")
    ap.add_argument("--out", default=f"{lib.GEA}/sv_snp_ld")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    founders = subprocess.run([BCF, "query", "-l", args.vcf],
                              capture_output=True, text=True).stdout.split()
    n = len(founders)
    # stream POS + GT for this chrom
    p = subprocess.Popen([BCF, "query", "-r", args.chrom_vcf,
                          "-f", "%POS[\t%GT]\n", args.vcf],
                         stdout=subprocess.PIPE, text=True, bufsize=1 << 20)
    pos = []
    rows = []
    cache = {}
    for line in p.stdout:
        f = line.rstrip("\n").split("\t")
        pos.append(int(f[0]))
        gts = f[1:]
        row = np.empty(n, np.float32)
        for j, g in enumerate(gts):
            v = cache.get(g)
            if v is None:
                v = GT2D.get(g, np.nan); cache[g] = v
            row[j] = v
        rows.append(row)
    p.wait()
    geno = np.array(rows, np.float32)
    pos = np.array(pos, np.int64)
    # mean-impute missing per SNP (column = founder; row = SNP)
    if np.isnan(geno).any():
        mu = np.nanmean(geno, axis=1, keepdims=True)
        mu[np.isnan(mu)] = 0.0
        miss = np.isnan(geno)
        geno[miss] = np.broadcast_to(mu, geno.shape)[miss]
    np.savez(f"{args.out}/shortread_geno_{args.chrom_out}.npz",
             snp_pos=pos, geno=geno, founders=np.array(founders))
    print(f"{args.chrom_out}: {len(pos):,} SNPs x {n} founders -> "
          f"{args.out}/shortread_geno_{args.chrom_out}.npz "
          f"(missing imputed; mean dosage {np.nanmean(geno):.3f})")


if __name__ == "__main__":
    main()
