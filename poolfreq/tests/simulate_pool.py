"""
Simulate a whole-Chr1 pool from N panel founders at known proportions.

Uses wgsim (samtools) to generate paired-end reads from each founder's Chr1,
scaled by that founder's pool weight. Output is a single combined paired-end
FASTQ that can be fed into our pool-seq pipeline.

Truth: h_true = pool weights (vector of length F = 82, summing to 1)

Usage:
    python simulate_pool.py --founders f1,f2,... --weights w1,w2,...
                            --coverage 30 --read-len 150 --out pool_sim
"""
from __future__ import annotations
import argparse, os, subprocess, sys, gzip
from pathlib import Path
import numpy as np

WGSIM = "/home/tbellagio/miniforge3/envs/sequencing_pipeline/bin/wgsim"
SAMTOOLS = "/home/tbellagio/miniforge3/envs/sequencing_pipeline/bin/samtools"

# Default founder FASTA dir + chromosome
FOUNDER_DIR = "/home/tbellagio/scratch/pang/pang_1001gplus/20260209_Exposito-Alonso/chr_only"


def simulate_one(founder_fa: str, n_reads: int, read_len: int, out_prefix: str,
                 chrom: str = "Chr1", err_rate: float = 0.005):
    """Simulate paired-end reads from a specific chromosome of a founder FASTA.

    wgsim outputs two FASTQs: <prefix>_1.fq and <prefix>_2.fq (uncompressed).
    """
    # Extract just Chr1 to a temp file
    chr_fa = f"{out_prefix}_{chrom}.fa"
    rc = subprocess.run(
        f"{SAMTOOLS} faidx {founder_fa} {chrom} > {chr_fa}",
        shell=True, check=True, capture_output=True, text=True,
    )
    # wgsim
    fq1 = f"{out_prefix}_1.fq"
    fq2 = f"{out_prefix}_2.fq"
    cmd = (f"{WGSIM} -N {n_reads} -1 {read_len} -2 {read_len} "
           f"-e {err_rate} -r 0 -R 0 -X 0 -d 400 -s 50 "
           f"{chr_fa} {fq1} {fq2}")
    rc = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if rc.returncode != 0:
        raise RuntimeError(f"wgsim failed for {founder_fa}: {rc.stderr}")
    os.remove(chr_fa)
    if os.path.exists(chr_fa + ".fai"):
        os.remove(chr_fa + ".fai")
    return fq1, fq2


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--founder-dir", default=FOUNDER_DIR)
    ap.add_argument("--founders", required=True,
                    help="comma-separated list of assembly IDs (e.g., 100042,100043,...)")
    ap.add_argument("--weights", required=True,
                    help="comma-separated weights (must sum to 1; will be normalized)")
    ap.add_argument("--coverage", type=float, default=30.0,
                    help="total target coverage")
    ap.add_argument("--read-len", type=int, default=150)
    ap.add_argument("--chrom", default="Chr1")
    ap.add_argument("--chrom-len-bp", type=int, default=30_427_671,
                    help="chromosome length in bp (Chr1 default)")
    ap.add_argument("--out", required=True, help="output prefix")
    ap.add_argument("--threads", type=int, default=4)
    args = ap.parse_args()

    founders = args.founders.split(",")
    weights = np.array([float(w) for w in args.weights.split(",")])
    assert len(founders) == len(weights)
    weights = weights / weights.sum()

    out_dir = Path(args.out).parent
    out_dir.mkdir(parents=True, exist_ok=True)

    # Compute n_reads per founder
    bp_per_pair = 2 * args.read_len  # paired
    total_bp = args.chrom_len_bp * args.coverage
    total_pairs = int(total_bp / bp_per_pair)
    print(f"Simulating {args.chrom} pool ({args.chrom_len_bp:,}bp) at {args.coverage}× → {total_pairs:,} read pairs")
    print(f"  founders: {len(founders)}")
    print(f"  weights: min={weights.min():.4f}, max={weights.max():.4f}, sum={weights.sum():.4f}")

    fastqs = []
    for f, w in zip(founders, weights):
        fa = os.path.join(args.founder_dir, f"{f}.chr.fa")
        if not os.path.exists(fa):
            print(f"  skipping {f}: FASTA not found at {fa}")
            continue
        n = max(1, int(total_pairs * w))
        prefix = f"{args.out}_per_founder/{f}"
        Path(prefix).parent.mkdir(parents=True, exist_ok=True)
        fq1, fq2 = simulate_one(fa, n, args.read_len, prefix, chrom=args.chrom)
        fastqs.append((fq1, fq2, f, w, n))
        print(f"  {f}: w={w:.4f}, n_pairs={n:,} → {os.path.basename(fq1)}, {os.path.basename(fq2)}")

    # Concatenate into pool FASTQs
    print(f"\nConcatenating {len(fastqs)} per-founder FASTQs into pool...")
    combined1 = f"{args.out}_pool_1.fq"
    combined2 = f"{args.out}_pool_2.fq"
    with open(combined1, "w") as o1, open(combined2, "w") as o2:
        for fq1, fq2, _, _, _ in fastqs:
            with open(fq1) as f1: o1.write(f1.read())
            with open(fq2) as f2: o2.write(f2.read())
    # Cleanup per-founder
    for fq1, fq2, _, _, _ in fastqs:
        os.remove(fq1); os.remove(fq2)
    # Compress
    for f in [combined1, combined2]:
        subprocess.run(["gzip", "-f", f], check=True)
    print(f"  {combined1}.gz")
    print(f"  {combined2}.gz")

    # Save truth (h_true vector)
    truth_path = f"{args.out}_truth.tsv"
    with open(truth_path, "w") as out:
        out.write("founder\tweight\tn_pairs\n")
        for fq1, fq2, fid, w, n in fastqs:
            out.write(f"{fid}\t{w:.6f}\t{n}\n")
    print(f"\nWrote truth: {truth_path}")
    print("DONE")


if __name__ == "__main__":
    main()
