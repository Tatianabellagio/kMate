#!/usr/bin/env python3
"""
Build TAIR10 -> raw-assembly coordinate liftover for each founder via minimap2.

Recombinant pools have crossover positions defined in TAIR10 (LD-block) coordinates,
but raw long-read assemblies are in their own coordinate frame (e.g. Chr1 = 34.6 Mb
vs TAIR10's 30.4 Mb). To stitch a mosaic from RAW assembly sequence we must project
each TAIR10 ancestry interval onto the founder's own coordinates. This script aligns
each founder's Chr to TAIR10 and stores the alignment blocks; sim_recomb_from_raw_
assemblies.py consumes them.

Output: <out-dir>/<asm_id>.<chrom>.blocks.tsv with columns
    t_start  t_end  q_start  q_end  strand     (0-based half-open; t=TAIR10, q=asm)
sorted by t_start. Reusable across all regimes/seeds (founder-only, sample-agnostic).

Run all founders (loop) or one (--asm-id, for SLURM array fan-out).
"""
import argparse, subprocess, sys, tempfile, os
from pathlib import Path


def extract_chrom(samtools, fa, chrom, out_fa):
    with open(out_fa, "w") as fh:
        rc = subprocess.run([samtools, "faidx", fa, chrom], stdout=fh)
    if rc.returncode != 0:
        raise RuntimeError(f"samtools faidx {fa} {chrom} failed")


def run_minimap2(minimap2, target_fa, query_fa, preset, threads):
    """target=TAIR10, query=assembly. Returns sorted list of primary PAF blocks."""
    cmd = [minimap2, "-c", "-x", preset, "-t", str(threads), target_fa, query_fa]
    p = subprocess.run(cmd, capture_output=True, text=True)
    if p.returncode != 0:
        raise RuntimeError(f"minimap2 failed:\n{p.stderr[:2000]}")
    blocks = []
    for line in p.stdout.splitlines():
        f = line.split("\t")
        if len(f) < 12:
            continue
        if not any(t == "tp:A:P" for t in f[12:]):  # primary alignments only
            continue
        qs, qe, strand = int(f[2]), int(f[3]), f[4]
        ts, te = int(f[7]), int(f[8])
        blocks.append((ts, te, qs, qe, strand))
    blocks.sort()
    return blocks


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--asm-dir", required=True)
    ap.add_argument("--rename", required=True, help="asm_id<TAB>ecotype_id per line")
    ap.add_argument("--tair10", required=True, help="TAIR10 multi-chrom FASTA")
    ap.add_argument("--chrom", default="Chr1")
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--preset", default="asm5", help="minimap2 -x preset (asm5/asm10/asm20)")
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--minimap2", default="/home/tbellagio/miniforge3/envs/sequencing_pipeline/bin/minimap2")
    ap.add_argument("--samtools", default="/home/tbellagio/miniforge3/envs/pang/bin/samtools")
    ap.add_argument("--asm-id", default=None, help="process only this asm_id (SLURM array)")
    args = ap.parse_args()

    out = Path(args.out_dir); out.mkdir(parents=True, exist_ok=True)
    asm_ids = []
    with open(args.rename) as fh:
        for line in fh:
            p = line.split()
            if p:
                asm_ids.append(p[0])
    if args.asm_id:
        asm_ids = [args.asm_id]

    tmp = Path(tempfile.mkdtemp(prefix="lift_"))
    tair_chr = tmp / f"tair10.{args.chrom}.fa"
    extract_chrom(args.samtools, args.tair10, args.chrom, str(tair_chr))

    for i, asm in enumerate(asm_ids):
        out_tsv = out / f"{asm}.{args.chrom}.blocks.tsv"
        if out_tsv.exists():
            print(f"[{i+1}/{len(asm_ids)}] {asm}: exists, skip", flush=True)
            continue
        asm_fa = Path(args.asm_dir) / f"{asm}.chr.fa"
        if not asm_fa.exists():
            sys.exit(f"missing {asm_fa}")
        asm_chr = tmp / f"{asm}.{args.chrom}.fa"
        extract_chrom(args.samtools, str(asm_fa), args.chrom, str(asm_chr))
        blocks = run_minimap2(args.minimap2, str(tair_chr), str(asm_chr),
                              args.preset, args.threads)
        cov = sum(te - ts for ts, te, *_ in blocks)
        minus = sum(1 for *_, st in blocks if st == "-")
        with open(out_tsv, "w") as fh:
            fh.write("t_start\tt_end\tq_start\tq_end\tstrand\n")
            for ts, te, qs, qe, st in blocks:
                fh.write(f"{ts}\t{te}\t{qs}\t{qe}\t{st}\n")
        os.remove(asm_chr)
        print(f"[{i+1}/{len(asm_ids)}] {asm}: {len(blocks)} blocks ({minus} on -), "
              f"TAIR10 {args.chrom} covered {cov:,} bp", flush=True)
    print("done", flush=True)


if __name__ == "__main__":
    main()
