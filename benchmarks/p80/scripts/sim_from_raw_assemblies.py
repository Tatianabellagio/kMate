#!/usr/bin/env python3
"""
Simulate pool-seq reads directly from raw long-read assemblies, bypassing VISOR.

Why this exists: VISOR's `bedfile` requirement extracts a region into a tmpfa
before calling pywgsim. With consensus FASTAs all aligned to TAIR10, that's a
no-op, but it forces raw assemblies to be truncated to TAIR10 coordinates.
This script calls `pywgsim.wgsim.core(ref=<raw.chr.fa>, ...)` directly on each
founder's full Chr1 sequence (29.5-34.6 Mb across the p80 panel) and concatenates
per-founder R1/R2 into a single pool R1/R2 FASTQ pair.

Read budget: match the existing consensus-FASTA g0 sim's per-individual read
count to keep the only varying factor the FASTA source.

Read parameters: 150 bp paired, err 0.001, frag mean 500 stdev 50 (pywgsim
defaults), to match VISOR's call signature in 06_run_sim_p80.sh.
"""

import argparse, os, sys, time, gzip, subprocess, csv
from pathlib import Path
from pywgsim import wgsim


def asm_chr1_len(fai_path: Path, chrom: str = "Chr1") -> int:
    for line in fai_path.read_text().splitlines():
        f = line.split("\t")
        if f[0] == chrom:
            return int(f[1])
    raise RuntimeError(f"{chrom} not in {fai_path}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ancestry", required=True,
                    help="Existing ancestry.tsv to reuse (5 cols, header: ind_id chrom start end founder)")
    ap.add_argument("--rename", required=True,
                    help="Tab-separated: asm_id<TAB>ecotype_id (samples_80_rename.tsv)")
    ap.add_argument("--asm-dir", required=True,
                    help="Directory holding <asm_id>.chr.fa raw assemblies")
    ap.add_argument("--chrom", default="Chr1")
    ap.add_argument("--coverage", type=float, default=10.0,
                    help="Pool coverage in × of TAIR10 Chr1 (matches VISOR --coverage)")
    ap.add_argument("--tair10-len", type=int, default=30427671,
                    help="TAIR10 Chr1 length for total-read-budget calibration")
    ap.add_argument("--read-len", type=int, default=150)
    ap.add_argument("--err-rate", type=float, default=0.001)
    ap.add_argument("--frag-mean", type=int, default=500)
    ap.add_argument("--frag-stdev", type=int, default=50)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--out", required=True, help="Output directory (will hold reads/r1.fq, reads/r2.fq)")
    args = ap.parse_args()

    out = Path(args.out); out.mkdir(parents=True, exist_ok=True)
    reads_dir = out / "reads"; reads_dir.mkdir(exist_ok=True)
    work = out / "per_ind"; work.mkdir(exist_ok=True)

    with open(args.ancestry) as fh:
        ancestry = list(csv.DictReader(fh, delimiter="\t"))
    # g0-only invariant: one row per individual. A recombinant (g>=1) ancestry has
    # multiple (start,end,founder) segments per ind_id and would be silently
    # mis-simulated here (each segment treated as a whole individual).
    ids = [r["ind_id"] for r in ancestry]
    if len(ids) != len(set(ids)):
        sys.exit("ancestry has >1 row per ind_id (recombinant g>=1); this g0-only "
                 "script needs one row per individual")
    eco_to_asm = {}
    with open(args.rename) as fh:
        for parts in csv.reader(fh, delimiter="\t"):
            if len(parts) >= 2:
                eco_to_asm[str(parts[1])] = str(parts[0])

    # Per-individual read budget: total pool read pairs = cov * TAIR10_len / (2 * read_len),
    # split equally across N individuals (mirrors VISOR with uniform clonefraction).
    n_indiv = len(ancestry)
    total_pairs = round(args.coverage * args.tair10_len / (2 * args.read_len))
    per_ind_pairs = total_pairs // n_indiv
    residual = total_pairs - per_ind_pairs * n_indiv
    print(f"[plan] N={n_indiv} indivs, total={total_pairs} pairs, per_ind={per_ind_pairs}, "
          f"first_ind_extra={residual} (absorbs fp residual)")

    r1_out = reads_dir / "r1.fq"
    r2_out = reads_dir / "r2.fq"
    if r1_out.exists(): r1_out.unlink()
    if r2_out.exists(): r2_out.unlink()

    # Per-individual prov record
    prov_rows = []
    t0 = time.time()
    for i, row in enumerate(ancestry):
        ind_id = row["ind_id"]
        founder = str(row["founder"])
        asm = eco_to_asm.get(founder)
        if asm is None:
            sys.exit(f"ecotype {founder} not in rename map")
        asm_fa = Path(args.asm_dir) / f"{asm}.chr.fa"
        fai = Path(str(asm_fa) + ".fai")
        if not fai.exists():
            subprocess.run(["samtools", "faidx", str(asm_fa)], check=True)
        chr_len = asm_chr1_len(fai, args.chrom)

        # Make a per-ind tmpfa with ONLY Chr1 (wgsim samples uniformly across the whole ref,
        # so multi-chrom would bleed reads). Use samtools faidx to extract.
        tmpfa = work / f"{ind_id}.chr1.fa"
        if not tmpfa.exists():
            subprocess.run(["samtools", "faidx", str(asm_fa), args.chrom, "-o", str(tmpfa)], check=True)

        n_pairs = per_ind_pairs + (residual if i == 0 else 0)
        r1_ind = work / f"{ind_id}.r1.fq"
        r2_ind = work / f"{ind_id}.r2.fq"

        wgsim.core(
            ref=str(tmpfa),
            r1=str(r1_ind),
            r2=str(r2_ind),
            err_rate=args.err_rate,
            mut_rate=0.0,        # we want clean reads — variants come from the assembly itself
            indel_frac=0.0,
            indel_ext=0.0,
            max_n=0.05,
            is_hap=1,             # haploid: each FASTA represents one haplotype
            N=n_pairs,
            dist=args.frag_mean,
            stdev=args.frag_stdev,
            size_l=args.read_len,
            size_r=args.read_len,
            is_fixed=0,
            seed=args.seed + i,    # unique seed per individual
        )

        # Append to pool R1/R2 (binary append — files are small enough)
        with open(r1_out, "ab") as fout, open(r1_ind, "rb") as fin:
            fout.write(fin.read())
        with open(r2_out, "ab") as fout, open(r2_ind, "rb") as fin:
            fout.write(fin.read())
        # Don't keep per-ind FASTQs once concatenated (saves disk)
        r1_ind.unlink(); r2_ind.unlink()

        prov_rows.append({
            "ind_id": ind_id, "founder": founder, "asm_id": asm,
            "chr1_len": chr_len, "n_pairs": n_pairs,
        })
        if (i + 1) % 5 == 0 or i == n_indiv - 1:
            dt = time.time() - t0
            print(f"  [{i+1}/{n_indiv}] {ind_id} <- {asm} (chr1_len={chr_len:,}, N={n_pairs:,}) [{dt:.0f}s]")

    with open(out / "provenance.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["ind_id", "founder", "asm_id", "chr1_len", "n_pairs"],
                           delimiter="\t")
        w.writeheader()
        w.writerows(prov_rows)
    print(f"[done] pool reads at {r1_out} / {r2_out}")
    print(f"       per-ind provenance at {out/'provenance.tsv'}")


if __name__ == "__main__":
    main()
