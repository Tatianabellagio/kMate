"""Haploidize a diploid VCF: one allele per GT, as kMate's panel builders require.

`build_kmer_pa` reconstructs ONE sequence per founder and `build_var_pa` applies a
carrier rule, so both refuse diploid genotypes rather than silently picking an
allele. This converts:

    0/0, 0|0  ->  0      homozygous reference
    1/1, 1|1  ->  1      homozygous alternate
    0/1, 1/0, 0|1, 1|0  ->  .    HETEROZYGOUS BECOMES MISSING (see below)
    ./., .|., ., anything else  ->  .

**Het -> missing, not a coin flip.** The panel represents inbred founders, where a
call that is heterozygous is more likely a genotyping artefact than real diploidy.
Marking it missing lets the AF projection exclude that founder at that record
(`var_called` is 0 there) instead of inventing a REF or ALT call; following the
precedent used for *Arabidopsis* panels (Arouisse et al. 2020, Plant J 102:872-882,
imputation for inbred Arabidopsis accessions). Ported from the awk step in
`panel/arch3/chr*/jobA2_*.sh`; the counts it reports let you check the rate.
"""
from __future__ import annotations

import gzip
import subprocess
import sys

_HOM_REF = {"0/0", "0|0"}
_HOM_ALT = {"1/1", "1|1"}
_HET = {"0/1", "1/0", "0|1", "1|0"}
_MISS = {"./.", ".|.", "."}


def haploidize_line(line: str, counts: dict) -> str:
    f = line.rstrip("\n").split("\t")
    if len(f) <= 9:
        return line
    f[8] = "GT"
    for i in range(9, len(f)):
        g = f[i].split(":")[0]
        if g in _HOM_REF:
            f[i] = "0"
        elif g in _HOM_ALT:
            f[i] = "1"
        elif g in _HET:
            f[i] = "."
            counts["het_to_missing"] += 1
        else:
            f[i] = "."
            if g not in _MISS:
                counts["other_to_missing"] += 1
        counts["cells"] += 1
    return "\t".join(f) + "\n"


def haploidize_vcf(in_vcf: str, out_vcf: str, bcftools: str = "bcftools",
                   threads: int = 4) -> dict:
    """Stream in_vcf -> out_vcf, haploidizing genotypes. Returns the counters."""
    counts = {"cells": 0, "het_to_missing": 0, "other_to_missing": 0}
    reader = subprocess.Popen([bcftools, "view", "--threads", str(max(1, threads // 2)),
                               "-Ov", in_vcf], stdout=subprocess.PIPE, text=True)
    writer = None
    try:
        if out_vcf.endswith(".gz"):
            writer = subprocess.Popen([bcftools, "view", "-Oz", "-o", out_vcf, "-"],
                                      stdin=subprocess.PIPE, text=True)
            sink = writer.stdin
        else:
            sink = open(out_vcf, "w")
        for line in reader.stdout:
            sink.write(line if line.startswith("#") else haploidize_line(line, counts))
        sink.close()
    finally:
        reader.wait()
        if writer is not None:
            writer.wait()

    pct = 100.0 * counts["het_to_missing"] / counts["cells"] if counts["cells"] else 0.0
    sys.stderr.write(
        f"[haploidize] {counts['cells']:,} genotype cells; "
        f"het->missing {counts['het_to_missing']:,} ({pct:.2f}%); "
        f"other->missing {counts['other_to_missing']:,}\n"
    )
    return counts
