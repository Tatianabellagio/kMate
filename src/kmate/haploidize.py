"""Make a diploid VCF haploid, because kMate's panel is a set of HAPLOTYPES.

`build_kmer_pa` reconstructs exactly ONE sequence per founder and `build_var_pa`
applies a carrier rule, so both refuse diploid genotypes rather than silently
picking an allele. That is a property of the model, not of any species: a panel
column is a haplotype.

How you get there depends on what your founders ARE. Two policies:

`--het missing` (default) — **for inbred founders** (e.g. *Arabidopsis* accessions,
MAGIC/RIL founders, maize NAM parents). These are expected to be homozygous, so a
heterozygous call is more likely a genotyping artefact than real diploidy:

    0/0, 0|0 -> 0      1/1, 1|1 -> 1      het -> .      ./., . -> .
    0, 1 (already haploid) -> unchanged

Marking het missing lets the AF projection *exclude* that founder at that record
(`var_called` is 0 there) instead of inventing a REF or ALT call. Precedent for
inbred *Arabidopsis* panels: Arouisse et al. (2020), Plant J 102:872-882. The
reported flip rate is the check — a few percent is normal, a lot means your
founders are not what this policy assumes.

`--het split` — **for outbred, PHASED founders** (e.g. HPRC-style assemblies, any
diploid individual with phased haplotypes). Discarding hets here would throw away
half the data and bias AF. Instead each sample becomes two haplotype columns:

    sample with 0|1   ->   sample.h1 = 0 ,  sample.h2 = 1

The founder axis doubles (N samples -> 2N haplotypes) and kMate estimates a
frequency per haplotype, which is what a founder-mixture model means for an
outbred panel. Requires phased genotypes: an unphased het (`0/1`) cannot be split
and is an error, because guessing the phase would fabricate haplotypes.

Ported from the awk step in `panel/arch3/chr*/jobA2_*.sh`, which implemented the
`missing` policy.
"""
from __future__ import annotations

import subprocess
import sys

_HOM_REF = {"0/0", "0|0"}
_HOM_ALT = {"1/1", "1|1"}
_HET = {"0/1", "1/0", "0|1", "1|0"}
_MISS = {"./.", ".|.", "."}
_HAPLOID = {"0", "1"}          # already one allele (e.g. Minigraph-Cactus on haploid assemblies)


def _gt_to_haploid(g: str, counts: dict) -> str:
    if g in _HAPLOID:
        return g
    if g in _HOM_REF:
        return "0"
    if g in _HOM_ALT:
        return "1"
    if g in _HET:
        counts["het_to_missing"] += 1
        return "."
    if g not in _MISS:
        counts["other_to_missing"] += 1
    return "."


def _gt_to_pair(g: str, counts: dict, where: str) -> tuple[str, str]:
    """Split one diploid GT into two haplotype alleles. Requires phasing for hets."""
    core = g.split(":")[0]
    if core in _MISS:
        return ".", "."
    if "|" in core:
        a, b = core.split("|", 1)
    elif "/" in core:
        a, b = core.split("/", 1)
        if a != b:
            raise ValueError(
                f"--het split needs PHASED genotypes; got unphased heterozygote "
                f"{core!r} at {where}. Phase the VCF, or use --het missing if your "
                f"founders are inbred lines."
            )
    else:
        # one allele: nothing to split. Copying it into h1 and h2 would invent a
        # second, identical founder haplotype.
        raise ValueError(
            f"--het split needs DIPLOID genotypes; got haploid {core!r} at {where}. "
            f"A haploid panel is already one column per haplotype: drop --haploidize."
        )
    counts["cells"] += 1
    return (a if a not in ("", ".") else "."), (b if b not in ("", ".") else ".")


def haploidize_line(line: str, counts: dict, het: str = "missing") -> str:
    f = line.rstrip("\n").split("\t")
    if len(f) <= 9:
        return line
    f[8] = "GT"
    if het == "split":
        where = f"{f[0]}:{f[1]}"
        out = f[:9]
        for i in range(9, len(f)):
            a, b = _gt_to_pair(f[i], counts, where)
            out.extend((a, b))
        return "\t".join(out) + "\n"
    for i in range(9, len(f)):
        f[i] = _gt_to_haploid(f[i].split(":")[0], counts)
        counts["cells"] += 1
    return "\t".join(f) + "\n"


def _split_header(line: str) -> str:
    f = line.rstrip("\n").split("\t")
    out = f[:9]
    for s in f[9:]:
        out.extend((f"{s}.h1", f"{s}.h2"))
    return "\t".join(out) + "\n"


def haploidize_vcf(in_vcf: str, out_vcf: str, bcftools: str = "bcftools",
                   threads: int = 4, het: str = "missing") -> dict:
    """Stream in_vcf -> out_vcf making genotypes haploid. Returns counters."""
    if het not in ("missing", "split"):
        raise ValueError("het must be 'missing' or 'split'")
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
            if line.startswith("#CHROM") and het == "split":
                sink.write(_split_header(line))
            elif line.startswith("#"):
                sink.write(line)
            else:
                sink.write(haploidize_line(line, counts, het=het))
        sink.close()
    finally:
        reader.wait()
        if writer is not None:
            writer.wait()

    if het == "split":
        sys.stderr.write(
            f"[haploidize] --het split: {counts['cells']:,} diploid genotypes -> "
            f"2 haplotype columns each (founder axis doubled)\n")
    else:
        pct = 100.0 * counts["het_to_missing"] / counts["cells"] if counts["cells"] else 0.0
        sys.stderr.write(
            f"[haploidize] {counts['cells']:,} genotype cells; "
            f"het->missing {counts['het_to_missing']:,} ({pct:.2f}%); "
            f"other->missing {counts['other_to_missing']:,}\n")
        if pct > 10:
            sys.stderr.write(
                "[haploidize] WARNING: >10% of calls were heterozygous. --het missing "
                "assumes INBRED founders; if yours are outbred and phased, use "
                "--het split instead or you are discarding real data.\n")
    return counts
