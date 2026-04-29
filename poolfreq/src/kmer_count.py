"""
K-mer counter for pool-seq BAMs against a fixed query k-mer set.

Wraps Jellyfish 2.x. Two-step:
  1. jellyfish count: build hash from BAM-derived FASTQ (samtools fastq pipe)
  2. jellyfish query: get counts for each k-mer in the query list

Returns a dict {kmer_str: count}. Reverse complements are handled via
jellyfish's canonical mode (--canonical, -C).

Usage:
    counts = count_kmers_in_bam(bam_path, query_kmers, k=31)

For batch use across many samples sharing the same query set, prefer
count_kmers_batch() which reuses the query file.
"""
from __future__ import annotations
import os, subprocess, tempfile, shutil
from pathlib import Path

JELLYFISH = "/home/tbellagio/miniforge3/envs/pangenie/bin/jellyfish"
SAMTOOLS = "/home/tbellagio/miniforge3/envs/sequencing_pipeline/bin/samtools"


def _write_query_fasta(kmers: list[str], path: str) -> None:
    """Write k-mers as a FASTA for jellyfish query input."""
    with open(path, "w") as f:
        for i, k in enumerate(kmers):
            f.write(f">k{i}\n{k}\n")


def count_kmers_in_bam(
    bam: str,
    query_kmers: list[str],
    k: int = 31,
    threads: int = 4,
    hash_size: str = "100M",
    workdir: str | None = None,
    keep_workdir: bool = False,
) -> dict[str, int]:
    """Count occurrences of each query k-mer in a BAM's reads.

    Args:
        bam:           path to indexed BAM
        query_kmers:   list of k-mer strings (length k)
        k:             k-mer length (default 31, PanGenie default)
        threads:       jellyfish threads
        hash_size:     initial Jellyfish hash size (e.g. "100M" or "1G")
        workdir:       optional working dir; if None, uses a tempdir
        keep_workdir:  if True, don't delete the working dir on exit (debug)

    Returns:
        dict {kmer: count}. Counts are sums over forward+reverse complement
        because jellyfish operates in canonical mode.
    """
    assert all(len(s) == k for s in query_kmers), f"all k-mers must have length {k}"
    assert os.path.exists(bam), f"missing BAM: {bam}"

    if workdir is None:
        workdir = tempfile.mkdtemp(prefix="kmer_count_")
    Path(workdir).mkdir(parents=True, exist_ok=True)

    jf_db = os.path.join(workdir, "counts.jf")
    query_fa = os.path.join(workdir, "query.fa")
    _write_query_fasta(query_kmers, query_fa)

    try:
        # 1. Stream BAM reads to jellyfish via samtools fastq.
        #    Use generator process with bash for the pipe; samtools fastq writes
        #    interleaved fastq to stdout, which jellyfish reads from /dev/fd/0.
        cmd_count = (
            f"{SAMTOOLS} fastq -F 0x900 {bam} 2>/dev/null | "
            f"{JELLYFISH} count -m {k} -s {hash_size} -t {threads} -C -o {jf_db} /dev/fd/0"
        )
        rc = subprocess.run(cmd_count, shell=True, executable="/bin/bash",
                             capture_output=True, text=True)
        if rc.returncode != 0:
            raise RuntimeError(f"jellyfish count failed:\n{rc.stderr}")

        # 2. Query each k-mer.
        cmd_query = f"{JELLYFISH} query {jf_db} -s {query_fa}"
        rc = subprocess.run(cmd_query, shell=True, capture_output=True, text=True)
        if rc.returncode != 0:
            raise RuntimeError(f"jellyfish query failed:\n{rc.stderr}")

        # Output format: "<kmer> <count>" per line
        result: dict[str, int] = {}
        for line in rc.stdout.strip().split("\n"):
            if not line: continue
            parts = line.split()
            if len(parts) >= 2:
                result[parts[0]] = int(parts[1])
        # Make sure every query k-mer is in the dict (jellyfish skips zeros sometimes)
        for kmer in query_kmers:
            result.setdefault(kmer, 0)
        return result
    finally:
        if not keep_workdir:
            shutil.rmtree(workdir, ignore_errors=True)


def count_kmers_in_fasta(
    fasta: str | list[str],
    query_kmers: list[str],
    k: int = 31,
    threads: int = 4,
    hash_size: str = "100M",
) -> dict[str, int]:
    """Same as count_kmers_in_bam but takes a FASTA/FASTQ (or list of them) directly.

    Supports gzipped input via jellyfish's `--generator` flag (we use a small
    bash wrapper to zcat-pipe gz files).
    """
    if isinstance(fasta, str):
        files = [fasta]
    else:
        files = list(fasta)
    has_gz = any(f.endswith(".gz") for f in files)
    with tempfile.TemporaryDirectory() as workdir:
        jf_db = os.path.join(workdir, "counts.jf")
        query_fa = os.path.join(workdir, "query.fa")
        _write_query_fasta(query_kmers, query_fa)

        if has_gz:
            # Stream all inputs through zcat (handles plain + gzipped)
            files_arg = " ".join(f'"{f}"' for f in files)
            cmd = (f"set -o pipefail; for f in {files_arg}; do "
                   f"  if [[ \"$f\" == *.gz ]]; then zcat \"$f\"; else cat \"$f\"; fi; "
                   f"done | {JELLYFISH} count -m {k} -s {hash_size} -t {threads} "
                   f"-C -o {jf_db} /dev/fd/0")
            rc = subprocess.run(cmd, shell=True, executable="/bin/bash",
                                 capture_output=True, text=True)
        else:
            files_arg = " ".join(f'"{f}"' for f in files)
            rc = subprocess.run(
                f"{JELLYFISH} count -m {k} -s {hash_size} -t {threads} -C "
                f"-o {jf_db} {files_arg}",
                shell=True, capture_output=True, text=True
            )
        if rc.returncode != 0:
            raise RuntimeError(f"jellyfish count failed:\n{rc.stderr}")
        rc = subprocess.run(
            f"{JELLYFISH} query {jf_db} -s {query_fa}",
            shell=True, capture_output=True, text=True
        )
        if rc.returncode != 0:
            raise RuntimeError(f"jellyfish query failed:\n{rc.stderr}")

        result = {}
        for line in rc.stdout.strip().split("\n"):
            parts = line.split()
            if len(parts) >= 2:
                result[parts[0]] = int(parts[1])
        for kmer in query_kmers:
            result.setdefault(kmer, 0)
        return result


def estimate_coverage(bam: str, genome_size_bp: int = 119_146_348) -> float:
    """Estimate per-base sequencing coverage (λ) from a BAM.

    λ = total mapped read bases / genome size.
    Default genome_size is TAIR10 nuclear (sum of Chr1-5).

    For k-mer-level coverage, multiply by (read_length - k + 1) / read_length.
    """
    cmd = f"{SAMTOOLS} stats -F 0x900 {bam} | grep -E '^SN\\s+(bases mapped|average length)' | head -2"
    rc = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    bases_mapped = None
    for line in rc.stdout.strip().split("\n"):
        if "bases mapped" in line:
            try:
                bases_mapped = int(line.split()[3])
            except (IndexError, ValueError):
                pass
    if bases_mapped is None:
        # fallback: count reads × average length
        cmd2 = f"{SAMTOOLS} view -c -F 0x900 {bam}"
        n_reads = int(subprocess.run(cmd2, shell=True, capture_output=True, text=True).stdout.strip())
        bases_mapped = n_reads * 150  # assumed read length
    return bases_mapped / genome_size_bp
