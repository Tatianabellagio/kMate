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
import os, sys, subprocess, tempfile, shutil
from pathlib import Path


def _resolve_tool(name: str) -> str:
    """Locate an external binary robustly across cluster/env changes.

    Order: (1) the same env as the running Python interpreter (the driver runs
    under `kmate`, which ships jellyfish + samtools), (2) PATH, (3) bare name
    (let the shell resolve). Replaces the previous hardcoded Carnegie-era paths
    (`envs/pangenie`, `envs/sequencing_pipeline`) that broke after the migration.
    Override with env vars KMATE_JELLYFISH / KMATE_SAMTOOLS if needed.
    """
    env_override = os.environ.get(f"KMATE_{name.upper()}")
    if env_override:
        return env_override
    cand = os.path.join(os.path.dirname(sys.executable), name)
    if os.path.exists(cand):
        return cand
    return shutil.which(name) or name


JELLYFISH = _resolve_tool("jellyfish")
SAMTOOLS = _resolve_tool("samtools")


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


def build_kmer_db(
    reads: str | list[str],
    jf_db: str,
    k: int = 31,
    threads: int = 4,
    hash_size: str = "3G",
) -> str:
    """Build a canonical Jellyfish k-mer-count DB from reads, ONCE.

    Scanning the full read pool is the expensive step (~minutes for a ~5 GB
    pool). When the same reads are queried against several disjoint k-mer sets
    (e.g. one per chromosome), build the DB once with this and query each set
    with `query_kmer_db`, instead of re-scanning the reads N times. The counts
    are byte-identical to what `count_kmers_in_fasta`/`count_kmers_in_bam`
    return, because the hash (canonical, same `-s`) is identical.

    Args:
        reads:      FASTA/FASTQ path or list of them (plain or .gz), or a single
                    .bam (streamed via `samtools fastq -F 0x900`).
        jf_db:      output DB path (persistent; the caller owns/cleans it).
        k:          k-mer length (default 31).
        threads:    jellyfish threads.
        hash_size:  initial Jellyfish hash size (auto-grows).

    Returns:
        jf_db (the path written).
    """
    files = [reads] if isinstance(reads, str) else list(reads)
    if len(files) == 1 and files[0].endswith(".bam"):
        cmd = (
            f"{SAMTOOLS} fastq -F 0x900 {files[0]} 2>/dev/null | "
            f"{JELLYFISH} count -m {k} -s {hash_size} -t {threads} -C -o {jf_db} /dev/fd/0"
        )
    else:
        # zcat -f transparently handles plain + gzipped inputs.
        files_arg = " ".join(f'"{f}"' for f in files)
        cmd = (
            f"set -o pipefail; zcat -f {files_arg} | "
            f"{JELLYFISH} count -m {k} -s {hash_size} -t {threads} -C -o {jf_db} /dev/fd/0"
        )
    rc = subprocess.run(cmd, shell=True, executable="/bin/bash",
                        capture_output=True, text=True)
    if rc.returncode != 0:
        raise RuntimeError(f"jellyfish count failed:\n{rc.stderr}")
    return jf_db


def query_kmer_db(jf_db: str, query_kmers: list[str], k: int = 31) -> dict[str, int]:
    """Query a prebuilt Jellyfish DB for each k-mer's canonical count.

    Pairs with `build_kmer_db` for the count-once / query-per-chrom pattern.
    Returns a dict {kmer: count} with every query k-mer present (0 if absent).
    """
    with tempfile.TemporaryDirectory() as workdir:
        query_fa = os.path.join(workdir, "query.fa")
        _write_query_fasta(query_kmers, query_fa)
        rc = subprocess.run(
            f"{JELLYFISH} query {jf_db} -s {query_fa}",
            shell=True, capture_output=True, text=True
        )
        if rc.returncode != 0:
            raise RuntimeError(f"jellyfish query failed:\n{rc.stderr}")
        result: dict[str, int] = {}
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
