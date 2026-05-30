#!/usr/bin/env python -u
"""
Standalone re-implementation of PanGenie-index's kmers.tsv.gz output.

Goal: produce the per-bubble unique-k-mer table that downstream tools
(`build_kmer_cn.py`, future hapwater) consume, WITHOUT depending on PanGenie
to build the rest of its 9 GB index.

Validation against PanGenie v4.2.1 on a 2-accession test panel (set_02_rep1,
Chr1:1-2Mb, 3252 variants → 1985 bubbles): 99.4% byte-equivalent (1973/1985
bubbles identical). The remaining 12 bubbles differ only because of a quirk
in PG's `stepwise_unique_kmers` (stepwiseuniquekmercomputer.cpp:9-32) where
`counts[current_kmer] += 1` after the enumeration loop emits a phantom k-mer
with a leading `A` from the uninitialized `jellyfish::mer_dna("")` high bits
when the overhang sequence is exactly k-1 chars. We do NOT reproduce that
bug — our output is semantically more correct.

Reference: Ebler et al. 2022, Nat Genet — "Pangenome-based genome inference
allows efficient and accurate genotyping across a wide spectrum of variant
classes". Algorithm reproduced from `external/pangenie-tools/src/`
under MIT License.

Algorithm (mirrors PanGenie-index):

  1. Parse VCF + reference FASTA.
  2. Cluster variants within (k-1) bp of each other into bubbles. Each bubble's
     "allele" is a haplotype path through the cluster — defined by the unique
     combinations of (allele_idx_at_var1, ..., allele_idx_at_varN) observed
     across founder paths (plus the all-REF path if --add-reference).
  3. Reconstruct each unique allele's sequence: left_flank (k-1 ref bases) +
     allele_v1 + inner_flank_1 (ref between v1 and v2) + ... + allele_vN +
     right_flank (k-1 ref bases). Write to a temp FASTA (path_segments.fasta).
  4. Build a genome-wide canonical-k-mer count via jellyfish on the FASTA.
  5. For each bubble: per-allele canonical k-mers; keep those with
     genomic_count == 1 AND local_count == 1 (unique to one allele within the
     bubble, and not elsewhere in the panel). Round-robin select up to
     --cap-biallelic (default 16) or --cap-multiallelic (default 32) k-mers
     per allele, with per-bubble total cap of max(nr_paths, 301).
  6. Overhang: per bubble, for each ±2k flanking region, select up to
     --overhang-cap (default 12) unique-genome-wide k-mers.
  7. Write `<prefix>_<chrom>_kmers.tsv.gz` with the same 5-column schema as
     PanGenie-index.

Parameters not in PG:
  --haploid             accept haploid GT (single integer, no `|`). Each sample
                        contributes 1 path instead of 2. Useful for cactus
                        pangenome output.
  --cap-biallelic       per-allele kmer cap in biallelic bubbles (PG: 16)
  --cap-multiallelic    per-allele kmer cap in multi-allelic bubbles (PG: 32)
  --overhang-cap        max kmers per flanking side (PG: 12)
  --no-add-reference    omit the all-REF synthetic path. PG includes it by default.

Validation strategy: byte-compare emitted TSV against PG's on a small panel
(see panel/pangenie_index/archive/test_vs_pangenie.py), and the full-scale
index-parity diff in panel/pangenie_index/scripts/diff_index_vs_pg.py.

Performance (2 Mb Chr1, 1985 bubbles, 4 threads, jellyfish hash 100M):
  PanGenie-index (C++):  ~7 s wall, ~878 MB RSS
  This script (Python):  ~21 s wall, ~846 MB RSS

The dominant bottleneck is Python-level per-kmer iteration in
select_unique_kmers and enumerate_fwd_kmers. A C extension (or cython) wrapping
those two loops could likely match PG's speed; not needed for current usage.
"""
from __future__ import annotations
import argparse
import gzip
import os
import subprocess
import sys
import tempfile
from collections import OrderedDict, defaultdict
from typing import Iterable, Iterator

import pysam


# ---------------------------------------------------------------------------
# k-mer utilities
# ---------------------------------------------------------------------------

_COMPLEMENT = bytes.maketrans(b"ACGTacgt", b"TGCAtgca")

def revcomp(s: bytes) -> bytes:
    return s.translate(_COMPLEMENT)[::-1]

def canonical(km: bytes) -> bytes:
    rc = revcomp(km)
    return km if km < rc else rc

def enumerate_fwd_kmers(seq: str, k: int) -> Iterator[bytes]:
    """Yield FORWARD k-mers from seq, skipping any window containing non-ACGT.

    PanGenie's `unique_kmers` stores fwd k-mers (built via jellyfish mer_dna
    shift_left, not auto-canonicalized). Canonicalization only happens at the
    moment of querying the genomic jellyfish hash."""
    seq = seq.upper().encode("ascii")
    n = len(seq)
    valid_set = set(b"ACGT")
    last_bad = -1
    for i in range(n):
        if seq[i] not in valid_set:
            last_bad = i
        if i >= k - 1 and last_bad <= i - k:
            yield seq[i - k + 1 : i + 1]


# ---------------------------------------------------------------------------
# VCF + bubble construction
# ---------------------------------------------------------------------------

class Variant:
    __slots__ = ("chrom", "start", "end", "ref", "alts", "paths")
    def __init__(self, chrom, start_0based, ref, alts, paths):
        self.chrom = chrom
        self.start = start_0based
        self.end = start_0based + len(ref)
        self.ref = ref.upper()
        self.alts = [a.upper() for a in alts]
        # paths: list of allele indices, one per haplotype path
        self.paths = paths

    def nr_alleles(self) -> int:
        return 1 + len(self.alts)

    def allele_sequence(self, idx: int) -> str:
        if idx == 0:
            return self.ref
        return self.alts[idx - 1]


def _open_vcf(vcf_path: str):
    """Open VCF for reading; auto-detect .vcf.gz / .bgz (gzip)."""
    if vcf_path.endswith((".gz", ".bgz")):
        return gzip.open(vcf_path, "rt")
    return open(vcf_path, "r")


def parse_vcf_variants(vcf_path: str, ref_fasta: pysam.FastaFile, k: int, haploid: bool, add_reference: bool):
    """Yield Variants from VCF. Validates REF against reference, skips invalid alts.

    Mirrors PanGenie GraphBuilder:
      * skips variants <2k from chrom ends
      * skips variants whose ALTs contain non-ACGT
      * stores paths as (n_samples * ploidy [+ 1 ref]) integer allele indices
    """
    # Map of chrom → length
    chrom_lens = {name: ref_fasta.get_reference_length(name) for name in ref_fasta.references}
    ploidy = 1 if haploid else 2
    n_samples = None
    with _open_vcf(vcf_path) as fh:
        for line in fh:
            if line.startswith("##"):
                continue
            if line.startswith("#"):
                cols = line.rstrip("\n").split("\t")
                n_samples = len(cols) - 9
                if n_samples < 1:
                    raise RuntimeError("VCF has no sample columns")
                continue
            cols = line.rstrip("\n").split("\t")
            chrom = cols[0]
            pos_1based = int(cols[1])
            ref = cols[3].upper()
            alts_field = cols[4]
            # 0-based start (PG semantics)
            start_0 = pos_1based - 1
            end_0 = start_0 + len(ref)

            # validate ALT characters
            valid = True
            for a in alts_field.split(","):
                for c in a.upper():
                    if c not in "ACGT,":
                        valid = False
                        break
                if not valid:
                    break
            if not valid:
                # PG prints a warning here; we mirror that behavior silently
                print(f"  skip variant at {chrom}:{start_0} (non-ACGT ALT)", file=sys.stderr)
                continue

            # chrom-end check (≥2k from start AND ≥2k from end)
            chrom_len = chrom_lens.get(chrom)
            if chrom_len is None:
                continue
            if start_0 < 2 * k or end_0 + 2 * k > chrom_len:
                print(f"  skip variant at {chrom}:{start_0} (too close to chrom end)", file=sys.stderr)
                continue

            # validate REF matches reference. PG itself does not validate REF↔FA;
            # we do, but we tolerate non-ACGT positions (IUPAC codes in FA, or N
            # in VCF where cactus has normalized IUPAC → N). A mismatch in ACGT
            # vs ACGT is still a hard error (real coordinate-system bug).
            ref_obs = ref_fasta.fetch(chrom, start_0, end_0).upper()
            if ref_obs != ref:
                acgt = set("ACGT")
                strict_mismatch = any(
                    (a in acgt) and (b in acgt) and (a != b)
                    for a, b in zip(ref, ref_obs)
                )
                if strict_mismatch or len(ref) != len(ref_obs):
                    raise RuntimeError(
                        f"REF mismatch at {chrom}:{pos_1based}: VCF={ref!r} FA={ref_obs!r}"
                    )

            alts = alts_field.split(",")
            # parse GT field → per-path allele indices
            paths = []
            if add_reference:
                paths.append(0)  # PG adds reference as the first path
            for i in range(9, 9 + n_samples):
                gt = cols[i]
                if haploid:
                    # haploid: gt is a single integer or '.'
                    a = _parse_haploid_allele(gt, alts)
                    paths.append(a)
                else:
                    if "/" in gt:
                        raise RuntimeError(f"Unphased GT at {chrom}:{pos_1based}: {gt}")
                    parts = gt.split("|")
                    if len(parts) != 2:
                        raise RuntimeError(f"GT not diploid at {chrom}:{pos_1based}: {gt}")
                    paths.append(_parse_haploid_allele(parts[0], alts))
                    paths.append(_parse_haploid_allele(parts[1], alts))

            yield Variant(chrom, start_0, ref, alts, paths)


def _parse_haploid_allele(s: str, alts: list[str]) -> int:
    """Return allele index. Missing '.' returns -1 (undefined).
    Caller treats -1 as 'this path passes through an N-allele' (PG behavior:
    adds an N-allele to the variant; we skip such paths for unique-k-mer
    computation since N-containing k-mers are filtered by enumerate_fwd_kmers)."""
    if s == ".":
        return -1
    return int(s)


def cluster_into_bubbles(variants: Iterable[Variant], k: int) -> list[list[Variant]]:
    """Group consecutive Variants whose distance < k-1 into the same bubble.
    Variants must be sorted by (chrom, start). Returns list of bubbles
    (each = list of Variants on the same chrom)."""
    bubbles: list[list[Variant]] = []
    cur: list[Variant] = []
    prev_chrom = None
    prev_end = 0
    for v in variants:
        if (prev_chrom is None) or (prev_chrom != v.chrom) or (v.start - prev_end >= k - 1):
            if cur:
                bubbles.append(cur)
            cur = [v]
        else:
            cur.append(v)
        prev_chrom = v.chrom
        prev_end = v.end
    if cur:
        bubbles.append(cur)
    return bubbles


def enumerate_bubble_alleles(bubble: list[Variant], add_reference: bool):
    """Given a bubble (list of variants whose path-vectors are aligned),
    return:
        unique_alleles: list of tuples (allele_idx_per_variant)
                       Always includes (0,0,...,0) as allele 0 (REF).
        path_to_allele: list[int] mapping each path index → unique-allele index.
    """
    n_paths = len(bubble[0].paths)
    # Sanity check: all variants in bubble must agree on number of paths
    for v in bubble:
        if len(v.paths) != n_paths:
            raise RuntimeError("path-count mismatch within bubble")

    # Build per-path allele combination tuple
    path_to_combo: list[tuple[int, ...]] = []
    for p in range(n_paths):
        combo = tuple(v.paths[p] for v in bubble)
        path_to_combo.append(combo)

    # PG includes (0,0,...,0) by default, even if no path actually has it
    # (it sets path_to_index[ref_path] = {})
    seen = OrderedDict()
    if add_reference:
        seen[tuple(0 for _ in bubble)] = None
    for combo in path_to_combo:
        # Skip combos that touch a missing allele (-1). PG materializes these
        # as a synthetic N-allele which yields zero k-mers anyway; we drop them
        # at enumeration to avoid IndexError in allele_sequence(-1) and to save
        # the wasted FASTA write.
        if -1 in combo:
            continue
        if combo not in seen:
            seen[combo] = None
    unique_alleles = list(seen.keys())

    # path_to_allele kept for parity with prior signature but no longer indexed
    # for missing combos (returns -1 for those). Caller discards it.
    allele_index = {c: i for i, c in enumerate(unique_alleles)}
    path_to_allele = [allele_index.get(c, -1) for c in path_to_combo]
    return unique_alleles, path_to_allele


def build_allele_sequence(
    bubble: list[Variant],
    allele_combo: tuple[int, ...],
    ref_fasta: pysam.FastaFile,
    k: int,
) -> str:
    """Construct the allele sequence: left_flank + variant_seq_1 + inner_flank_1 + ... + right_flank."""
    chrom = bubble[0].chrom
    bubble_start = bubble[0].start
    bubble_end = bubble[-1].end
    # left flank: k-1 bases immediately before bubble start
    left_flank = ref_fasta.fetch(chrom, bubble_start - (k - 1), bubble_start).upper()
    # right flank: k-1 bases immediately after bubble end
    right_flank = ref_fasta.fetch(chrom, bubble_end, bubble_end + (k - 1)).upper()

    parts = [left_flank]
    for i, v in enumerate(bubble):
        parts.append(v.allele_sequence(allele_combo[i]))
        if i < len(bubble) - 1:
            # inner flank between bubble[i].end and bubble[i+1].start
            inner = ref_fasta.fetch(chrom, v.end, bubble[i + 1].start).upper()
            parts.append(inner)
    parts.append(right_flank)
    return "".join(parts)


# ---------------------------------------------------------------------------
# Per-bubble unique-k-mer selection
# ---------------------------------------------------------------------------

def select_unique_kmers(
    bubble_alleles: list[str],
    genomic_kmer_count,  # callable: bytes -> int
    is_biallelic: bool,
    cap_biallelic: int,
    cap_multiallelic: int,
    k: int,
) -> dict[int, list[bytes]]:
    """Mirror PG's stepwise_unique_kmers + select_kmers logic.

    Returns: { allele_idx → ordered list of canonical k-mers chosen }
    """
    # Step 1: for each allele, count FORWARD k-mer occurrences (within-allele).
    # PG enumerates k-mers via jellyfish mer_dna shift_left which does NOT
    # canonicalize. Canonicalization happens only when looking up the genomic
    # count (jellyfish was built with -C). We mirror that exactly.
    occurrences: dict[bytes, list[int]] = defaultdict(list)
    for a_idx, seq in enumerate(bubble_alleles):
        counts: dict[bytes, int] = defaultdict(int)
        for km in enumerate_fwd_kmers(seq, k):
            counts[km] += 1
        # kmer is "unique to this allele" if it appears exactly once in counts;
        # PG keeps k-mers that occur EXACTLY ONCE in the allele AND appear in
        # exactly one allele across the bubble.
        for km, c in counts.items():
            if c == 1:
                occurrences[km].append(a_idx)

    # Step 2: filter to k-mers that (genomic_count == local_count) (i.e. not
    # elsewhere in the panel) AND occur in exactly one allele of this bubble.
    # PG iterates the occurrences map in canonical-k-mer-sorted order
    # (std::map<jellyfish::mer_dna, ...>); we mirror that by sorting keys
    # lexicographically before populating allele_to_kmers, so the per-allele
    # cap picks the same subset of k-mers.
    allele_to_kmers: dict[int, list[bytes]] = defaultdict(list)
    for km in sorted(occurrences.keys()):
        alleles_carrying = occurrences[km]
        if len(alleles_carrying) != 1:
            continue   # appears in multiple alleles → skip
        local_count = 1
        # query jellyfish on the CANONICAL form (jellyfish was built with -C)
        genomic = genomic_kmer_count(canonical(km))
        if genomic - local_count != 0:
            continue   # k-mer occurs elsewhere → skip
        allele_to_kmers[alleles_carrying[0]].append(km)

    # Step 3: round-robin selection up to per-allele cap, with per-bubble total
    n_alleles = len(bubble_alleles)
    max_kmers = cap_biallelic if is_biallelic else cap_multiallelic
    # Per-bubble total cap: max(nr_paths, 301). nr_paths ≈ n_alleles (we don't
    # have the original path count here; use n_alleles which is conservative)
    max_alleles_total = max(n_alleles, 301)

    result: dict[int, list[bytes]] = defaultdict(list)
    nr_selected = 0
    # PG iterates allele_to_kmers in std::map (canonical k-mer-sorted order).
    # We can't reproduce that exactly without jellyfish; for now use insertion
    # order by allele_idx. We'll address strict ordering in validation phase.
    sorted_alleles = sorted(allele_to_kmers.keys())
    keep_adding = True
    while keep_adding and nr_selected < max_alleles_total:
        kmer_added = False
        for a in sorted_alleles:
            if allele_to_kmers[a] and len(result[a]) < max_kmers:
                result[a].append(allele_to_kmers[a].pop(0))
                kmer_added = True
                nr_selected += 1
            if nr_selected >= max_alleles_total:
                break
        keep_adding = kmer_added
    return result


# ---------------------------------------------------------------------------
# Genomic k-mer count via jellyfish
# ---------------------------------------------------------------------------

class JellyfishCounter:
    """Wraps `jellyfish count` + on-disk hash. Uses dna_jellyfish Python bindings
    (SWIG wrappers around libjellyfish) for direct C++-level queries against the
    on-disk hash. Constant memory; no Python dict copy of all genomic k-mers.

    Requires the `dna_jellyfish` package (ships with bioconda `kmer-jellyfish`).
    """

    def __init__(self, fasta_path: str, k: int, hash_size: int = 100_000_000, threads: int = 4, workdir: str | None = None):
        import dna_jellyfish as _jf  # local import so the module loads without the
                                     # binding; it is REQUIRED here at construction
                                     # (no fallback) and raises if missing.
        self._jf = _jf
        self.k = k
        self.workdir = workdir or tempfile.mkdtemp(prefix="jf_")
        self.jf_path = os.path.join(self.workdir, "panel.jf")
        self._count(fasta_path, hash_size, threads)
        # configure mer_dna global k-size; required before any MerDNA op
        self._jf.MerDNA.k(self.k)
        # QueryMerFile: opens the hash file with mmap, queries in C++ time.
        self.qf = self._jf.QueryMerFile(self.jf_path)
        # Reusable MerDNA scratch object — avoids per-query allocation in the
        # inner loop.
        self._scratch = self._jf.MerDNA()

    def _count(self, fasta_path, hash_size, threads):
        cmd = [
            "jellyfish", "count",
            "-m", str(self.k),
            "-s", str(hash_size),
            "-t", str(threads),
            "-C",  # canonical
            "-o", self.jf_path,
            fasta_path,
        ]
        print(f"  [jellyfish] {' '.join(cmd)}", file=sys.stderr)
        subprocess.run(cmd, check=True)

    def get_count(self, kmer: bytes) -> int:
        """Look up genomic count for kmer. `kmer` is a fwd-direction bytes object;
        QueryMerFile auto-canonicalizes internally.

        Fast path: reuse a scratch MerDNA via .set(str) — avoids per-call
        allocation of MerDNA + string conversion.
        """
        self._scratch.set(kmer.decode("ascii") if isinstance(kmer, (bytes, bytearray)) else kmer)
        return self.qf[self._scratch]


# ---------------------------------------------------------------------------
# Overhang flanking-kmer selection
# ---------------------------------------------------------------------------

def select_overhang_kmers(
    bubble_chrom: str,
    bubble_start: int,
    bubble_end: int,
    prev_bubble_end: int,    # 0 if no prev bubble on this chrom
    next_bubble_start: int,  # chrom length if no next bubble on this chrom
    ref_fasta: pysam.FastaFile,
    k: int,
    overhang_size: int,
    overhang_cap: int,
    genomic_kmer_count,
) -> list[bytes]:
    """For each side (left, right) flank, select up to `overhang_cap` k-mers
    unique genome-wide. PG clips the overhang to the adjacent bubble's
    boundary so the same ref region is never counted in two bubbles."""
    # left: max(bubble_start - overhang_size, prev_bubble_end) .. bubble_start
    left_start = max(bubble_start - overhang_size, prev_bubble_end)
    left_seq = ref_fasta.fetch(bubble_chrom, left_start, bubble_start).upper()
    # right: bubble_end .. min(bubble_end + overhang_size, next_bubble_start)
    right_end = min(bubble_end + overhang_size, next_bubble_start)
    right_seq = ref_fasta.fetch(bubble_chrom, bubble_end, right_end).upper()

    out: list[bytes] = []
    for seq in (left_seq, right_seq):
        # PG enumerates fwd k-mers; canonicalize only at jellyfish lookup time.
        counts: dict[bytes, int] = defaultdict(int)
        for km in enumerate_fwd_kmers(seq, k):
            counts[km] += 1
        selected = 0
        for km in sorted(counts.keys()):   # PG uses std::map iteration order
            if selected >= overhang_cap:
                break
            if counts[km] != 1:
                continue
            if genomic_kmer_count(canonical(km)) != 1:
                continue
            out.append(km)
            selected += 1
    return out


# ---------------------------------------------------------------------------
# Orchestration
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("--vcf", required=True, help="Input VCF (uncompressed)")
    ap.add_argument("--ref", required=True, help="Reference FASTA")
    ap.add_argument("--out", required=True, help="Output prefix")
    ap.add_argument("-k", "--kmer-size", type=int, default=31)
    ap.add_argument("--haploid", action="store_true", help="Accept haploid GTs (e.g. cactus output)")
    ap.add_argument("--no-add-reference", action="store_true", help="Omit synthetic REF path (PG includes by default)")
    ap.add_argument("--cap-biallelic", type=int, default=16)
    ap.add_argument("--cap-multiallelic", type=int, default=32)
    ap.add_argument("--overhang-cap", type=int, default=12)
    ap.add_argument("--jellyfish-threads", type=int, default=4)
    ap.add_argument("--jellyfish-hash", type=int, default=100_000_000)
    ap.add_argument("--keep-tempfiles", action="store_true")
    args = ap.parse_args()

    add_reference = not args.no_add_reference
    k = args.kmer_size
    overhang_size = 2 * k

    print(f"[init] vcf={args.vcf} ref={args.ref} k={k} haploid={args.haploid} add_ref={add_reference}", file=sys.stderr)

    ref_fasta = pysam.FastaFile(args.ref)

    # 1. parse + cluster
    print("[stage 1] parsing VCF + building bubbles", file=sys.stderr)
    variants = list(parse_vcf_variants(args.vcf, ref_fasta, k, args.haploid, add_reference))
    print(f"  {len(variants)} variants accepted", file=sys.stderr)
    bubbles = cluster_into_bubbles(variants, k)
    print(f"  {len(bubbles)} bubbles", file=sys.stderr)

    # 2. enumerate alleles + write path_segments.fasta (input to jellyfish)
    print("[stage 2] enumerating alleles + writing path_segments FASTA", file=sys.stderr)
    workdir = args.out + "_tmp"
    os.makedirs(workdir, exist_ok=True)
    segments_fasta = os.path.join(workdir, "path_segments.fasta")
    # Pre-compute all bubble alleles so we can stream once for jellyfish
    # and once again for unique-kmer selection.
    bubble_alleles_all: list[list[str]] = []     # per bubble: list of allele sequences
    bubble_unique_combos_all: list[list[tuple]] = []  # per bubble: list of allele-combo tuples (in same order)

    with open(segments_fasta, "w") as fa:
        prev_end = {}  # per chrom
        for bidx, bubble in enumerate(bubbles):
            chrom = bubble[0].chrom
            bubble_start = bubble[0].start
            bubble_end = bubble[-1].end
            # write reference unitig from prev_end to bubble_start (inter-bubble)
            if chrom not in prev_end:
                prev_end[chrom] = 0
            ref_unitig = ref_fasta.fetch(chrom, prev_end[chrom], bubble_start).upper()
            fa.write(f">{chrom}_reference_{bubble_start}\n{ref_unitig}\n")
            unique_combos, _ = enumerate_bubble_alleles(bubble, add_reference)
            allele_seqs = []
            for a_idx, combo in enumerate(unique_combos):
                seq = build_allele_sequence(bubble, combo, ref_fasta, k)
                allele_seqs.append(seq)
                fa.write(f">{chrom}_{bubble_start}_{a_idx}\n{seq}\n")
            bubble_alleles_all.append(allele_seqs)
            bubble_unique_combos_all.append(unique_combos)
            prev_end[chrom] = bubble_end
        # tail unitigs per chrom
        for c, end in prev_end.items():
            chrom_len = ref_fasta.get_reference_length(c)
            ref_tail = ref_fasta.fetch(c, end, chrom_len).upper()
            fa.write(f">{c}_reference_end\n{ref_tail}\n")
    print(f"  {segments_fasta}", file=sys.stderr)

    # 3. jellyfish count + open via dna_jellyfish bindings (constant memory)
    print("[stage 3] jellyfish count + open query handle", file=sys.stderr)
    jf = JellyfishCounter(segments_fasta, k, args.jellyfish_hash, args.jellyfish_threads, workdir)

    # 4. per-bubble unique kmer + overhang → emit TSV per chrom
    print("[stage 4] selecting unique kmers + writing TSV", file=sys.stderr)
    # one TSV per chrom (matching PG's layout)
    per_chrom_out: dict[str, gzip.GzipFile] = {}
    header = b"#chromosome\tstart\tend\tunique_kmers\tunique_kmers_overhang\n"

    # Precompute per-chrom bubble index for prev/next lookups in overhang clip
    chrom_bubble_indices: dict[str, list[int]] = defaultdict(list)
    for bidx, bubble in enumerate(bubbles):
        chrom_bubble_indices[bubble[0].chrom].append(bidx)

    for bidx, bubble in enumerate(bubbles):
        chrom = bubble[0].chrom
        out_fn = f"{args.out}_{chrom}_kmers.tsv.gz"
        if chrom not in per_chrom_out:
            f = gzip.open(out_fn, "wb")
            f.write(header)
            per_chrom_out[chrom] = f
        bstart = bubble[0].start
        bend = bubble[-1].end

        # Determine prev/next bubble bounds ON THIS CHROM for overhang clipping
        bidx_in_chrom = chrom_bubble_indices[chrom]
        pos_in_chrom = bidx_in_chrom.index(bidx)
        if pos_in_chrom == 0:
            prev_bubble_end = 0
        else:
            prev_b = bubbles[bidx_in_chrom[pos_in_chrom - 1]]
            prev_bubble_end = prev_b[-1].end
        if pos_in_chrom == len(bidx_in_chrom) - 1:
            next_bubble_start = ref_fasta.get_reference_length(chrom)
        else:
            next_b = bubbles[bidx_in_chrom[pos_in_chrom + 1]]
            next_bubble_start = next_b[0].start

        unique_combos = bubble_unique_combos_all[bidx]
        is_biallelic = (len(unique_combos) == 2)
        allele_seqs = bubble_alleles_all[bidx]
        selected = select_unique_kmers(
            allele_seqs, jf.get_count, is_biallelic,
            args.cap_biallelic, args.cap_multiallelic, k
        )
        # flatten selection (per-allele list of k-mers) into a CSV
        # PG order: iterate allele_to_kmers ascending by allele idx; within
        # each allele, the order is queue-pop (which preserves jellyfish/map
        # canonical-k-mer order from `select_kmers`).
        flat = []
        for a in sorted(selected.keys()):
            flat.extend(selected[a])
        kmers_str = b",".join(flat).decode("ascii") if flat else "nan"

        overhang_kmers = select_overhang_kmers(
            chrom, bstart, bend, prev_bubble_end, next_bubble_start,
            ref_fasta, k, overhang_size, args.overhang_cap, jf.get_count
        )
        oh_str = b",".join(overhang_kmers).decode("ascii") if overhang_kmers else "nan"

        line = f"{chrom}\t{bstart}\t{bend}\t{kmers_str}\t{oh_str}\n"
        per_chrom_out[chrom].write(line.encode("ascii"))
    for f in per_chrom_out.values():
        f.close()

    # cleanup
    if not args.keep_tempfiles:
        import shutil
        shutil.rmtree(workdir, ignore_errors=True)

    print("[done]", file=sys.stderr)


if __name__ == "__main__":
    main()
