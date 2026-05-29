"""
Build the founder × k-mer copy-number matrix from PanGenie-index output +
**per-founder consensus FASTAs** (path A).

Same panel-k-mer schema as `build_kmer_cn.py` (PanGenie's bubble + unique-kmer
list), but instead of reconstructing each founder's bubble sequence on the fly
from TAIR10 + VCF GTs, this builder reads the founder's bubble slice straight
out of their pre-built consensus FASTA. The FASTAs are produced upstream by
`bcftools consensus -H A -s <eco> -f TAIR10 founders_231_chr.haploid.vcf.gz`
(see panel/pangenie_genotyping/scripts/build_consensus_fastas_one.sh).

Why path A: making the founder sequences explicit on disk makes the build
deterministic and inspectable — no compound-variant edge cases in custom code,
no dependence on bcftools consensus + per-bubble VCF queries running in lock
step. The output schema is identical to build_kmer_cn.py so this is a drop-in
replacement for cn_full_231_v2 → cn_full_231_v3.

Input layout:
  <fastas-dir>/<eco>.chr.fa   (one per founder; chrom-keyed)
  <fastas-dir>/<eco>.chr.fa.fai

Each .chr.fa must hold contigs named exactly as `chrom` ("Chr1" etc.) — the
output of `bcftools consensus` with the production TAIR10.chr.numeric.iupacN.fa.

Output (matches build_kmer_cn.py):
  <out>.cn.npz       sparse CSR (F × K_total) int8
  <out>.meta.npz     kmer_index, bubble_id, bubble_chrom/start/end, founders
"""
from __future__ import annotations
import argparse, gzip, sys, time
from pathlib import Path
import numpy as np
import pysam
from scipy.sparse import csr_matrix, save_npz


_COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(s: str) -> str:
    return s.translate(_COMPLEMENT)[::-1]


def canonical(kmer: str) -> str:
    rc = revcomp(kmer)
    return kmer if kmer < rc else rc


def canonical_kmer_set(seq: str, k: int) -> set:
    s = set()
    seq = seq.upper()
    for i in range(len(seq) - k + 1):
        km = seq[i : i + k]
        if "N" in km:
            continue
        s.add(canonical(km))
    return s


def parse_kmer_file(kmers_tsv_gz: str):
    bubbles = []
    with gzip.open(kmers_tsv_gz, "rt") as f:
        next(f)  # header
        for line in f:
            p = line.rstrip("\n").split("\t")
            chrom, start, end = p[0], int(p[1]), int(p[2])
            kmers = p[3].split(",") if p[3] else []
            overhang = p[4] if len(p) > 4 else ""
            bubbles.append({
                "chrom": chrom, "start": start, "end": end,
                "kmers": kmers, "overhang": overhang,
            })
    return bubbles


def list_founders(fastas_dir: Path, founders_order: list[str] | None = None) -> list[str]:
    """Return ordered founder list. If founders_order is given, validate and use it.
    Else list *.chr.fa in fastas_dir lexically."""
    if founders_order is not None:
        for f in founders_order:
            p = fastas_dir / f"{f}.chr.fa"
            if not p.exists():
                sys.exit(f"[fatal] expected FASTA missing: {p}")
        return founders_order
    return sorted(p.stem.replace(".chr", "") for p in fastas_dir.glob("*.chr.fa"))


def build_cn_for_chrom(
    kmers_tsv_gz: str,
    fastas_dir: str,
    chrom: str,
    founders_order: list[str] | None = None,
    k: int = 31,
    flank: int = 100,
    max_bubbles: int | None = None,
    verbose: bool = True,
):
    fastas_dir_p = Path(fastas_dir)
    bubbles = parse_kmer_file(kmers_tsv_gz)
    bubbles = [b for b in bubbles if b["chrom"] == chrom]
    if max_bubbles is not None:
        bubbles = bubbles[:max_bubbles]
    if verbose:
        print(f"[build_cn] {chrom}: {len(bubbles):,} bubbles", flush=True)

    founders = list_founders(fastas_dir_p, founders_order)
    F = len(founders)
    if verbose:
        print(f"[build_cn] founders: F={F}", flush=True)

    # Open every founder's FASTA. pysam.FastaFile lazy-mmaps; F=231 is fine.
    fasta_handles = {f: pysam.FastaFile(str(fastas_dir_p / f"{f}.chr.fa")) for f in founders}

    # Pre-compute canonical k-mers for each bubble + total K
    bubble_kmer_canon = []
    total_K = 0
    for b in bubbles:
        canon = [canonical(km) for km in b["kmers"]]
        bubble_kmer_canon.append(canon)
        total_K += len(canon)
    if verbose:
        print(f"[build_cn] total k-mers across bubbles: {total_K:,}", flush=True)

    cn_rows: list[int] = []
    cn_cols: list[int] = []
    bubble_id = np.zeros(total_K, dtype=np.int32)
    kmer_index = [""] * total_K

    t0 = time.time()
    k_offset = 0
    for b_idx, b in enumerate(bubbles):
        canon = bubble_kmer_canon[b_idx]
        K_b = len(canon)
        if K_b == 0:
            continue

        # Region with flanks. bcftools consensus output is 1-based and has
        # the same chrom names as TAIR10. pysam.fetch is 0-based half-open.
        rs = max(0, b["start"] - 1 - flank)
        re_ = b["end"] + flank

        for j, km in enumerate(canon):
            bubble_id[k_offset + j] = b_idx
            kmer_index[k_offset + j] = km

        for f_idx, founder in enumerate(founders):
            fa = fasta_handles[founder]
            try:
                hap = fa.fetch(chrom, rs, re_).upper()
            except Exception as e:
                if verbose and b_idx < 3:
                    print(f"  WARN: {founder} fetch failed at {chrom}:{rs}-{re_}: {e}",
                          flush=True)
                continue
            hap_kmers = canonical_kmer_set(hap, k)
            for j, km in enumerate(canon):
                if km in hap_kmers:
                    cn_rows.append(f_idx)
                    cn_cols.append(k_offset + j)

        k_offset += K_b
        if verbose and (b_idx + 1) % 1000 == 0:
            print(f"  ... {b_idx+1:,}/{len(bubbles):,} bubbles, "
                  f"nnz={len(cn_rows):,}, elapsed={time.time()-t0:.0f}s",
                  flush=True)

    cn = csr_matrix(
        (np.ones(len(cn_rows), dtype=np.int8), (cn_rows, cn_cols)),
        shape=(F, total_K),
        dtype=np.int8,
    )
    bubble_meta = [(b["chrom"], b["start"], b["end"]) for b in bubbles]
    if verbose:
        print(f"[build_cn] DONE. cn shape={cn.shape}, nnz={cn.nnz:,}, "
              f"density={cn.nnz / (cn.shape[0] * cn.shape[1]):.4%}",
              flush=True)

    for h in fasta_handles.values():
        h.close()

    return cn, kmer_index, bubble_id, bubble_meta, founders


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--kmers", required=True, help="PanGenie kmers.tsv.gz for one chrom")
    ap.add_argument("--fastas-dir", required=True,
                    help="Dir of per-founder FASTAs (<eco>.chr.fa each)")
    ap.add_argument("--chrom", required=True, help="e.g. Chr1")
    ap.add_argument("--out", required=True, help="Output prefix")
    ap.add_argument("--founders-file",
                    help="Optional: file with one founder name per line "
                         "to fix order. Defaults to lex-sorted FASTA basenames")
    ap.add_argument("--max-bubbles", type=int, default=None)
    ap.add_argument("--k", type=int, default=31)
    ap.add_argument("--flank", type=int, default=100)
    args = ap.parse_args()

    founders_order = None
    if args.founders_file:
        with open(args.founders_file) as f:
            founders_order = [ln.strip() for ln in f if ln.strip()]

    cn, kmer_index, bubble_id, bubble_meta, founders = build_cn_for_chrom(
        args.kmers, args.fastas_dir, args.chrom,
        founders_order=founders_order,
        k=args.k, flank=args.flank, max_bubbles=args.max_bubbles,
    )
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    save_npz(args.out + ".cn.npz", cn)
    np.savez(
        args.out + ".meta.npz",
        kmer_index=np.array(kmer_index),
        bubble_id=bubble_id,
        bubble_chrom=np.array([m[0] for m in bubble_meta]),
        bubble_start=np.array([m[1] for m in bubble_meta], dtype=np.int64),
        bubble_end=np.array([m[2] for m in bubble_meta], dtype=np.int64),
        founders=np.array(founders),
    )
    print(f"Wrote {args.out}.cn.npz and {args.out}.meta.npz")


if __name__ == "__main__":
    main()
