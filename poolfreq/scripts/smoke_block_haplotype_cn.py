"""Smoke test for haplotype-level cn matrix on a SINGLE BigLD block.

Goal: validate the haplotype-level k-mer matrix idea on one block before
scaling up. Picks a SNP-poor block (worst case for bigld_panel) and
reports:
  - n_uniq_haps in the block
  - n_kmers from full-haplotype reconstruction (per-rep)
  - n_kmers after within-block discriminative filter
  - sparsity of cn_block
  - average #carriers per kept k-mer

If n_kept_kmers >> n_uniq_haps (over-determined), the haplotype-level EM
should be statistically viable.

Usage:
    python smoke_block_haplotype_cn.py --block-idx 50  # try one block
"""
from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np
import pysam


def canonical_kmers(seq: str, k: int = 31) -> set[str]:
    """Return the set of canonical k-mers (lex-min of fwd/rc) in seq.

    Skips k-mers containing N or other ambiguous bases.
    """
    seq = seq.upper()
    rc_table = str.maketrans("ACGT", "TGCA")
    out = set()
    for i in range(len(seq) - k + 1):
        kmer = seq[i:i + k]
        if any(c not in "ACGT" for c in kmer):
            continue
        rc = kmer.translate(rc_table)[::-1]
        out.add(min(kmer, rc))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--block-index",
                    default="sims/visor_freqk/chr1_only_panel/hapfire_block_index_chr1.npz")
    ap.add_argument("--fastas-dir",
                    default="sims/visor_freqk/founder_fastas_231")
    ap.add_argument("--block-idx", type=int, default=None,
                    help="block index to test; if omitted picks a snp-poor block")
    ap.add_argument("--k", type=int, default=31)
    args = ap.parse_args()

    z = np.load(args.block_index, allow_pickle=True)
    ecotypes = np.asarray(z["ecotypes"]).astype(str)
    block_chrom = np.asarray(z["block_chrom"]).astype(str)
    block_pos_start = np.asarray(z["block_pos_start"]).astype(np.int64)
    block_pos_end = np.asarray(z["block_pos_end"]).astype(np.int64)
    block_n_uniq = np.asarray(z["block_n_uniq"]).astype(np.int32)
    block_n_snps = np.asarray(z["block_n_snps"]).astype(np.int32)
    hap_idx_d1 = np.asarray(z["hap_idx_d1"]).astype(np.int32)
    # selfing -> d1 == d2 in panel; we use d1 throughout for haploid-mode

    if args.block_idx is None:
        # Pick a block in the middle of the SNP-distribution (median snp count)
        # but with reasonable size — gives us a representative case
        bp_size = block_pos_end - block_pos_start
        ok = (bp_size >= 5_000) & (bp_size <= 20_000)
        snp_med = int(np.median(block_n_snps[ok]))
        cand = np.where(ok & (block_n_snps <= snp_med))[0]
        # Pick the one closest to median size
        sizes = bp_size[cand]
        med_sz = int(np.median(sizes))
        b = cand[int(np.argmin(np.abs(sizes - med_sz)))]
    else:
        b = args.block_idx
    chrom_short = block_chrom[b]
    chrom = f"Chr{chrom_short}" if not chrom_short.startswith("Chr") else chrom_short
    pos_start = int(block_pos_start[b])
    pos_end = int(block_pos_end[b])
    n_uniq = int(block_n_uniq[b])
    n_snps = int(block_n_snps[b])
    bp_size = pos_end - pos_start + 1

    print(f"Block {b}: {chrom}:{pos_start:,}-{pos_end:,} "
          f"({bp_size:,} bp), n_snps={n_snps}, n_uniq_haps={n_uniq}")

    # For each unique-hap class, pick a representative founder (the lowest-
    # indexed founder mapping to that class — deterministic).
    hap_d1 = hap_idx_d1[b]  # n_eco
    rep_founder_for_hap = -np.ones(n_uniq, dtype=np.int32)
    for f in range(len(hap_d1)):
        h = hap_d1[f]
        if rep_founder_for_hap[h] == -1:
            rep_founder_for_hap[h] = f
    print(f"  representatives picked for all {n_uniq} unique haps")

    # Slice the block from each rep's FASTA, hash canonical k-mers, build
    # hap × kmer presence-set
    fdir = Path(args.fastas_dir)
    hap_kmers: list[set[str]] = []
    for h in range(n_uniq):
        f = int(rep_founder_for_hap[h])
        eco = ecotypes[f]
        fa_path = fdir / f"{eco}.chr.fa"
        fa = pysam.FastaFile(str(fa_path))
        # 1-based pos_start, pysam fetch is 0-based half-open
        seq = fa.fetch(chrom, pos_start - 1, pos_end)
        fa.close()
        kmers = canonical_kmers(seq, k=args.k)
        hap_kmers.append(kmers)
    sizes = [len(s) for s in hap_kmers]
    print(f"  per-hap n_kmers: median={int(np.median(sizes))}, "
          f"min={min(sizes)}, max={max(sizes)} (block bp ~{bp_size:,})")

    # Build the discriminative kmer set: k-mers carried by 1..(n_uniq-1) haps
    union = set()
    for s in hap_kmers:
        union |= s
    print(f"  union k-mers across haps: {len(union):,}")

    carrier_count = np.zeros(len(union), dtype=np.int32)
    kmer_list = sorted(union)
    kmer_to_idx = {km: i for i, km in enumerate(kmer_list)}
    cn = np.zeros((n_uniq, len(union)), dtype=np.uint8)
    for h, s in enumerate(hap_kmers):
        for km in s:
            j = kmer_to_idx[km]
            cn[h, j] = 1
            carrier_count[j] += 1
    discrim_mask = (carrier_count > 0) & (carrier_count < n_uniq)
    n_discrim = int(discrim_mask.sum())
    cn_d = cn[:, discrim_mask]

    print(f"  discriminative k-mers (1<=carriers<={n_uniq-1}): "
          f"{n_discrim:,} ({n_discrim/max(1,len(union))*100:.1f}% of union)")
    if n_discrim > 0:
        cc = carrier_count[discrim_mask]
        print(f"    carrier counts: median={int(np.median(cc))}, "
              f"min={cc.min()}, max={cc.max()}")
        print(f"    mean k-mers/hap (in discrim set): "
              f"{cn_d.sum(axis=1).mean():.1f}")

    # Identifiability check: K_discrim vs n_uniq
    print(f"  EM identifiability: K_discrim={n_discrim} vs n_uniq={n_uniq} "
          f"(over-determined factor: {n_discrim/max(1,n_uniq):.1f}×)")

    # Compare with what bigld_panel currently has: cn_full_231 k-mers within
    # this block range (only PanGenie bubble-local-unique k-mers)
    try:
        cn_full_meta = np.load(
            "poolfreq/data/cn_full_231_v2/cn_Chr1.meta.npz",
            allow_pickle=True)
        bub_chrom = np.asarray(cn_full_meta["bubble_chrom"]).astype(str)
        bub_start = np.asarray(cn_full_meta["bubble_start"]).astype(np.int64)
        # bubble_id might give the kmer-to-bubble mapping
        # We want k-mers whose bubble_pos overlaps the block window
        in_block = (bub_chrom == chrom) & (bub_start >= pos_start) & (bub_start <= pos_end)
        n_pangenie_kmers = int(in_block.sum())
        print(f"  bigld_panel comparison: PanGenie cn_full has "
              f"{n_pangenie_kmers:,} k-mers in this block "
              f"(vs {n_discrim:,} for haplotype-level)")
    except Exception as e:
        print(f"  (could not compare with cn_full_231: {e})")


if __name__ == "__main__":
    main()
