#!/usr/bin/env python
"""Per-block PANEL k-mer count for the fine CLQ0.9 blocks (Chr1), cross-referenced with
the benchmark's local-fit status, to find the k-mer-coverage threshold for a local h fit.

panel k-mers = # allele-specific k-mers in the filtered index whose bubble centroid falls
in the block (coverage-independent; assign_kmers_to_blocks scheme). status (from the base
window run) = 0 local-fit (>=50 OBSERVED at 10x) / 1 fallback / 2 empty.

Output: results/grenenet_gea/blocks_mcf90/chr1_block_kmer_coverage.csv
"""
import numpy as np, pandas as pd

META = "data/kmer_pa_231_arch3_filt2inv/kmer_pa_Chr1.meta.npz"
BLOCKS = "results/grenenet_gea/blocks_mcf90/chr1_clq0.9_blocks_clq0.9.tsv"
STATUS = "benchmarks/ldblock_window_test/sweep_s42_base.h_blocks_per_chrom.npz"


def main():
    b = pd.read_csv(BLOCKS, sep="\t")
    starts = b.start_pos.values; ends = b.end_pos.values
    print(f"{len(b)} fine blocks; loading k-mer index bubble coords ...", flush=True)
    m = np.load(META, allow_pickle=True)
    bid = m["bubble_id"]                      # per-kmer bubble index
    bstart = m["bubble_start"]; bend = m["bubble_end"]
    cent = (bstart + bend) // 2               # bubble centroid (all Chr1)
    # assign bubble -> block (blocks sorted, disjoint): searchsorted on starts
    i = np.clip(np.searchsorted(starts, cent, side="right") - 1, 0, len(starts) - 1)
    inb = (cent >= starts[i]) & (cent <= ends[i])
    bubble_block = np.where(inb, i, -1)
    kmer_block = bubble_block[bid]
    v = kmer_block >= 0
    panel = np.bincount(kmer_block[v], minlength=len(b))
    b["panel_kmers"] = panel
    print(f"  total panel k-mers {len(bid):,}; assigned to a block {v.sum():,} "
          f"({100*v.mean():.0f}%)", flush=True)

    st = np.load(STATUS, allow_pickle=True)["Chr1_status"]
    assert len(st) == len(b), (len(st), len(b))
    b["status"] = st
    b.to_csv("results/grenenet_gea/blocks_mcf90/chr1_block_kmer_coverage.csv", index=False)

    print("\npanel k-mers by benchmark status (0=local,1=fallback,2=empty):")
    print(b.groupby("status").panel_kmers.describe()[["count", "25%", "50%", "75%"]].round(0).to_string())
    print(f"\nblocks with panel_kmers==0 (true deserts): {int((b.panel_kmers==0).sum())} "
          f"({100*(b.panel_kmers==0).mean():.1f}%); these are status-2={int(((b.panel_kmers==0)&(b.status==2)).sum())}")
    print("\n% local-fit (status==0) by panel-kmer bin  ->  the coverage threshold:")
    bins = [0, 20, 50, 100, 150, 200, 300, 500, 1e9]
    b["pk_bin"] = pd.cut(b.panel_kmers, bins, right=False)
    g = b.groupby("pk_bin", observed=True)
    for iv, gg in g:
        print(f"  panel_kmers {str(iv):>14}: n={len(gg):6d}  %local {100*(gg.status==0).mean():5.0f}%  "
              f"med n_var {int(gg.n_variants.median())}")


if __name__ == "__main__":
    main()
