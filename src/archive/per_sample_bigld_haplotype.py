"""Per-sample driver for cactus_em --block-mode bigld_haplotype.

Drop-in replacement for per_sample_per_chrom.py for the new mode. Counts
k-mers, runs per-block haplotype-level EM, projects to per-record AFs.

Inputs:
  --block-hap-cn-prefix    e.g. data/block_haplotype_cn/chr1
                           expects <prefix>.npz or <prefix>_<chrom>.npz
  --cn-var, --cn-var-meta  same as per_sample_per_chrom.py
  --reads                  FASTQ pair or BAM

Output: per-record alt-freq TSV (same schema as per_sample_per_chrom.py).
"""
from __future__ import annotations
import argparse
import gc
import os
import sys
import time
from pathlib import Path

import numpy as np
from scipy.sparse import load_npz

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from block_haplotype_em import run_chrom_bigld_haplotype


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--block-hap-cn-prefix", required=True,
                    help="prefix for block_haplotype_cn npz; either"
                         " <prefix>.npz (single file) or <prefix>_<chrom>.npz"
                         " per chrom.")
    ap.add_argument("--cn-var", required=True)
    ap.add_argument("--cn-var-meta", required=True)
    ap.add_argument("--reads", required=True, nargs="+")
    ap.add_argument("--sample", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--chroms", nargs="+", default=["Chr1"])
    ap.add_argument("--em-max-iter", type=int, default=200)
    ap.add_argument("--min-block-kmer-count", type=int, default=30)
    ap.add_argument("--inference-unit", default="class",
                    choices=["class", "founder"],
                    help="class (default): EM on per-block sequence-unique "
                         "haplotype simplex, then scatter to founders by "
                         "founder_to_class. founder: scatter cn rows class->"
                         "founder before EM, run EM directly on the "
                         "231-founder simplex (more evidence per dim).")
    ap.add_argument("--save-h-blocks", default=None,
                    help="optional path to save per-block h_hap arrays")
    ap.add_argument("--smooth-passes", type=int, default=0,
                    help="Li-Stephens-style HMM smoothing across blocks; "
                         "0 = no smoothing (default)")
    ap.add_argument("--smooth-alpha", type=float, default=0.7,
                    help="weight on per-block evidence vs neighbors (0=full smooth, "
                         "1=no smooth)")
    ap.add_argument("--smooth-recomb-rate", type=float, default=4e-8,
                    help="A. thaliana recomb rate per bp (default 4e-8 = 4cM/Mb)")
    ap.add_argument("--dirichlet-alpha", type=float, default=1.0,
                    help="Dirichlet prior strength on h (1.0=MLE, 1.01-1.5 = mild "
                         "regularization to prevent winner-takes-all)")
    ap.add_argument("--fallback-to-global", action="store_true",
                    help="blend per-block h with chrom-global h for low-confidence blocks")
    ap.add_argument("--fallback-eff-n-threshold", type=float, default=0.0,
                    help="blend toward global if block's eff_n (=1/sum(h^2)) < threshold")
    args = ap.parse_args()

    print(f"=== {args.sample} (bigld_haplotype mode) ===", flush=True)
    print(f"  reads: {args.reads}", flush=True)
    t0 = time.time()

    cn_var = load_npz(args.cn_var)
    var_meta = np.load(args.cn_var_meta, allow_pickle=True)
    n_records = cn_var.shape[1]
    print(f"  cn_var: {cn_var.shape}, n_records: {n_records:,}", flush=True)

    reads_input = args.reads if len(args.reads) > 1 else args.reads[0]

    freqs_global = np.full(n_records, np.nan, dtype=np.float32)
    h_save = {}

    for chrom in args.chroms:
        # Resolve block_hap_cn npz path
        # Try <prefix>_<chrom>.npz first (per-chrom), fall back to <prefix>.npz
        candidates = [
            f"{args.block_hap_cn_prefix}_{chrom}.npz",
            f"{args.block_hap_cn_prefix}.npz",
            args.block_hap_cn_prefix,
        ]
        npz_path = next((p for p in candidates if os.path.exists(p)), None)
        if npz_path is None:
            print(f"  [{chrom}] no block_hap_cn npz found at "
                  f"{candidates}; skipping", flush=True)
            continue

        idx, freqs, pack, _ = run_chrom_bigld_haplotype(
            chrom=chrom,
            block_hap_npz_path=npz_path,
            cn_var=cn_var, var_meta=var_meta,
            reads_input=reads_input,
            threads=args.threads,
            em_max_iter=args.em_max_iter,
            min_block_kmer_count=args.min_block_kmer_count,
            inference_unit=args.inference_unit,
            smooth_passes=args.smooth_passes,
            smooth_alpha=args.smooth_alpha,
            smooth_recomb_rate=args.smooth_recomb_rate,
            dirichlet_alpha=args.dirichlet_alpha,
            fallback_to_global=args.fallback_to_global,
            fallback_eff_n_threshold=args.fallback_eff_n_threshold,
        )
        if idx is None:
            continue
        freqs_global[idx] = freqs
        if args.save_h_blocks:
            h_save[f"{chrom}_block_status"] = pack["block_status"]
            h_save[f"{chrom}_block_chrom"] = pack["block_chrom"]
            h_save[f"{chrom}_block_pos_start"] = pack["block_pos_start"]
            h_save[f"{chrom}_block_pos_end"] = pack["block_pos_end"]
            # Raw per-block h (class mode → n_uniq vector; founder mode → n_eco).
            # h_per_block_post is the founder-level h actually used for projection,
            # post fallback-blend + smoothing — what we want for identity diagnostics.
            h_arrays = np.empty(len(pack["h_per_block"]), dtype=object)
            for i, h in enumerate(pack["h_per_block"]):
                h_arrays[i] = h
            h_save[f"{chrom}_h_per_block"] = h_arrays
            if "h_per_block_post" in pack:
                h_post = np.empty(len(pack["h_per_block_post"]), dtype=object)
                for i, h in enumerate(pack["h_per_block_post"]):
                    h_post[i] = h
                h_save[f"{chrom}_h_per_block_post"] = h_post
        gc.collect()

    chrom_arr = np.asarray(var_meta["chrom"]).astype(str)
    pos_arr = np.asarray(var_meta["pos"]).astype(np.int64)
    ref_arr = np.asarray(var_meta["ref_len"]).astype(np.int64)
    alt_arr = np.asarray(var_meta["alt_len"]).astype(np.int64)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    print(f"writing TSV to {args.out}...", flush=True)
    with open(args.out, "w") as f:
        f.write("chrom\tpos\tref_len\talt_len\talt_freq\n")
        for i in range(n_records):
            af = freqs_global[i]
            af_str = "" if np.isnan(af) else f"{float(af):.6f}"
            f.write(f"{chrom_arr[i]}\t{pos_arr[i]}\t{ref_arr[i]}\t{alt_arr[i]}\t{af_str}\n")

    if args.save_h_blocks:
        np.savez_compressed(args.save_h_blocks, **h_save, allow_pickle=True)
        print(f"saved h_blocks diagnostics to {args.save_h_blocks}", flush=True)

    print(f"=== done. wall: {time.time()-t0:.0f}s ===", flush=True)


if __name__ == "__main__":
    main()
