"""
Multi-sample batch runner: process N pool-seq samples in parallel.

Each sample produces a per-VCF-record alt-allele frequency TSV. Output paths
follow a deterministic scheme so the runner is restartable (skip samples whose
output already exists).

Usage:
    python batch_runner.py --manifest samples.tsv --out-dir results/ \
        --cn-kmer-prefix data/cn_full --cn-var data/cn_var_82.cn_var.npz \
        --cn-var-meta data/cn_var_82.meta.npz --threads 8

Manifest format (TSV with header):
    sample_id  reads_path                    [reads_path2]
    MLFH001    /path/to/MLFH001.fq.gz        /path/to/MLFH001_R2.fq.gz
    MLFH002    /path/to/MLFH002.fq.gz        /path/to/MLFH002_R2.fq.gz
    ...

For the GrENE-Net 2,415 evolved samples, the manifest is constructed from
freqk_gr/data/samples.tsv.
"""
from __future__ import annotations
import argparse, os, sys, time, multiprocessing as mp
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import pandas as pd
from scipy.sparse import load_npz
from per_sample_driver import load_cn_kmer_genomewide, run_one_sample


def process_one(args):
    """Worker: process a single sample. Returns (sample_id, status, runtime)."""
    (sample_id, reads_paths, cn_kmer_data, cn_var, var_meta, bubble_meta,
     founders, out_path, threads, block_mode, window_bp) = args
    t0 = time.time()
    if os.path.exists(out_path):
        return (sample_id, "skipped (exists)", 0)
    try:
        cn_kmer, kmer_index = cn_kmer_data
        run_one_sample(
            cn_kmer, kmer_index, cn_var, var_meta, founders,
            reads_paths, sample_id, out_path,
            em_max_iter=200, threads=threads,
            block_mode=block_mode, window_bp=window_bp,
            bubble_meta=bubble_meta,
        )
        return (sample_id, "ok", time.time() - t0)
    except Exception as e:
        return (sample_id, f"FAILED: {e}", time.time() - t0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", required=True,
                    help="TSV with header: sample_id, reads_path, [reads_path2]")
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--cn-kmer-prefix", required=True,
                    help="prefix for per-chrom cn matrices (e.g., data/cn_full)")
    ap.add_argument("--cn-var", required=True)
    ap.add_argument("--cn-var-meta", required=True)
    ap.add_argument("--threads", type=int, default=4,
                    help="threads for jellyfish k-mer counting")
    ap.add_argument("--workers", type=int, default=1,
                    help="parallel samples (each uses --threads). Sample-level "
                         "parallelism is what actually scales; block-level "
                         "threading inside one sample is BW-bound and doesn't "
                         "help — prefer larger SLURM arrays over multi-worker.")
    ap.add_argument("--limit", type=int, default=None,
                    help="optional: process only first N samples (for testing)")
    ap.add_argument("--block-mode", default="window",
                    choices=["window", "global"],
                    help="window: HAFpipe-style per-window EM (default). "
                         "global: single genome-wide h.")
    ap.add_argument("--window-bp", type=int, default=200_000,
                    help="window size for block_mode=window (default 200 kb)")
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    # Load shared data once
    print(f"[{time.strftime('%H:%M:%S')}] Loading shared cn matrices...")
    (cn_kmer, kmer_index, bubble_id, bubble_chrom,
     bubble_start, bubble_end, founders) = load_cn_kmer_genomewide(
        args.cn_kmer_prefix)
    cn_var = load_npz(args.cn_var)
    meta_raw = np.load(args.cn_var_meta, allow_pickle=True)
    var_meta = {k: meta_raw[k] for k in ["chrom", "pos", "ref_len", "alt_len"]}
    bubble_meta = {
        "bubble_id": bubble_id,
        "bubble_chrom": bubble_chrom,
        "bubble_start": bubble_start,
        "bubble_end": bubble_end,
    }
    print(f"  cn_kmer: {cn_kmer.shape}, cn_var: {cn_var.shape}, "
          f"founders: {len(founders)}, bubbles: {len(bubble_chrom):,}")
    print(f"  block_mode={args.block_mode}, window_bp={args.window_bp:,}")

    # Read manifest
    manifest = pd.read_csv(args.manifest, sep="\t")
    if args.limit:
        manifest = manifest.head(args.limit)
    print(f"  manifest: {len(manifest)} samples")

    # Build job list (skip already-done)
    jobs = []
    for _, row in manifest.iterrows():
        sample_id = str(row["sample_id"])
        out_path = os.path.join(args.out_dir, f"{sample_id}.tsv")
        reads = [row["reads_path"]]
        if "reads_path2" in row and pd.notna(row["reads_path2"]):
            reads.append(row["reads_path2"])
        jobs.append((sample_id, reads, (cn_kmer, kmer_index), cn_var, var_meta,
                     bubble_meta, founders, out_path, args.threads,
                     args.block_mode, args.window_bp))

    n_already = sum(1 for j in jobs if os.path.exists(j[7]))
    print(f"  {n_already} samples already done, {len(jobs) - n_already} to process")

    # Run
    t0 = time.time()
    if args.workers == 1:
        for job in jobs:
            sid, status, rt = process_one(job)
            print(f"[{time.strftime('%H:%M:%S')}] {sid}: {status} ({rt:.0f}s)")
    else:
        with mp.Pool(args.workers) as pool:
            for sid, status, rt in pool.imap_unordered(process_one, jobs):
                print(f"[{time.strftime('%H:%M:%S')}] {sid}: {status} ({rt:.0f}s)")
    print(f"\nTotal time: {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
