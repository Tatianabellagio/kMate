"""
Per-sample driver: BAM or paired FASTQs in → per-VCF-record alt-allele
frequency table out.

Pipeline:
  1. Count k-mers in reads against the precomputed unique-k-mer set
  2. Estimate coverage λ from total counts
  3. EM solve on cn_kmer matrix → founder frequency vector h
  4. Project h to per-record alt-allele frequencies via cn_var:
        f_alt(record) = h · cn_var[:, record]
        where cn_var[f, r] = 1 if founder f's GT at record r is the alt allele
  5. Write TSV: chrom, pos, ref, alt, alt_freq

Input prerequisites (precomputed once):
  - cn_kmer matrix  (data/cn_full_<chrom>.cn.npz)  — F × K_kmer
  - cn_var matrix   (built from the biallelic.norm VCF + sample list) — F × N_records
"""
from __future__ import annotations
import sys, os, time, argparse
sys.path.insert(0, os.path.join(os.path.dirname(__file__)))
import numpy as np
from scipy.sparse import load_npz, hstack
from em_solver import solve_em
from kmer_count import count_kmers_in_bam, count_kmers_in_fasta, estimate_coverage
from block_em import (define_windows, assign_kmers_to_blocks,
                      assign_records_to_blocks, solve_em_per_block,
                      project_blocks_to_records)


def load_cn_kmer_genomewide(prefix_template: str, chroms=("Chr1","Chr2","Chr3","Chr4","Chr5")):
    """Load and concat per-chrom cn_kmer matrices + bubble position metadata.

    Returns:
        cn:           F × K sparse cn matrix, all chroms concatenated
        kmer_index:   K-vec of canonical k-mer strings
        bubble_id:    K-vec of global bubble id (offset across chroms)
        bubble_chrom: B-vec of chromosome per bubble (B = total bubbles)
        bubble_start: B-vec of bubble start position
        bubble_end:   B-vec of bubble end position
        founders:     F-vec of founder names
    """
    cn_blocks, kmer_indices, bubble_ids = [], [], []
    bubble_chrom_list, bubble_start_list, bubble_end_list = [], [], []
    founders = None
    bubble_offset = 0
    for c in chroms:
        cn_path = prefix_template + f"_{c}.cn.npz"
        meta_path = prefix_template + f"_{c}.meta.npz"
        if not os.path.exists(cn_path):
            print(f"  WARN: {c} cn missing")
            continue
        cn_b = load_npz(cn_path)
        meta = np.load(meta_path, allow_pickle=True)
        cn_blocks.append(cn_b)
        kmer_indices.append(meta["kmer_index"])
        bubble_ids.append(meta["bubble_id"] + bubble_offset)
        bubble_chrom_list.append(meta["bubble_chrom"])
        bubble_start_list.append(meta["bubble_start"])
        bubble_end_list.append(meta["bubble_end"])
        bubble_offset += int(meta["bubble_id"].max()) + 1
        if founders is None:
            founders = meta["founders"]
    cn = hstack(cn_blocks, format="csr")
    return (cn,
            np.concatenate(kmer_indices),
            np.concatenate(bubble_ids),
            np.concatenate(bubble_chrom_list),
            np.concatenate(bubble_start_list),
            np.concatenate(bubble_end_list),
            founders)


def run_one_sample(
    cn_kmer,                  # F × K_kmer (sparse or dense)
    kmer_index,               # K_kmer-vector of canonical k-mer strings
    cn_var,                   # F × N_records (sparse) — for projection
    var_meta,                 # dict with chrom, pos, ref, alt arrays (length N_records)
    founders,                 # F-vector of founder names
    reads_input,              # str or list[str] — BAM or paired FASTQs
    sample_name: str,
    output_path: str,
    em_max_iter: int = 200,
    threads: int = 8,
    block_mode: str = "window",          # "window" | "global"
    window_bp: int = 200_000,
    bubble_meta: dict | None = None,     # required for block_mode="window":
                                         #   keys bubble_id, bubble_chrom,
                                         #   bubble_start, bubble_end
    global_anchor_weight: float = 0.0,   # NEW: λ for per-window EM anchor toward h_global
):
    """Run the full pipeline on one sample.

    block_mode="window": HAFpipe-style per-window EM. Each genomic window
        (default 200 kb) gets its own founder freq vector h. Per-record
        alt freqs use the local h of their window. Captures mosaic ancestry
        in evolved samples.

    block_mode="global": single genome-wide h. Use for SEEDMIX-like F0 pools
        or as a comparison baseline.
    """
    F = cn_kmer.shape[0]
    K = cn_kmer.shape[1]

    # 1. Count k-mers
    print(f"[{sample_name}] Counting {K:,} k-mers in reads...", flush=True)
    t = time.time()
    if isinstance(reads_input, str) and reads_input.endswith(".bam"):
        cd = count_kmers_in_bam(reads_input, list(kmer_index), k=31, threads=threads, hash_size="3G")
    else:
        cd = count_kmers_in_fasta(reads_input, list(kmer_index), k=31, threads=threads, hash_size="3G")
    counts = np.array([cd[km] for km in kmer_index], dtype=np.int64)
    print(f"  [{time.time()-t:.0f}s] nonzero {(counts>0).sum():,}/{K:,}", flush=True)

    # 2. Densify cn to float32. Float32 throughout avoids the silent dtype
    # upcast in numpy matmul that was the EM bottleneck (>10× slowdown).
    import gc
    cn_dense = np.asarray(cn_kmer.todense() if hasattr(cn_kmer, "todense") else cn_kmer).astype(np.float32)
    gc.collect()
    ac = cn_dense.sum(axis=0)
    cov = counts.sum() * F / max(1, ac.sum())
    print(f"  coverage estimate (k-mer): {cov:.1f}×, cn_dense: {cn_dense.nbytes/1e9:.1f} GB", flush=True)

    if block_mode == "window":
        # Per-window EM
        if bubble_meta is None:
            raise ValueError("block_mode='window' requires bubble_meta")
        blocks = define_windows(bubble_meta["bubble_chrom"],
                                bubble_meta["bubble_start"],
                                bubble_meta["bubble_end"],
                                window_bp=window_bp)
        n_blocks = len(blocks)
        print(f"  defined {n_blocks} windows of {window_bp:,} bp", flush=True)
        kmer_block = assign_kmers_to_blocks(
            bubble_meta["bubble_id"],
            bubble_meta["bubble_chrom"],
            bubble_meta["bubble_start"],
            bubble_meta["bubble_end"],
            blocks,
        )
        # K-mers can be filtered to nonzero counts globally before solve
        # (the per-block EM also filters internally, but pre-filtering
        # avoids constructing the full block subset). Keep all for now —
        # block_em.solve_em_per_block does the right thing.
        t = time.time()
        h_blocks, status, global_h = solve_em_per_block(
            counts.astype(np.float32), cn_dense, kmer_block, n_blocks,
            cov, em_max_iter=em_max_iter, tol=1e-7,
            min_kmers_per_block=200, verbose=True,
            global_anchor_weight=global_anchor_weight,
        )
        print(f"  block-EM total: {time.time()-t:.0f}s "
              f"({(status==0).sum()}/{n_blocks} local fits, "
              f"{(status==1).sum()} global fallbacks)", flush=True)
        del cn_dense
        gc.collect()

        # Project: smooth interpolation between adjacent blocks on same chrom.
        # Removes step discontinuities at block boundaries (HAFpipe-style).
        rec_chrom = np.asarray(var_meta["chrom"])
        rec_pos = np.asarray(var_meta["pos"])
        rec_block = assign_records_to_blocks(rec_chrom, rec_pos, blocks)
        freqs = project_blocks_to_records(
            h_blocks, status, global_h, cn_var, rec_block,
            record_chrom=rec_chrom, record_pos=rec_pos, blocks=blocks,
            smooth=True,
        )

        # Optional: dump per-block h
        h_path = output_path.replace(".tsv", ".h_blocks.npz")
        if h_path != output_path:
            np.savez(h_path,
                     h_blocks=h_blocks,
                     status=status,
                     global_h=global_h,
                     block_chrom=np.array([b.chrom for b in blocks]),
                     block_start=np.array([b.start for b in blocks]),
                     block_end=np.array([b.end for b in blocks]),
                     founders=founders)
            print(f"  per-block h → {h_path}", flush=True)

    elif block_mode == "global":
        # Filter k-mers to those with c>0
        nz = counts > 0
        n_nz = int(nz.sum())
        print(f"  filtering EM to nonzero-count k-mers: "
              f"{n_nz:,}/{K:,} ({n_nz/K:.1%})", flush=True)
        cn_em = np.ascontiguousarray(cn_dense[:, nz])
        counts_em = counts[nz].astype(np.float32)
        del cn_dense
        gc.collect()

        t = time.time()
        h, info = solve_em(counts_em, cn_em, cov, max_iter=em_max_iter, tol=1e-7)
        print(f"  EM solved in {info['iterations']} iters [{time.time()-t:.0f}s]", flush=True)
        print(f"  effective n founders (1/Σh²): {1/np.sum(h**2):.1f}", flush=True)

        # Project: cn_var.T @ h (uniform h across genome)
        freqs = cn_var.T @ h
        if hasattr(freqs, "toarray"):
            freqs = np.asarray(freqs).flatten()
    else:
        raise ValueError(f"unknown block_mode: {block_mode!r}")

    # Write per-record TSV (same schema for both modes)
    n_var = len(var_meta["chrom"])
    with open(output_path, "w") as f:
        f.write("chrom\tpos\tref_len\talt_len\talt_freq\n")
        for i in range(n_var):
            f.write(f"{var_meta['chrom'][i]}\t{var_meta['pos'][i]}\t"
                    f"{var_meta['ref_len'][i]}\t{var_meta['alt_len'][i]}\t"
                    f"{freqs[i]:.5f}\n")
    print(f"  wrote {n_var:,} records to {output_path}", flush=True)
    return None, freqs


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--cn-kmer-prefix", required=True,
                    help="prefix used by build_kmer_cn (e.g., data/cn_full)")
    ap.add_argument("--cn-var", required=True, help="path to cn_var .npz file")
    ap.add_argument("--cn-var-meta", required=True, help="path to cn_var meta .npz")
    ap.add_argument("--reads", required=True, nargs="+", help="BAM or list of FASTQs")
    ap.add_argument("--sample", required=True, help="sample name")
    ap.add_argument("--out", required=True, help="output TSV path")
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--block-mode", default="window",
                    choices=["window", "global"],
                    help="window: HAFpipe-style per-window EM (default). "
                         "global: single genome-wide h.")
    ap.add_argument("--window-bp", type=int, default=200_000,
                    help="window size for block_mode=window (default 200 kb)")
    ap.add_argument("--global-anchor-weight", type=float, default=0.0,
                    help="λ for per-window EM Dirichlet anchor toward chrom-wide "
                         "h_global. 0 = legacy MLE (default). 0.05–0.5 = mild "
                         "anchor — helps low-evidence windows match the global "
                         "solution while letting high-evidence windows adapt.")
    args = ap.parse_args()

    (cn_kmer, kmer_index, bubble_id, bubble_chrom,
     bubble_start, bubble_end, founders) = load_cn_kmer_genomewide(
        args.cn_kmer_prefix)
    cn_var = load_npz(args.cn_var)
    var_meta_raw = np.load(args.cn_var_meta, allow_pickle=True)
    var_meta = {
        "chrom": var_meta_raw["chrom"],
        "pos": var_meta_raw["pos"],
        "ref_len": var_meta_raw["ref_len"],
        "alt_len": var_meta_raw["alt_len"],
    }
    bubble_meta = {
        "bubble_id": bubble_id,
        "bubble_chrom": bubble_chrom,
        "bubble_start": bubble_start,
        "bubble_end": bubble_end,
    }
    reads = args.reads if len(args.reads) > 1 else args.reads[0]
    run_one_sample(cn_kmer, kmer_index, cn_var, var_meta, founders,
                    reads, args.sample, args.out, threads=args.threads,
                    block_mode=args.block_mode, window_bp=args.window_bp,
                    bubble_meta=bubble_meta,
                    global_anchor_weight=args.global_anchor_weight)
