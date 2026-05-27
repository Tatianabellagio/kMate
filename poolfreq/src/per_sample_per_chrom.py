"""
Per-chromosome kMate driver — founder-mixture (h) estimation + AF projection.

Two estimators only:
  * global  — one h per chromosome (selfing / inbred / F0 pools, e.g. SEEDMIX).
  * window  — per-window h for recombinant pools (the production "star2" recipe:
              --block-mode window --window-bp 10000 --global-anchor-weight 0.3
              --hmm-smooth-passes 5 --hmm-smooth-alpha 0.5). These are the
              defaults below, so `--block-mode window` alone reproduces it.

Memory note: instead of loading the genome-wide cn_kmer matrix (~74 GB dense
float32 for 80M k-mers × 231 founders), we process one chromosome at a time
(~5× lower peak). Per-chrom h agrees with a genome-wide solve to ~0.1%.

Usage:
    python per_sample_per_chrom.py \\
        --cn-kmer-prefix data/cn_full_231/cn \\
        --cn-var       arch3/chr1/cn_var_231_arch3_chr1.cn_var.npz \\
        --cn-var-called arch3/chr1/cn_var_231_arch3_chr1.cn_var_called.npz \\
        --cn-var-meta  arch3/chr1/cn_var_231_arch3_chr1.meta.npz \\
        --reads R1.fq R2.fq --sample <name> --out <name>.tsv \\
        --threads 8 --chroms Chr1 --block-mode global

Output: per-record TSV (chrom, pos, ref_len, alt_len, alt_freq, info, n_called, se).
"""
from __future__ import annotations
import argparse, gc, os, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from scipy.sparse import load_npz
from em_solver import solve_em
from kmer_count import count_kmers_in_bam, count_kmers_in_fasta
from block_em import (define_windows, assign_kmers_to_blocks,
                      assign_records_to_blocks, solve_em_per_block,
                      project_blocks_to_records, BlockSpec)
from block_haplotype_em import smooth_h_across_blocks


def _count_and_load_cn_dense(chrom, cn_prefix, reads_input, threads):
    """Load one chrom's cn_kmer + meta, count k-mers in reads, densify to float32.

    Shared between global and window modes. Returns
    (cn_dense, counts, meta, cov, F, K) or
    (None, None, None, 0.0, 0, 0) if the chrom is missing.
    """
    cn_path = cn_prefix + f"_{chrom}.cn.npz"
    meta_path = cn_prefix + f"_{chrom}.meta.npz"
    if not os.path.exists(cn_path):
        print(f"  [{chrom}] cn missing — skip")
        return None, None, None, 0.0, 0, 0

    cn_kmer = load_npz(cn_path)
    meta = np.load(meta_path, allow_pickle=True)
    kmer_index = meta["kmer_index"]
    F, K = cn_kmer.shape
    print(f"  [{chrom}] cn_kmer F={F}, K={K:,}", flush=True)

    t = time.time()
    if isinstance(reads_input, str) and reads_input.endswith(".bam"):
        cd = count_kmers_in_bam(reads_input, list(kmer_index), k=31, threads=threads, hash_size="3G")
    else:
        cd = count_kmers_in_fasta(reads_input, list(kmer_index), k=31, threads=threads, hash_size="3G")
    counts = np.array([cd[km] for km in kmer_index], dtype=np.int64)
    print(f"  [{chrom}] {time.time()-t:.0f}s count: nonzero {(counts>0).sum():,}/{K:,}", flush=True)

    cn_dense = np.asarray(cn_kmer.todense() if hasattr(cn_kmer, "todense") else cn_kmer).astype(np.float32)
    del cn_kmer
    gc.collect()
    ac = cn_dense.sum(axis=0)
    cov = counts.sum() * F / max(1, ac.sum())
    print(f"  [{chrom}] cov estimate: {cov:.1f}×, cn_dense: {cn_dense.nbytes/1e9:.1f} GB", flush=True)
    return cn_dense, counts, meta, cov, F, K


# Module-level globals populated by main(): cn_var_called sparse matrix.
# _project_with_called_mask reads this so we don't plumb it through every
# function signature.
_CN_VAR_CALLED = None  # scipy.sparse, founder × variant; 1 if GT != ./.

def _project_with_called_mask(cn_var_chrom, h, idx):
    """Project h through cn_var with per-record renormalization by the called mask.

    AF_est[r] = (h @ cn_var)[r] / (h @ cn_var_called)[r]

    Handles missing GTs: at records where a subset of founders is ./., their
    h-mass is excluded from BOTH numerator and denominator. Equals AC/AN when
    h is uniform. Falls back to plain (h @ cn_var) if cn_var_called is absent.
    """
    freqs = cn_var_chrom.T @ h
    if hasattr(freqs, "toarray"):
        freqs = np.asarray(freqs).flatten()
    if _CN_VAR_CALLED is None:
        return freqs
    called_chrom = _CN_VAR_CALLED[:, idx]
    called_weight = called_chrom.T @ h
    if hasattr(called_weight, "toarray"):
        called_weight = np.asarray(called_weight).flatten()
    safe = np.maximum(called_weight, np.array(1e-12, dtype=freqs.dtype))
    return (freqs / safe).astype(freqs.dtype)


def run_one_chrom_global(chrom, cn_prefix, cn_var, var_meta, reads_input, threads,
                         em_max_iter=200):
    """Global-mode (single h per chrom) EM + projection."""
    t_chrom = time.time()
    cn_dense, counts, meta, cov, F, K = _count_and_load_cn_dense(
        chrom, cn_prefix, reads_input, threads)
    if cn_dense is None:
        return None, None, None, 0.0

    # Filter to nonzero-count k-mers (the EM only needs those)
    nz = counts > 0
    cn_em = np.ascontiguousarray(cn_dense[:, nz])
    counts_em = counts[nz].astype(np.float32)
    del cn_dense
    gc.collect()

    t = time.time()
    h, info = solve_em(counts_em, cn_em, cov, max_iter=em_max_iter, tol=1e-7)
    print(f"  [{chrom}] EM solved in {info['iterations']} iters [{time.time()-t:.0f}s]; "
          f"eff_n_founders = {1/np.sum(h**2):.1f}", flush=True)
    del cn_em, counts_em
    gc.collect()

    rec_chrom = np.asarray(var_meta["chrom"]).astype(str)
    idx = np.where(rec_chrom == str(chrom))[0]
    cn_var_chrom = cn_var[:, idx]
    freqs = _project_with_called_mask(cn_var_chrom, h, idx)
    # Per-record info = h-mass on called founders at record r (the projection
    # denominator). info ∈ [0, 1]; small info → low-confidence AF.
    if _CN_VAR_CALLED is not None:
        info = np.asarray(_CN_VAR_CALLED[:, idx].T @ h).flatten().astype(np.float32)
    else:
        info = np.ones(len(idx), dtype=np.float32)
    elapsed = time.time() - t_chrom
    print(f"  [{chrom}] {elapsed:.0f}s total — {len(idx):,} records projected", flush=True)
    return idx, freqs, info, h, elapsed


def run_one_chrom_window(chrom, cn_prefix, cn_var, var_meta, reads_input, threads,
                         window_bp=10_000, em_max_iter=200,
                         global_anchor_weight=0.3,
                         hmm_smooth_passes=5,
                         hmm_smooth_alpha=0.5,
                         hmm_smooth_recomb_rate=4e-8,
                         projection_smooth="auto"):
    """Window-mode per-chrom EM + smooth projection (production "star2" recipe).

    Fixed-bp windows; per-window EM optionally anchored toward the chrom-wide
    h_global (global_anchor_weight) and post-smoothed across windows
    (Li-Stephens-style, hmm_smooth_*). Each window's h projects its records
    through cn_var with the missing-aware (called-mask) normalization.
    """
    t_chrom = time.time()
    cn_dense, counts, meta, cov, F, K = _count_and_load_cn_dense(
        chrom, cn_prefix, reads_input, threads)
    if cn_dense is None:
        return None, None, None, None, 0.0

    bubble_id = meta["bubble_id"]
    bubble_chrom = meta["bubble_chrom"]
    bubble_start = meta["bubble_start"]
    bubble_end = meta["bubble_end"]

    blocks = define_windows(bubble_chrom, bubble_start, bubble_end,
                            window_bp=window_bp)
    print(f"  [{chrom}] {len(blocks)} fixed windows of {window_bp:,} bp", flush=True)
    n_blocks = len(blocks)

    kmer_block = assign_kmers_to_blocks(bubble_id, bubble_chrom,
                                        bubble_start, bubble_end, blocks)

    t = time.time()
    h_blocks, status, global_h = solve_em_per_block(
        counts.astype(np.float32), cn_dense, kmer_block, n_blocks,
        cov, em_max_iter=em_max_iter, tol=1e-7,
        min_kmers_per_block=200, verbose=False,
        global_anchor_weight=global_anchor_weight,
    )
    print(f"  [{chrom}] block-EM {time.time()-t:.0f}s "
          f"({(status==0).sum()}/{n_blocks} local fits, "
          f"{(status==1).sum()} fallbacks)", flush=True)
    del cn_dense
    gc.collect()

    # Optional post-EM HMM smoothing across blocks (Li-Stephens style).
    if hmm_smooth_passes > 0:
        block_pos_start = np.array([b.start for b in blocks], dtype=np.int64)
        block_pos_end = np.array([b.end for b in blocks], dtype=np.int64)
        n_eco = h_blocks.shape[1]
        # Only smooth blocks with a successful local fit (status == 0)
        h_dict = {b: h_blocks[b].copy() for b in range(n_blocks) if status[b] == 0}
        n_in = len(h_dict)
        t = time.time()
        h_dict = smooth_h_across_blocks(
            h_dict, block_pos_start, block_pos_end, n_eco,
            recomb_rate=hmm_smooth_recomb_rate,
            alpha=hmm_smooth_alpha,
            n_passes=hmm_smooth_passes)
        for b, h in h_dict.items():
            h_blocks[b] = h
        print(f"  [{chrom}] HMM smoothing: {hmm_smooth_passes} passes α={hmm_smooth_alpha} "
              f"on {n_in}/{n_blocks} blocks in {time.time()-t:.0f}s", flush=True)

    rec_chrom = np.asarray(var_meta["chrom"]).astype(str)
    rec_pos = np.asarray(var_meta["pos"])
    idx = np.where(rec_chrom == str(chrom))[0]
    cn_var_chrom = cn_var[:, idx]
    # MAR projection: pass cn_var_called sliced to this chrom so window-mode
    # matches global-mode's (h@cn_var)/(h@cn_var_called) semantics.
    cn_var_called_chrom = _CN_VAR_CALLED[:, idx] if _CN_VAR_CALLED is not None else None

    # If HMM smoothing already ran, skip projection-time linear interpolation
    # (double-smoothing over-attenuates local signal).
    if projection_smooth == "auto":
        do_smooth = (hmm_smooth_passes == 0)
    elif projection_smooth in ("on", "true", "1"):
        do_smooth = True
    elif projection_smooth in ("off", "false", "0"):
        do_smooth = False
    else:
        raise ValueError(f"projection_smooth must be auto/on/off, got {projection_smooth!r}")

    rec_block = assign_records_to_blocks(rec_chrom[idx], rec_pos[idx], blocks)
    freqs = project_blocks_to_records(
        h_blocks, status, global_h, cn_var_chrom, rec_block,
        record_chrom=rec_chrom[idx], record_pos=rec_pos[idx],
        blocks=blocks, smooth=do_smooth,
        cn_var_called=cn_var_called_chrom,
    )
    # Info = same projection but with cn_var_called as the carrier matrix and
    # no called-mask normalization (raw h-weighted called mass per record).
    if cn_var_called_chrom is not None:
        info = project_blocks_to_records(
            h_blocks, status, global_h, cn_var_called_chrom, rec_block,
            record_chrom=rec_chrom[idx], record_pos=rec_pos[idx],
            blocks=blocks, smooth=do_smooth,
            cn_var_called=None,
        ).astype(np.float32)
    else:
        info = np.ones(len(idx), dtype=np.float32)
    elapsed = time.time() - t_chrom
    print(f"  [{chrom}] {elapsed:.0f}s total — {len(idx):,} records projected", flush=True)
    return idx, freqs, info, (h_blocks, status, global_h, blocks), elapsed


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cn-kmer-prefix", required=True)
    ap.add_argument("--cn-var", required=True)
    ap.add_argument("--cn-var-meta", required=True)
    ap.add_argument("--cn-var-called", default=None,
                    help="Path to cn_var_called.npz (founder × variant 1/0 mask "
                         "of called genotypes). When provided, AF projection "
                         "uses (h@cn_var)/(h@cn_var_called) to correctly handle "
                         "./. cells. Auto-detected next to --cn-var if not given.")
    ap.add_argument("--reads", required=True, nargs="+")
    ap.add_argument("--sample", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--chroms", nargs="+", default=["Chr1","Chr2","Chr3","Chr4","Chr5"])
    ap.add_argument("--block-mode", default="global", choices=["global", "window"],
                    help="global: one h per chrom (selfing / inbred / F0 pools). "
                         "window: per-window h for recombinant pools (production "
                         "'star2' recipe; window defaults below reproduce it).")
    ap.add_argument("--window-bp", type=int, default=10_000,
                    help="fixed-window size for --block-mode window (production: 10 kb)")
    ap.add_argument("--global-anchor-weight", type=float, default=0.3,
                    help="λ for per-window EM Dirichlet anchor toward chrom-wide "
                         "h_global. 0 = pure per-window MLE; production: 0.3.")
    ap.add_argument("--hmm-smooth-passes", type=int, default=5,
                    help="Post-EM Li-Stephens-style smoothing passes on h_blocks. "
                         "0 = off; production: 5.")
    ap.add_argument("--hmm-smooth-alpha", type=float, default=0.5,
                    help="α for HMM smoothing: h_b' = α·h_b + (1-α)·neighbor_avg. "
                         "Smaller = more smoothing. Production: 0.5.")
    ap.add_argument("--hmm-smooth-recomb-rate", type=float, default=4e-8,
                    help="Recomb rate per bp for distance-weighted neighbor averaging "
                         "(default 4e-8 = 4 cM/Mb, A. thaliana).")
    ap.add_argument("--projection-smooth", default="auto",
                    choices=["auto", "on", "off"],
                    help="Linear-interpolation smoothing between adjacent windows "
                         "at projection time. 'auto' (default): off when HMM "
                         "smoothing is on (avoids double-smoothing), on otherwise.")
    args = ap.parse_args()

    print(f"=== {args.sample} (per-chrom {args.block_mode} mode) ===", flush=True)
    print(f"  reads: {args.reads}", flush=True)
    t0 = time.time()

    cn_var = load_npz(args.cn_var)
    var_meta = np.load(args.cn_var_meta, allow_pickle=True)
    n_records = cn_var.shape[1]
    print(f"  cn_var: {cn_var.shape}, n_records: {n_records:,}", flush=True)

    # Optional: load cn_var_called mask for missing-aware AF projection.
    global _CN_VAR_CALLED
    called_path = args.cn_var_called
    if called_path is None:
        guess = args.cn_var.replace(".cn_var.npz", ".cn_var_called.npz")
        if guess != args.cn_var and os.path.exists(guess):
            called_path = guess
    if called_path and os.path.exists(called_path):
        _CN_VAR_CALLED = load_npz(called_path)
        print(f"  cn_var_called: {_CN_VAR_CALLED.shape}, nnz={_CN_VAR_CALLED.nnz:,} "
              f"(loaded from {called_path}; AF projection will divide by h@called)",
              flush=True)
    else:
        _CN_VAR_CALLED = None
        print(f"  cn_var_called: NOT FOUND — AF projection treats ./. as REF "
              f"(legacy behavior). Rebuild with build_cn_var.py for "
              f"missing-aware projection.", flush=True)

    reads_input = args.reads if len(args.reads) > 1 else args.reads[0]

    freqs_global = np.full(n_records, np.nan, dtype=np.float32)
    info_global  = np.full(n_records, np.nan, dtype=np.float32)
    h_save = {}  # global: h_per_chrom[chrom] = h. window: per-chrom block packs.

    # Panel-level per-record called counts (h-independent QC metric). When
    # cn_var_called is unavailable, assume the panel size F (all called).
    F_total = cn_var.shape[0]
    if _CN_VAR_CALLED is not None:
        n_called_per_rec = np.asarray(_CN_VAR_CALLED.sum(axis=0)).flatten().astype(np.int32)
    else:
        n_called_per_rec = np.full(n_records, F_total, dtype=np.int32)

    for chrom in args.chroms:
        if args.block_mode == "global":
            idx, freqs, info, h, _ = run_one_chrom_global(
                chrom, args.cn_kmer_prefix, cn_var, var_meta,
                reads_input, args.threads)
            if idx is None:
                continue
            h_save[chrom] = h
        else:  # window
            idx, freqs, info, pack, _ = run_one_chrom_window(
                chrom, args.cn_kmer_prefix, cn_var, var_meta,
                reads_input, args.threads,
                window_bp=args.window_bp,
                global_anchor_weight=args.global_anchor_weight,
                hmm_smooth_passes=args.hmm_smooth_passes,
                hmm_smooth_alpha=args.hmm_smooth_alpha,
                hmm_smooth_recomb_rate=args.hmm_smooth_recomb_rate,
                projection_smooth=args.projection_smooth,
            )
            if idx is None:
                continue
            h_blocks, status, global_h, blocks = pack
            h_save[f"{chrom}_h_blocks"] = h_blocks
            h_save[f"{chrom}_status"] = status
            h_save[f"{chrom}_global_h"] = global_h
            h_save[f"{chrom}_block_chrom"] = np.array([b.chrom for b in blocks])
            h_save[f"{chrom}_block_start"] = np.array([b.start for b in blocks])
            h_save[f"{chrom}_block_end"] = np.array([b.end for b in blocks])
        freqs_global[idx] = freqs
        info_global[idx]  = info
        gc.collect()

    # SE per record (Wald, using n_called as effective N). NaN if p is NaN or
    # n_called == 0. n_called is the integer count of called founders at each
    # record (h-independent panel QC); info is the h-weighted version.
    with np.errstate(invalid='ignore', divide='ignore'):
        p_clip = np.clip(freqs_global, 0.0, 1.0)
        se_global = np.sqrt(p_clip * (1.0 - p_clip) / np.maximum(n_called_per_rec, 1))
        se_global = np.where(np.isfinite(freqs_global) & (n_called_per_rec > 0),
                             se_global, np.nan).astype(np.float32)

    chrom_arr = np.asarray(var_meta["chrom"])
    pos_arr = np.asarray(var_meta["pos"])
    ref_arr = np.asarray(var_meta["ref_len"])
    alt_arr = np.asarray(var_meta["alt_len"])
    with open(args.out, "w") as f:
        f.write("chrom\tpos\tref_len\talt_len\talt_freq\tinfo\tn_called\tse\n")
        for i in range(n_records):
            af = freqs_global[i]
            af_str = f"{af:.5f}" if np.isfinite(af) else "NaN"
            inf_str = f"{info_global[i]:.5f}" if np.isfinite(info_global[i]) else "NaN"
            nc_str = f"{int(n_called_per_rec[i])}"
            se_str = f"{se_global[i]:.5f}" if np.isfinite(se_global[i]) else "NaN"
            f.write(f"{chrom_arr[i]}\t{pos_arr[i]}\t{ref_arr[i]}\t{alt_arr[i]}\t"
                    f"{af_str}\t{inf_str}\t{nc_str}\t{se_str}\n")
    print(f"  wrote {n_records:,} records to {args.out}", flush=True)

    suffix = ".h_per_chrom.npz" if args.block_mode == "global" else ".h_blocks_per_chrom.npz"
    h_path = args.out.replace(".tsv", suffix)
    if h_path != args.out:
        founders = np.asarray(var_meta.get("founders", np.array([], dtype=object)))
        np.savez(h_path, founders=founders, **h_save)

    print(f"  TOTAL: {time.time()-t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
