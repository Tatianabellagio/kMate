"""
Per-chromosome kMate driver — founder-mixture (h) estimation + AF projection.

Two estimators only:
  * global  — one h per chromosome (selfing / inbred / F0 pools, e.g. SEEDMIX).
  * window  — per-window h for recombinant pools (the production "star2" recipe:
              --block-mode window --window-bp 10000 --global-anchor-weight 0.3
              --hmm-smooth-passes 5 --hmm-smooth-alpha 0.5). These are the
              defaults below, so `--block-mode window` alone reproduces it.

Memory note: instead of loading the genome-wide kmer_pa matrix (~74 GB dense
float32 for 80M k-mers × 231 founders), we process one chromosome at a time
(~5× lower peak). Per-chrom h agrees with a genome-wide solve to ~0.1%.

Usage:
    python per_sample_per_chrom.py \\
        --kmer-pa-prefix data/kmer_pa_231/kmer_pa \\
        --var-pa       panel/arch3/chr1/var_pa_231_arch3_chr1.var_pa.npz \\
        --var-called panel/arch3/chr1/var_pa_231_arch3_chr1.var_called.npz \\
        --var-meta  panel/arch3/chr1/var_pa_231_arch3_chr1.meta.npz \\
        --reads R1.fq R2.fq --sample <name> --out <name>.tsv \\
        --threads 8 --chroms Chr1 --block-mode global

Output: per-record TSV (chrom, pos, ref_len, alt_len, alt_freq, info, n_called, se).
"""
from __future__ import annotations
import argparse, gc, os, sys, time
import numpy as np
from scipy.sparse import load_npz
from .em_solver import solve_em
from .kmer_count import count_kmers_in_bam, count_kmers_in_fasta, query_kmer_db
from .block_em import (define_windows, assign_kmers_to_blocks,
                      assign_records_to_blocks, solve_em_per_block,
                      project_blocks_to_records, BlockSpec)
from .block_haplotype_em import smooth_h_across_blocks


def _count_and_load_kmer_pa_dense(chrom, kmer_pa_prefix, reads_input, threads,
                                  kmer_db=None, hash_size="3G"):
    """Load one chrom's kmer_pa + meta, count k-mers in reads, densify to float32.

    Shared between global and window modes. Returns
    (kmer_pa_dense, counts, meta, cov, F, K) or
    (None, None, None, 0.0, 0, 0) if the chrom is missing.

    kmer_db: optional path to a prebuilt Jellyfish DB (from build_kmer_db over
    the full read pool). When given, this chrom's k-mers are *queried* against
    it instead of re-scanning the reads — the count-once / query-per-chrom path
    that avoids the ~5x redundant read scan across chroms. Counts are identical
    to the per-chrom count path (same canonical hash). When None, falls back to
    counting the reads directly (legacy behavior, unchanged).
    """
    cn_path = kmer_pa_prefix + f"_{chrom}.kmer_pa.npz"
    meta_path = kmer_pa_prefix + f"_{chrom}.meta.npz"
    if not os.path.exists(cn_path):
        print(f"  [{chrom}] kmer_pa missing — skip")
        return None, None, None, 0.0, 0, 0

    kmer_pa = load_npz(cn_path)
    meta = np.load(meta_path, allow_pickle=True)
    kmer_index = meta["kmer_index"]
    F, K = kmer_pa.shape
    print(f"  [{chrom}] kmer_pa F={F}, K={K:,}", flush=True)

    t = time.time()
    if kmer_db is not None:
        cd = query_kmer_db(kmer_db, list(kmer_index), k=31)
        verb = "query"
    elif isinstance(reads_input, str) and reads_input.endswith(".bam"):
        cd = count_kmers_in_bam(reads_input, list(kmer_index), k=31, threads=threads, hash_size=hash_size)
        verb = "count"
    else:
        cd = count_kmers_in_fasta(reads_input, list(kmer_index), k=31, threads=threads, hash_size=hash_size)
        verb = "count"
    counts = np.array([cd[km] for km in kmer_index], dtype=np.int64)
    print(f"  [{chrom}] {time.time()-t:.0f}s {verb}: nonzero {(counts>0).sum():,}/{K:,}", flush=True)

    kmer_pa_dense = np.asarray(kmer_pa.todense() if hasattr(kmer_pa, "todense") else kmer_pa).astype(np.float32)
    del kmer_pa
    gc.collect()
    ac = kmer_pa_dense.sum(axis=0)
    cov = counts.sum() * F / max(1, ac.sum())
    print(f"  [{chrom}] cov estimate: {cov:.1f}×, kmer_pa_dense: {kmer_pa_dense.nbytes/1e9:.1f} GB", flush=True)
    return kmer_pa_dense, counts, meta, cov, F, K


# Module-level globals populated by main(): var_called sparse matrix.
# _project_with_called_mask reads this so we don't plumb it through every
# function signature.
_VAR_CALLED = None  # scipy.sparse, founder × variant; 1 if GT != ./.

def _project_with_called_mask(var_pa_chrom, h, idx):
    """Project h through var_pa with per-record renormalization by the called mask.

    Returns (alt_freq, info): AF = (h@var_pa)/(h@var_called) and info = the
    projection denominator, i.e. the h-weighted called mass per record (ones when
    no called mask is loaded). At records where some founders are ./., their
    h-mass is excluded from both numerator and denominator; equals AC/AN under
    uniform h.
    """
    freqs = var_pa_chrom.T @ h
    if hasattr(freqs, "toarray"):
        freqs = np.asarray(freqs).flatten()
    if _VAR_CALLED is None:
        return freqs, np.ones(len(idx), dtype=freqs.dtype)
    called_weight = _VAR_CALLED[:, idx].T @ h
    if hasattr(called_weight, "toarray"):
        called_weight = np.asarray(called_weight).flatten()
    safe = np.maximum(called_weight, np.array(1e-12, dtype=freqs.dtype))
    return (freqs / safe).astype(freqs.dtype), called_weight.astype(freqs.dtype)


def run_one_chrom_global(chrom, kmer_pa_prefix, var_pa, var_meta, reads_input, threads,
                         em_max_iter=200, kmer_weight="uniform", kmer_db=None,
                         hash_size="3G", h_only=False):
    """Global-mode (single h per chrom) EM + projection.

    kmer_weight: "uniform" (ω_k=1, MLE) or "inv_mb" (ω_k=1/m_b per-bubble
    de-replication; m_b = #k-mers sharing the k-mer's bubble_id).
    kmer_db: optional prebuilt Jellyfish DB (count-once path; see
    _count_and_load_kmer_pa_dense).
    """
    t_chrom = time.time()
    kmer_pa_dense, counts, meta, cov, F, K = _count_and_load_kmer_pa_dense(
        chrom, kmer_pa_prefix, reads_input, threads, kmer_db=kmer_db,
        hash_size=hash_size)
    if kmer_pa_dense is None:
        return None, None, None, None, 0.0

    # Filter to nonzero-count k-mers (the EM only needs those)
    nz = counts > 0
    kmer_pa_em = np.ascontiguousarray(kmer_pa_dense[:, nz])
    counts_em = counts[nz].astype(np.float32)
    del kmer_pa_dense
    gc.collect()

    # Optional per-bubble de-replication weight ω_k = 1/m_b (on nz k-mers).
    omega = None
    if kmer_weight == "inv_mb":
        bid = np.asarray(meta["bubble_id"]).astype(np.int64)
        m_b = np.bincount(bid)[bid].astype(np.float32)        # #k-mers per bubble
        omega = (1.0 / m_b[nz]).astype(np.float32)
        print(f"  [{chrom}] ω_k=1/m_b weighting: m_b median={np.median(m_b[nz]):.0f} "
              f"max={m_b[nz].max():.0f}", flush=True)

    t = time.time()
    h, info = solve_em(counts_em, kmer_pa_em, cov, max_iter=em_max_iter, tol=1e-7,
                       omega=omega)
    print(f"  [{chrom}] EM solved in {info['iterations']} iters [{time.time()-t:.0f}s]; "
          f"eff_n_founders = {1/np.sum(h**2):.1f}", flush=True)
    del kmer_pa_em, counts_em
    gc.collect()

    if h_only:
        # MOI / founder-mixture use: skip the AF projection (no var_pa load, no
        # per-record output) — we only need h.
        print(f"  [{chrom}] {time.time()-t_chrom:.0f}s total — h-only (projection skipped)",
              flush=True)
        return None, None, None, h, time.time() - t_chrom

    rec_chrom = np.asarray(var_meta["chrom"]).astype(str)
    idx = np.where(rec_chrom == str(chrom))[0]
    var_pa_chrom = var_pa[:, idx]
    freqs, info = _project_with_called_mask(var_pa_chrom, h, idx)
    info = info.astype(np.float32)
    elapsed = time.time() - t_chrom
    print(f"  [{chrom}] {elapsed:.0f}s total — {len(idx):,} records projected", flush=True)
    return idx, freqs, info, h, elapsed


def load_blocks_tsv(path, chrom):
    """BlockSpec list for one chrom from an LD-block TSV
    (cols: chrom start_pos end_pos n_variants; TAIR10 1-based)."""
    import pandas as pd
    bt = pd.read_csv(path, sep="\t")
    bt = bt[bt["chrom"].astype(str) == str(chrom)]
    return [BlockSpec(chrom=str(chrom), start=int(s), end=int(e))
            for s, e in zip(bt["start_pos"], bt["end_pos"])]


def run_one_chrom_window(chrom, kmer_pa_prefix, var_pa, var_meta, reads_input, threads,
                         window_bp=10_000, em_max_iter=200,
                         global_anchor_weight=0.3,
                         hmm_smooth_passes=5,
                         hmm_smooth_alpha=0.5,
                         hmm_smooth_recomb_rate=4e-8,
                         kmer_weight="uniform", kmer_db=None, hash_size="3G",
                         blocks_tsv=None, min_kmers_per_block=200):
    """Window-mode per-chrom EM + smooth projection (production "star2" recipe).

    Fixed-bp windows; per-window EM optionally anchored toward the chrom-wide
    h_global (global_anchor_weight) and post-smoothed across windows
    (Li-Stephens-style, hmm_smooth_*). Each window's h projects its records
    through var_pa with the missing-aware (called-mask) normalization.
    """
    t_chrom = time.time()
    kmer_pa_dense, counts, meta, cov, F, K = _count_and_load_kmer_pa_dense(
        chrom, kmer_pa_prefix, reads_input, threads, kmer_db=kmer_db,
        hash_size=hash_size)
    if kmer_pa_dense is None:
        return None, None, None, None, 0.0

    bubble_id = meta["bubble_id"]
    bubble_chrom = meta["bubble_chrom"]
    bubble_start = meta["bubble_start"]
    bubble_end = meta["bubble_end"]

    if blocks_tsv:
        blocks = load_blocks_tsv(blocks_tsv, chrom)
        print(f"  [{chrom}] {len(blocks)} LD blocks from {blocks_tsv}", flush=True)
    else:
        blocks = define_windows(bubble_chrom, bubble_start, bubble_end,
                                window_bp=window_bp)
        print(f"  [{chrom}] {len(blocks)} fixed windows of {window_bp:,} bp", flush=True)
    n_blocks = len(blocks)

    kmer_block = assign_kmers_to_blocks(bubble_id, bubble_chrom,
                                        bubble_start, bubble_end, blocks)

    # Optional per-bubble de-replication weight ω_k = 1/m_b for window-mode EM.
    omega = None
    if kmer_weight == "inv_mb":
        m_b = np.bincount(bubble_id)[bubble_id].astype(np.float32)
        omega = (1.0 / m_b).astype(np.float32)
        print(f"  [{chrom}] window ω_k=1/m_b weighting: m_b median={np.median(m_b):.0f} "
              f"max={m_b.max():.0f}", flush=True)

    t = time.time()
    h_blocks, status, global_h = solve_em_per_block(
        counts.astype(np.float32), kmer_pa_dense, kmer_block, n_blocks,
        cov, em_max_iter=em_max_iter, tol=1e-7,
        min_kmers_per_block=min_kmers_per_block, verbose=False,
        global_anchor_weight=global_anchor_weight,
        omega=omega,
    )
    print(f"  [{chrom}] block-EM {time.time()-t:.0f}s "
          f"({(status==0).sum()}/{n_blocks} local fits, "
          f"{(status==1).sum()} fallbacks)", flush=True)
    del kmer_pa_dense
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
    var_pa_chrom = var_pa[:, idx]
    # MAR projection: pass var_called sliced to this chrom so window-mode
    # matches global-mode's (h@var_pa)/(h@var_called) semantics.
    var_called_chrom = _VAR_CALLED[:, idx] if _VAR_CALLED is not None else None

    # Hard window assignment: each record gets the h of its window (already
    # HMM-smoothed across windows above). The projection returns the AF and the
    # per-record info (h-weighted called mass) in one pass.
    rec_block = assign_records_to_blocks(rec_chrom[idx], rec_pos[idx], blocks)
    freqs, info = project_blocks_to_records(
        h_blocks, global_h, var_pa_chrom, rec_block,
        var_called=var_called_chrom,
    )
    info = info.astype(np.float32)
    elapsed = time.time() - t_chrom
    print(f"  [{chrom}] {elapsed:.0f}s total — {len(idx):,} records projected", flush=True)
    return idx, freqs, info, (h_blocks, status, global_h, blocks), elapsed


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--kmer-pa-prefix", required=True)
    ap.add_argument("--var-pa", required=False,
                    help="Required unless --h-only (the AF projection target).")
    ap.add_argument("--var-meta", required=False,
                    help="Required unless --h-only.")
    ap.add_argument("--h-only", action="store_true",
                    help="Estimate only the founder mixture h; skip the per-record "
                         "AF projection and TSV. Writes just <out>.h_per_chrom.npz "
                         "(founder order = kmer_pa). For MOI / outcrossing or any use "
                         "that needs founder frequencies, not per-variant AF — faster, "
                         "no var_pa load. Global mode only.")
    ap.add_argument("--var-called", default=None,
                    help="Path to var_called.npz (founder × variant 1/0 mask "
                         "of called genotypes). When provided, AF projection "
                         "uses (h@var_pa)/(h@var_called) to correctly handle "
                         "./. cells. Auto-detected next to --var-pa if not given.")
    ap.add_argument("--reads", required=True, nargs="+")
    ap.add_argument("--kmer-db", default=None,
                    help="Path to a prebuilt Jellyfish DB (build_kmer_db) over the "
                         "full read pool. When given, each chrom QUERIES it instead of "
                         "re-counting the reads — the count-once / query-per-chrom "
                         "speedup. Counts are identical to the per-chrom count path. "
                         "When omitted, the reads in --reads are counted per chrom "
                         "(legacy behavior).")
    ap.add_argument("--sample", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--hash-size", default="3G",
                    help="Initial Jellyfish hash size for k-mer counting "
                         "(e.g. '3G', '100M'). Auto-grows, so this only sets the "
                         "starting allocation; lower it (e.g. '100M') for small "
                         "pools or memory-capped jobs. Production: 3G.")
    ap.add_argument("--chroms", nargs="+", default=["Chr1","Chr2","Chr3","Chr4","Chr5"])
    ap.add_argument("--kmer-weight", default="uniform", choices=["uniform", "inv_mb"],
                    help="Per-k-mer EM weight ω_k. 'uniform' = MLE; 'inv_mb' = 1/m_b "
                         "per-bubble de-replication (removes imbalanced-design over-credit). "
                         "Applied in both global and window modes.")
    ap.add_argument("--block-mode", default="global", choices=["global", "window"],
                    help="global: one h per chrom (selfing / inbred / F0 pools). "
                         "window: per-window h for recombinant pools (production "
                         "'star2' recipe; window defaults below reproduce it).")
    ap.add_argument("--window-bp", type=int, default=10_000,
                    help="fixed-window size for --block-mode window (production: 10 kb)")
    ap.add_argument("--blocks-tsv", default=None,
                    help="LD-block TSV (chrom start_pos end_pos n_variants; TAIR10 coords). "
                         "In window mode, use these blocks as the windows instead of "
                         "fixed --window-bp windows. Records/k-mers outside any block "
                         "get no local fit (NaN).")
    ap.add_argument("--min-kmers-per-block", type=int, default=200,
                    help="min k-mers for a local per-block EM fit; below this the block "
                         "falls back to the chrom-wide h. Lower it to probe thin LD blocks.")
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
    args = ap.parse_args()

    print(f"=== {args.sample} (per-chrom {args.block_mode} mode"
          f"{', h-only' if args.h_only else ''}) ===", flush=True)
    print(f"  reads: {args.reads}", flush=True)
    t0 = time.time()

    global _VAR_CALLED
    if args.h_only:
        if args.block_mode != "global":
            sys.exit("ERROR: --h-only is supported in global mode only")
        var_pa = None
        var_meta = None
        n_records = None
        _VAR_CALLED = None
        print("  h-only mode: skipping var_pa load + AF projection "
              "(writing only the h vector)", flush=True)
    else:
        if not args.var_pa or not args.var_meta:
            sys.exit("ERROR: --var-pa and --var-meta are required unless --h-only")
        var_pa = load_npz(args.var_pa)
        var_meta = np.load(args.var_meta, allow_pickle=True)
        n_records = var_pa.shape[1]
        print(f"  var_pa: {var_pa.shape}, n_records: {n_records:,}", flush=True)

    # Optional: load var_called mask for missing-aware AF projection.
    called_path = args.var_called if not args.h_only else None
    if called_path is None and not args.h_only:
        guess = args.var_pa.replace(".var_pa.npz", ".var_called.npz")
        if guess != args.var_pa and os.path.exists(guess):
            called_path = guess
    if called_path and os.path.exists(called_path):
        _VAR_CALLED = load_npz(called_path)
        print(f"  var_called: {_VAR_CALLED.shape}, nnz={_VAR_CALLED.nnz:,} "
              f"(loaded from {called_path}; AF projection will divide by h@called)",
              flush=True)
    else:
        _VAR_CALLED = None
        if not args.h_only:
            print(f"  var_called: NOT FOUND — AF projection treats ./. as REF "
                  f"(legacy behavior). Rebuild with build_var_pa.py for "
                  f"missing-aware projection.", flush=True)

    reads_input = args.reads if len(args.reads) > 1 else args.reads[0]

    kmer_db = args.kmer_db
    if kmer_db is not None:
        if not os.path.exists(kmer_db):
            sys.exit(f"ERROR: --kmer-db {kmer_db} does not exist")
        print(f"  kmer_db: {kmer_db} (count-once: querying prebuilt DB per chrom)",
              flush=True)

    h_save = {}  # global: h_per_chrom[chrom] = h. window: per-chrom block packs.
    if not args.h_only:
        freqs_global = np.full(n_records, np.nan, dtype=np.float32)
        info_global  = np.full(n_records, np.nan, dtype=np.float32)
        # Panel-level per-record called counts (h-independent QC metric). When
        # var_called is unavailable, assume the panel size F (all called).
        F_total = var_pa.shape[0]
        if _VAR_CALLED is not None:
            n_called_per_rec = np.asarray(_VAR_CALLED.sum(axis=0)).flatten().astype(np.int32)
        else:
            n_called_per_rec = np.full(n_records, F_total, dtype=np.int32)

    for chrom in args.chroms:
        if args.block_mode == "global":
            idx, freqs, info, h, _ = run_one_chrom_global(
                chrom, args.kmer_pa_prefix, var_pa, var_meta,
                reads_input, args.threads, kmer_weight=args.kmer_weight,
                kmer_db=kmer_db, hash_size=args.hash_size, h_only=args.h_only)
            if h is None:
                continue
            h_save[chrom] = h
            if args.h_only:
                gc.collect()
                continue
        else:  # window
            idx, freqs, info, pack, _ = run_one_chrom_window(
                chrom, args.kmer_pa_prefix, var_pa, var_meta,
                reads_input, args.threads,
                window_bp=args.window_bp,
                global_anchor_weight=args.global_anchor_weight,
                hmm_smooth_passes=args.hmm_smooth_passes,
                hmm_smooth_alpha=args.hmm_smooth_alpha,
                hmm_smooth_recomb_rate=args.hmm_smooth_recomb_rate,
                kmer_weight=args.kmer_weight,
                kmer_db=kmer_db, hash_size=args.hash_size,
                blocks_tsv=args.blocks_tsv,
                min_kmers_per_block=args.min_kmers_per_block,
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

    if not args.h_only:
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
    h_path = args.out.replace(".tsv", suffix) if args.out.endswith(".tsv") else args.out + suffix
    if h_path != args.out:
        if args.h_only:
            # founder order = kmer_pa (the EM's founder axis); var_meta not loaded
            km = np.load(f"{args.kmer_pa_prefix}_{args.chroms[0]}.meta.npz", allow_pickle=True)
            founders = np.asarray(km["founders"]) if "founders" in km.files else np.array([], dtype=object)
        else:
            founders = np.asarray(var_meta.get("founders", np.array([], dtype=object)))
        np.savez(h_path, founders=founders, **h_save)
        print(f"  wrote h -> {h_path}", flush=True)

    print(f"  TOTAL: {time.time()-t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
