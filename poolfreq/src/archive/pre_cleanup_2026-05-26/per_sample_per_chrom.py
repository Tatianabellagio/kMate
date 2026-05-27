"""
Per-chromosome cactus_em driver.

Memory optimization for production scale-out: instead of loading the genome-wide
cn_kmer matrix (80M k-mers × 231 founders → ~74 GB dense float32), process one
chromosome at a time. Each chromosome is ~1/5 of the matrix → peak memory drops
~5×, fits in 32-64 GB SLURM allocations and unlocks the 128 GB memex nodes for
parallel scale-out.

Tradeoff: each chrom's EM uses only that chrom's k-mers as evidence (vs the
joint genome-wide EM in per_sample_driver.py). With ~16 M k-mers per chrom and
a 231-founder simplex, the EM is still massively over-determined — empirically
the per-chrom h vectors agree to within ~0.1% of the genome-wide h.

Usage (drop-in replacement for per_sample_driver.py):

    python per_sample_per_chrom.py \\
        --cn-kmer-prefix data/cn_full_231_v2/cn \\
        --cn-var data/cn_var_231_v2.cn_var.npz \\
        --cn-var-meta data/cn_var_231_v2.meta.npz \\
        --reads R1.fq R2.fq \\
        --sample <name> \\
        --out <name>.tsv \\
        --threads 4

Output: per-VCF-record alt-allele frequency TSV (same schema as
per_sample_driver.py: chrom, pos, ref_len, alt_len, alt_freq).
"""
from __future__ import annotations
import argparse, gc, os, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from scipy.sparse import load_npz
from em_solver import solve_em
from kmer_count import count_kmers_in_bam, count_kmers_in_fasta
from block_em import (define_windows, assign_kmers_to_blocks,
                      assign_kmers_to_blocks_multi,
                      assign_records_to_blocks, solve_em_per_block,
                      project_blocks_to_records,
                      project_blocks_to_records_overlap,
                      BlockSpec)
from block_haplotype_em import smooth_h_across_blocks
from ld_blocks import compute_ld_blocks_gabriel, compute_ld_blocks


def load_bigld_panel_blocks(npz_path, chrom_filter=None):
    """Load fine BigLD blocks from a hapfire_block_index.npz (built by
    build_hapfire_block_index.py from xwu's panel partition). Returns a list
    of BlockSpec.

    The npz uses chromosome IDs '1', '2', ... (panel VCF convention); we
    convert to 'ChrN' to match cactus_em's internal usage. If chrom_filter
    is given, only keep blocks on that chrom.
    """
    npz = np.load(npz_path, allow_pickle=True)
    block_chrom = np.asarray(npz['block_chrom']).astype(str)
    block_pos_start = np.asarray(npz['block_pos_start']).astype(np.int64)
    block_pos_end = np.asarray(npz['block_pos_end']).astype(np.int64)
    blocks = []
    for c, s, e in zip(block_chrom, block_pos_start, block_pos_end):
        target = f'Chr{c}' if not str(c).startswith('Chr') else str(c)
        if chrom_filter is not None and target != chrom_filter:
            continue
        blocks.append(BlockSpec(chrom=target, start=int(s), end=int(e)))
    return blocks


def _count_and_load_cn_dense(chrom, cn_prefix, reads_input, threads,
                              row_normalize: bool = False):
    """Load one chrom's cn_kmer + meta, count k-mers in reads, densify to float32.

    Shared between global and window modes. Returns
    (cn_dense, counts, meta, cov, F, K, K_f) or
    (None, None, None, 0.0, 0, 0, None) if the chrom is missing.

    When row_normalize=True, divides each founder row of cn_dense by K_f
    (= per-founder unique-k-mer count) before returning. K_f is always returned
    in its pre-normalization form so the projection step can divide h by K_f
    to recover mass-domain founder weights.
    """
    cn_path = cn_prefix + f"_{chrom}.cn.npz"
    meta_path = cn_prefix + f"_{chrom}.meta.npz"
    if not os.path.exists(cn_path):
        print(f"  [{chrom}] cn missing — skip")
        return None, None, None, 0.0, 0, 0, None

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

    # K_f = per-founder unique-k-mer count (sum of cn rows). Saved BEFORE
    # any normalization, so callers can use it to correct h post-EM.
    K_f = cn_dense.sum(axis=1).astype(np.float32)  # F-vector
    if row_normalize:
        # Equalize per-founder evidence budget: each row sums to 1 after this.
        # Math: replace cn[f,k]=1 with cn[f,k]=1/K_f. Multiplicative-update EM
        # accepts non-binary cn. Recovery for mass-domain projection:
        #   h_proj = (h_em / K_f); h_proj /= h_proj.sum()
        K_f_safe = np.maximum(K_f, np.float32(1.0)).reshape(-1, 1)
        cn_dense = cn_dense / K_f_safe
        kf_cv = float(K_f.std() / max(K_f.mean(), 1.0))
        print(f"  [{chrom}] row-normalized cn_dense; K_f mean={K_f.mean():.0f} "
              f"min={int(K_f.min())} max={int(K_f.max())} CV={kf_cv:.3f}", flush=True)
    return cn_dense, counts, meta, cov, F, K, K_f


# Module-level globals populated by main(): cn_var_called sparse matrix +
# whether it's available. _project_with_called_mask reads these so we don't
# need to plumb the called matrix through every function signature.
_CN_VAR_CALLED = None  # scipy.sparse, founder × variant; 1 if GT != ./.

def _project_with_called_mask(cn_var_chrom, h, idx):
    """Project h through cn_var with per-record renormalization by the called mask.

    AF_est[r] = (h @ cn_var)[r] / (h @ cn_var_called)[r]

    This handles missing GTs correctly: at records where a subset of founders
    is ./., their h-mass is excluded from BOTH the numerator and the denominator.
    Equivalent to AC/AN when h is uniform. Falls back to plain (h @ cn_var) if
    cn_var_called is not available.
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
    # avoid 0/0 for records where no called founder has any h
    safe = np.maximum(called_weight, np.array(1e-12, dtype=freqs.dtype))
    return (freqs / safe).astype(freqs.dtype)


def _apply_kf_correction(h, K_f, alpha: float = 1.0):
    """Per-founder K_f-based correction on EM-domain h.

    Generalizes the rownorm post-correction with a tunable exponent α:
        h_proj[f] ∝ h_em[f] / K_f[f]^α

    α=1.0 (default): theoretical mass-domain inverse — undoes the rownorm
    transform exactly under the Poisson model. This is the established
    "row-norm + recovery" pipeline.
    α=0.0: skip correction — h_proj = h_em (use after rownorm to see the
    untreated EM output).
    α>1.0: stronger penalty on rich-fingerprint founders. Empirically used
    to suppress residual cactus/PG h-bias when rownorm-only leaves a
    structural asymmetry (per-side cactus-only-shared k-mer dominance).

    Works on a 1-D h (global) or 2-D h_blocks (windows). Treats the LAST axis
    as the founder axis. Computation is done in float64 to avoid underflow
    when K_f^α ≫ 1 (e.g. K_f≈2e6, α=5 → K_f^α≈3e31).
    """
    if h is None:
        return h
    h64 = h.astype(np.float64)
    if alpha == 0.0:
        s = h64.sum(axis=-1, keepdims=True)
        out = h64 / np.maximum(s, 1e-300)
        return out.astype(h.dtype)
    K_f_safe = np.maximum(K_f.astype(np.float64), 1.0)
    h_corr = h64 / (K_f_safe ** float(alpha))
    s = h_corr.sum(axis=-1, keepdims=True)
    out = h_corr / np.maximum(s, 1e-300)
    return out.astype(h.dtype)


def run_one_chrom_global(chrom, cn_prefix, cn_var, var_meta, reads_input, threads,
                         em_max_iter=200, row_normalize_cn: bool = False,
                         ac_weight_counts: bool = False,
                         kf_correction_alpha: float = 1.0):
    """Global-mode (single h per chrom) EM + projection.

    ac_weight_counts: when True, multiplies each k-mer's count by its
    carrier count ac_k (number of founders carrying that k-mer in cn) before
    running EM. This cancels the 1/ac_k amplification that singletons get in
    the EM denominator (the "singleton voice" problem). Math:
        em_term[f] = h[f] · Σ_k ac_k · cn[f,k] · c[k] / (h · cn[:,k])
    Equivalent to running standard EM with c_new[k] = c[k] · ac_k.
    """
    t_chrom = time.time()
    cn_dense, counts, meta, cov, F, K, K_f = _count_and_load_cn_dense(
        chrom, cn_prefix, reads_input, threads,
        row_normalize=row_normalize_cn)
    if cn_dense is None:
        return None, None, None, 0.0

    # Filter to nonzero-count k-mers (the EM only needs those)
    nz = counts > 0
    n_nz = int(nz.sum())
    cn_em = np.ascontiguousarray(cn_dense[:, nz])
    counts_em = counts[nz].astype(np.float32)
    del cn_dense
    gc.collect()

    if ac_weight_counts:
        # carrier-weighted EM: scale c[k] by ac_k (number of founders carrying k)
        # This cancels the 1/ac_k amplification singletons get in the M-step.
        ac_k = cn_em.sum(axis=0).astype(np.float32)
        counts_em = counts_em * ac_k
        print(f"  [{chrom}] ac-weighted: counts sum {counts.sum():.0f} -> {counts_em.sum():.0f} "
              f"(ratio {counts_em.sum()/max(counts.sum(),1):.2f}); "
              f"ac_k: min={int(ac_k.min())} median={int(np.median(ac_k))} max={int(ac_k.max())}",
              flush=True)

    t = time.time()
    h, info = solve_em(counts_em, cn_em, cov, max_iter=em_max_iter, tol=1e-7)
    print(f"  [{chrom}] EM solved in {info['iterations']} iters [{time.time()-t:.0f}s]; "
          f"eff_n_founders = {1/np.sum(h**2):.1f}", flush=True)
    del cn_em, counts_em
    gc.collect()

    if row_normalize_cn:
        h_for_proj = _apply_kf_correction(h, K_f, alpha=kf_correction_alpha)
        print(f"  [{chrom}] K_f-corrected h (α={kf_correction_alpha}): "
              f"eff_n_founders = {1/np.sum(h_for_proj**2):.1f} "
              f"(was {1/np.sum(h**2):.1f} pre-correction)",
              flush=True)
    else:
        h_for_proj = h

    rec_chrom = np.asarray(var_meta["chrom"]).astype(str)
    idx = np.where(rec_chrom == str(chrom))[0]
    cn_var_chrom = cn_var[:, idx]
    freqs = _project_with_called_mask(cn_var_chrom, h_for_proj, idx)
    # Per-record info = h-mass on called founders at record r (the denominator
    # of the MAR projection). info ∈ [0, 1]; small info → low-confidence AF.
    if _CN_VAR_CALLED is not None:
        info = np.asarray(_CN_VAR_CALLED[:, idx].T @ h_for_proj).flatten().astype(np.float32)
    else:
        info = np.ones(len(idx), dtype=np.float32)
    elapsed = time.time() - t_chrom
    print(f"  [{chrom}] {elapsed:.0f}s total — {len(idx):,} records projected", flush=True)
    return idx, freqs, info, h_for_proj, elapsed


def run_one_chrom_window(chrom, cn_prefix, cn_var, var_meta, reads_input, threads,
                         window_bp=200_000, em_max_iter=200,
                         block_partition="fixed",
                         r2_threshold=0.5, smooth_window=50,
                         min_block_records=100, max_block_bp=500_000,
                         bigld_panel_npz=None,
                         global_anchor_weight: float = 0.0,
                         window_step: int | None = None,
                         hmm_smooth_passes: int = 0,
                         hmm_smooth_alpha: float = 0.2,
                         hmm_smooth_recomb_rate: float = 4e-8,
                         projection_smooth: str = "auto",
                         row_normalize_cn: bool = False,
                         kf_correction_alpha: float = 1.0):
    """Window/block-mode per-chrom EM + smooth projection.

    block_partition options:
      "fixed":      hard-cut every window_bp (default, HAFpipe-like)
      "ld_gabriel": Gabriel-style adjacent-r² block boundaries (LD-informed)
      "ld_complete": CompleteLDPartition (hapFIRE-style fully-independent blocks)

    LD modes use cn_var (founder × record) to derive blocks where founders are
    in LD; recombination breakpoints between blocks are where founder identity
    can shift independently. Block sizes adapt to the local recombination
    landscape, so they're large where LD spans far (chromosome arms) and short
    where LD breaks (centromeres, hotspots).
    """
    t_chrom = time.time()
    cn_dense, counts, meta, cov, F, K, K_f = _count_and_load_cn_dense(
        chrom, cn_prefix, reads_input, threads,
        row_normalize=row_normalize_cn)
    if cn_dense is None:
        return None, None, None, None, 0.0

    bubble_id = meta["bubble_id"]
    bubble_chrom = meta["bubble_chrom"]
    bubble_start = meta["bubble_start"]
    bubble_end = meta["bubble_end"]

    if block_partition == "fixed":
        blocks = define_windows(bubble_chrom, bubble_start, bubble_end,
                                window_bp=window_bp, window_step=window_step)
        if window_step is not None and window_step < window_bp:
            print(f"  [{chrom}] {len(blocks)} OVERLAPPING windows of {window_bp:,} bp "
                  f"(step={window_step:,} bp, ~{window_bp//window_step}× cover)",
                  flush=True)
        else:
            print(f"  [{chrom}] {len(blocks)} fixed windows of {window_bp:,} bp",
                  flush=True)
    elif block_partition == "bigld_panel":
        # Load pre-computed BigLD blocks from hapfire_block_index.npz, filter
        # to this chrom. Identical block boundaries to hapFIRE's per-fine-block
        # output → the strictest apples-to-apples for "same blocks, different
        # evidence (k-mers vs SNP-haplotype)".
        if bigld_panel_npz is None:
            raise ValueError("block_partition='bigld_panel' requires --bigld-panel-npz")
        blocks = load_bigld_panel_blocks(bigld_panel_npz, chrom_filter=str(chrom))
        sizes_bp = [(b.end - b.start + 1) for b in blocks] if blocks else [0]
        print(f"  [{chrom}] {len(blocks)} BigLD-panel blocks; "
              f"sizes bp: median={int(np.median(sizes_bp)):,}, max={max(sizes_bp):,}", flush=True)
    elif block_partition in ("ld_gabriel", "ld_complete"):
        # Restrict cn_var to this chrom's biallelic records (cn_var rows are founders,
        # cols are records on the projection axis). LD blocks are panel-derived,
        # one-time work per chromosome.
        rec_chrom_all = np.asarray(var_meta["chrom"]).astype(str)
        rec_pos_all = np.asarray(var_meta["pos"]).astype(np.int64)
        rec_idx = np.where(rec_chrom_all == str(chrom))[0]
        cn_var_chrom_records = cn_var[:, rec_idx]
        rec_chrom_c = rec_chrom_all[rec_idx]
        rec_pos_c = rec_pos_all[rec_idx]
        t = time.time()
        if block_partition == "ld_gabriel":
            blocks = compute_ld_blocks_gabriel(
                cn_var_chrom_records, rec_chrom_c, rec_pos_c,
                r2_threshold=r2_threshold, smooth_window=smooth_window,
                min_block_records=min_block_records, max_block_bp=max_block_bp,
                verbose=False,
            )
        else:  # ld_complete
            blocks = compute_ld_blocks(
                cn_var_chrom_records, rec_chrom_c, rec_pos_c,
                r2_threshold=r2_threshold, search_window=100,
                min_block_records=min_block_records, verbose=False,
            )
        sizes_bp = [(b.end - b.start + 1) for b in blocks]
        print(f"  [{chrom}] {len(blocks)} LD blocks ({block_partition}, "
              f"r²={r2_threshold}); block sizes bp: median={int(np.median(sizes_bp)):,}, "
              f"max={max(sizes_bp):,}; {time.time()-t:.0f}s", flush=True)
    else:
        raise ValueError(f"unknown block_partition: {block_partition!r}")
    n_blocks = len(blocks)

    overlapping = (block_partition == "fixed"
                   and window_step is not None and window_step < window_bp)
    if overlapping:
        # Each k-mer can be in multiple windows — pass per-block kmer indices
        block_kmer_idx = assign_kmers_to_blocks_multi(
            bubble_id, bubble_chrom, bubble_start, bubble_end, blocks)
        kmer_block = None
    else:
        kmer_block = assign_kmers_to_blocks(bubble_id, bubble_chrom,
                                             bubble_start, bubble_end, blocks)
        block_kmer_idx = None

    t = time.time()
    h_blocks, status, global_h = solve_em_per_block(
        counts.astype(np.float32), cn_dense, kmer_block, n_blocks,
        cov, em_max_iter=em_max_iter, tol=1e-7,
        min_kmers_per_block=200, verbose=False,
        global_anchor_weight=global_anchor_weight,
        block_kmer_idx_override=block_kmer_idx,
    )
    print(f"  [{chrom}] block-EM {time.time()-t:.0f}s "
          f"({(status==0).sum()}/{n_blocks} local fits, "
          f"{(status==1).sum()} fallbacks)", flush=True)
    del cn_dense
    gc.collect()

    # Optional post-EM HMM smoothing across blocks (Li-Stephens style).
    # Recipe matches clean_smooth's best params on bigld blocks (passes=10, α=0.2).
    if hmm_smooth_passes > 0:
        block_pos_start = np.array([b.start for b in blocks], dtype=np.int64)
        block_pos_end = np.array([b.end for b in blocks], dtype=np.int64)
        n_eco = h_blocks.shape[1]
        # Only smooth blocks with successful local fit (status == 0)
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

    if row_normalize_cn:
        h_blocks_proj = _apply_kf_correction(h_blocks, K_f, alpha=kf_correction_alpha)
        global_h_proj = _apply_kf_correction(global_h, K_f, alpha=kf_correction_alpha)
        eff_pre = 1.0 / np.sum(global_h ** 2) if global_h is not None and global_h.sum() > 0 else float('nan')
        eff_post = 1.0 / np.sum(global_h_proj ** 2) if global_h_proj is not None and global_h_proj.sum() > 0 else float('nan')
        print(f"  [{chrom}] K_f-corrected h_blocks (α={kf_correction_alpha}) for projection; "
              f"global eff_n_founders {eff_pre:.1f} → {eff_post:.1f}", flush=True)
    else:
        h_blocks_proj = h_blocks
        global_h_proj = global_h

    rec_chrom = np.asarray(var_meta["chrom"]).astype(str)
    rec_pos = np.asarray(var_meta["pos"])
    idx = np.where(rec_chrom == str(chrom))[0]
    cn_var_chrom = cn_var[:, idx]

    # MAR projection (2026-05-21): pass cn_var_called sliced to this chrom so
    # window-mode matches global-mode's (h@cn_var)/(h@cn_var_called) semantics.
    cn_var_called_chrom = _CN_VAR_CALLED[:, idx] if _CN_VAR_CALLED is not None else None

    if overlapping:
        # Multi-cover inverse-distance averaging
        freqs = project_blocks_to_records_overlap(
            h_blocks_proj, status, global_h_proj, cn_var_chrom, blocks,
            record_chrom=rec_chrom[idx], record_pos=rec_pos[idx],
            cn_var_called=cn_var_called_chrom,
        )
        # Info = same projection but with cn_var_called as the carrier matrix
        # and no called-mask normalization (so we get the raw h-weighted called
        # mass per record).
        if cn_var_called_chrom is not None:
            info = project_blocks_to_records_overlap(
                h_blocks_proj, status, global_h_proj, cn_var_called_chrom, blocks,
                record_chrom=rec_chrom[idx], record_pos=rec_pos[idx],
                cn_var_called=None,
            ).astype(np.float32)
        else:
            info = np.ones(len(idx), dtype=np.float32)
    else:
        # Decide whether to apply linear-interpolation smoothing at projection.
        # If we already applied HMM smoothing on h_blocks, doing linear
        # interpolation on top double-smooths and over-attenuates local signal.
        if projection_smooth == "auto":
            do_smooth = (hmm_smooth_passes == 0)
        elif projection_smooth in ("on", "true", "1"):
            do_smooth = True
        elif projection_smooth in ("off", "false", "0"):
            do_smooth = False
        else:
            raise ValueError(f"projection_smooth must be auto/on/off, got {projection_smooth!r}")
        if hmm_smooth_passes > 0:
            print(f"  [{chrom}] projection smooth={do_smooth} "
                  f"(HMM smoothing {'already' if hmm_smooth_passes>0 else 'not'} applied)",
                  flush=True)
        rec_block = assign_records_to_blocks(rec_chrom[idx], rec_pos[idx], blocks)
        freqs = project_blocks_to_records(
            h_blocks_proj, status, global_h_proj, cn_var_chrom, rec_block,
            record_chrom=rec_chrom[idx], record_pos=rec_pos[idx],
            blocks=blocks, smooth=do_smooth,
            cn_var_called=cn_var_called_chrom,
        )
        # Info = same projection but with cn_var_called as the carrier matrix
        # and no called-mask normalization.
        if cn_var_called_chrom is not None:
            info = project_blocks_to_records(
                h_blocks_proj, status, global_h_proj, cn_var_called_chrom, rec_block,
                record_chrom=rec_chrom[idx], record_pos=rec_pos[idx],
                blocks=blocks, smooth=do_smooth,
                cn_var_called=None,
            ).astype(np.float32)
        else:
            info = np.ones(len(idx), dtype=np.float32)
    elapsed = time.time() - t_chrom
    print(f"  [{chrom}] {elapsed:.0f}s total — {len(idx):,} records projected", flush=True)
    return idx, freqs, info, (h_blocks_proj, status, global_h_proj, blocks), elapsed


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
    ap.add_argument("--block-mode", default="global",
                    choices=["global", "window", "ld_gabriel", "ld_complete", "bigld_panel"],
                    help="global: one h per chrom (F0 pools). "
                         "window: fixed-bp per-window h. "
                         "ld_gabriel: Gabriel-style adjacent-r² LD blocks (cn_var-derived). "
                         "ld_complete: hapFIRE CompleteLDPartition. "
                         "bigld_panel: use pre-computed BigLD blocks from a panel "
                         "(hapfire_block_index.npz) — strictest head-to-head with hapFIRE.")
    ap.add_argument("--bigld-panel-npz", default=None,
                    help="path to hapfire_block_index.npz for --block-mode bigld_panel")
    ap.add_argument("--window-bp", type=int, default=200_000,
                    help="fixed-window size for --block-mode window (default 200 kb)")
    ap.add_argument("--ld-r2", type=float, default=0.5,
                    help="r² threshold for LD-block partitioning (default 0.5)")
    ap.add_argument("--ld-smooth-window", type=int, default=50,
                    help="smoothing window (records) for ld_gabriel (default 50)")
    ap.add_argument("--ld-min-block-records", type=int, default=100,
                    help="merge blocks with fewer than N records into prior (default 100)")
    ap.add_argument("--ld-max-block-bp", type=int, default=500_000,
                    help="hard cap on LD block size in bp (default 500 kb)")
    ap.add_argument("--global-anchor-weight", type=float, default=0.0,
                    help="λ for per-window EM Dirichlet anchor toward chrom-wide "
                         "h_global. 0 = legacy MLE (default). 0.05–0.5 = mild "
                         "anchor — helps low-evidence windows match the global "
                         "solution while letting high-evidence windows adapt.")
    ap.add_argument("--window-step", type=int, default=None,
                    help="step (bp) between adjacent fixed windows. Default = "
                         "window-bp (disjoint windows). Set < window-bp to use "
                         "OVERLAPPING windows (Route 2): each k-mer enters "
                         "multiple windows, per-record AF is inverse-distance "
                         "averaged over all covering windows. E.g. window-bp "
                         "10000 + window-step 3000 → ~3-4× cover.")
    ap.add_argument("--hmm-smooth-passes", type=int, default=0,
                    help="Post-EM Li-Stephens-style HMM smoothing on h_blocks. "
                         "0 = off (default). clean_smooth's best params: passes=10, α=0.2.")
    ap.add_argument("--hmm-smooth-alpha", type=float, default=0.2,
                    help="α for HMM smoothing: h_b' = α·h_b + (1-α)·neighbor_avg. "
                         "Smaller α = more smoothing. Default 0.2.")
    ap.add_argument("--hmm-smooth-recomb-rate", type=float, default=4e-8,
                    help="Recomb rate per bp for distance-weighted neighbor averaging "
                         "(default 4e-8 = 4 cM/Mb, A. thaliana).")
    ap.add_argument("--projection-smooth", default="auto",
                    choices=["auto", "on", "off"],
                    help="Whether to apply HAFpipe-style linear-interpolation "
                         "smoothing between adjacent windows at projection time. "
                         "'auto' (default): smooth=False if HMM smoothing is on "
                         "(avoids double-smoothing), True otherwise. "
                         "'on' / 'off' force the behaviour.")
    ap.add_argument("--row-normalize-cn", action="store_true",
                    help="Row-normalize cn_full at EM-load time so each founder "
                         "contributes equal total k-mer mass (each row sums to 1). "
                         "Aims to remove the +41%% cactus-vs-PG h-bias caused by "
                         "per-founder unique-k-mer-count asymmetry in cn_full_v3. "
                         "Post-EM, h is divided by K_f and renormalized to recover "
                         "the mass-domain founder fractions before projection "
                         "through cn_var. See BALANCING_KMERS.md approach #1.")
    ap.add_argument("--kf-correction-alpha", type=float, default=1.0,
                    help="Exponent α for the post-EM K_f correction "
                         "h_proj ∝ h_em / K_f^α (only applies with "
                         "--row-normalize-cn). α=1.0 (default) = standard "
                         "rownorm mass-domain recovery. α=0 = skip correction. "
                         "α>1 = stronger penalty on rich-fingerprint founders, "
                         "used to suppress residual cactus/PG h-bias. "
                         "α≈5 expected to ~fully balance per-founder h.")
    ap.add_argument("--ac-weight-counts", action="store_true",
                    help="Carrier-weighted EM: multiply each k-mer count by its "
                         "ac_k (number of founders carrying it) before EM. "
                         "Cancels the 1/ac_k amplification that singletons receive "
                         "in the EM denominator, putting every k-mer on equal "
                         "evidence footing. Analogous to freqk's per-allele "
                         "normalization but applied per k-mer. Currently only "
                         "wired through for --block-mode global.")
    args = ap.parse_args()

    print(f"=== {args.sample} (per-chrom {args.block_mode} mode) ===", flush=True)
    print(f"  reads: {args.reads}", flush=True)
    t0 = time.time()

    cn_var = load_npz(args.cn_var)
    var_meta = np.load(args.cn_var_meta, allow_pickle=True)
    n_records = cn_var.shape[1]
    print(f"  cn_var: {cn_var.shape}, n_records: {n_records:,}", flush=True)

    # Optional: load cn_var_called mask for proper missing-aware AF projection.
    # See SESSION_2026-05-18.md cn_var bug.
    global _CN_VAR_CALLED
    called_path = args.cn_var_called
    if called_path is None:
        # auto-detect: same prefix + .cn_var_called.npz
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
              f"(legacy behavior). Rebuild with the updated build_cn_var.py for "
              f"missing-aware projection.", flush=True)

    reads_input = args.reads if len(args.reads) > 1 else args.reads[0]

    freqs_global = np.full(n_records, np.nan, dtype=np.float32)
    info_global  = np.full(n_records, np.nan, dtype=np.float32)
    h_save = {}  # global mode: h_per_chrom[chrom] = h. window mode: per-chrom block packs.

    # Precompute panel-level per-record called counts (h-independent QC metric).
    # When cn_var_called is unavailable, assume the panel size F (all called).
    F_total = cn_var.shape[0]
    if _CN_VAR_CALLED is not None:
        n_called_per_rec = np.asarray(_CN_VAR_CALLED.sum(axis=0)).flatten().astype(np.int32)
    else:
        n_called_per_rec = np.full(n_records, F_total, dtype=np.int32)

    if args.row_normalize_cn:
        print(f"  --row-normalize-cn ON: each cn_full row will be divided by K_f "
              "before EM; h will be divided by K_f and renormalized before "
              "projection. See BALANCING_KMERS.md approach #1.", flush=True)

    if args.ac_weight_counts:
        print(f"  --ac-weight-counts ON: per-k-mer counts will be multiplied by "
              "their carrier count ac_k before EM. See em_solver.py and "
              "BALANCING_KMERS.md for rationale.", flush=True)

    for chrom in args.chroms:
        if args.block_mode == "global":
            idx, freqs, info, h, _ = run_one_chrom_global(
                chrom, args.cn_kmer_prefix, cn_var, var_meta,
                reads_input, args.threads,
                row_normalize_cn=args.row_normalize_cn,
                ac_weight_counts=args.ac_weight_counts,
                kf_correction_alpha=args.kf_correction_alpha,
            )
            if idx is None:
                continue
            h_save[chrom] = h
        else:
            partition = "fixed" if args.block_mode == "window" else args.block_mode
            idx, freqs, info, pack, _ = run_one_chrom_window(
                chrom, args.cn_kmer_prefix, cn_var, var_meta,
                reads_input, args.threads,
                window_bp=args.window_bp,
                block_partition=partition,
                r2_threshold=args.ld_r2,
                smooth_window=args.ld_smooth_window,
                min_block_records=args.ld_min_block_records,
                max_block_bp=args.ld_max_block_bp,
                bigld_panel_npz=args.bigld_panel_npz,
                global_anchor_weight=args.global_anchor_weight,
                window_step=args.window_step,
                hmm_smooth_passes=args.hmm_smooth_passes,
                hmm_smooth_alpha=args.hmm_smooth_alpha,
                hmm_smooth_recomb_rate=args.hmm_smooth_recomb_rate,
                projection_smooth=args.projection_smooth,
                row_normalize_cn=args.row_normalize_cn,
                kf_correction_alpha=args.kf_correction_alpha,
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

    # SE per record (Wald, using n_called as effective N). Set to NaN if p is
    # NaN or n_called == 0. Downstream consumers can filter/IVW on this.
    # n_called is the integer count of called founders at each record
    # (h-independent panel QC). info is the h-weighted version (varies with EM).
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
        # Schema (2026-05-21): added info, n_called, se columns. Old parsers
        # that read only the first 5 columns still work; new code can read
        # the per-record uncertainty metrics.
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
