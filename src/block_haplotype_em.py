"""Per-block haplotype-level EM driver for cactus_em (bigld_haplotype mode).

Input:
  - block_haplotype_cn_<chrom>.npz (output of build_block_haplotype_cn.py)
    holds per-block (unique-hap × k-mer) presence matrices, founder→hap
    mappings, and kept k-mer strings.
  - reads (FASTQ or BAM)
  - cn_var (founder × record), used for the final AF projection

Inference unit: per BigLD block, ~30-60-dim unique-haplotype simplex.
EM uses k-mer evidence DIRECTLY at the fine-block level — no propagation
from coarser blocks (this is the structural difference vs hapFIRE_perblock,
which propagates HARP estimates from coarse independent LD blocks
analytically into fine blocks under the assumption of intact within-coarse
ancestry).

Output: per-record alt_freq for chrom, plus h_hap arrays for diagnostics.
"""
from __future__ import annotations
import os
import sys
import time
from pathlib import Path
import numpy as np
from scipy.sparse import csr_matrix

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from em_solver import solve_em
from kmer_count import count_kmers_in_bam, count_kmers_in_fasta


def load_block_hap_cn(npz_path: str):
    """Load a chromosome's block_haplotype_cn npz and reconstruct the
    big sparse CSR + per-block slicing arrays.

    Returns dict with:
      ecotypes:           (N_eco,) str
      block_idx:          (n_blocks,) int — original block index in source npz
      block_chrom:        (n_blocks,) str
      block_pos_start:    (n_blocks,) int
      block_pos_end:      (n_blocks,) int
      block_n_haps:       (n_blocks,) int
      block_n_kmers:      (n_blocks,) int
      block_kmer_offsets: (n_blocks+1,) int — k-mer string range per block
      block_hap_offsets:  (n_blocks+1,) int — hap-row range per block
      kmer_strings:       (total_kmers,) str — all kept k-mer strings (concat per block)
      cn:                 sparse CSR of shape (total_haps, total_kmers); rows are
                          block-haps (block b = rows [hap_offsets[b]:hap_offsets[b+1]]),
                          cols are global-kmer indices (block b uses cols
                          [kmer_offsets[b]:kmer_offsets[b+1]]).
      founder_to_class:         (n_blocks, N_eco) int — for each (block, founder) which
                          unique-hap-within-block it maps to.
      k:                  k-mer length
    """
    z = np.load(npz_path, allow_pickle=True)
    ecotypes = np.asarray(z["ecotypes"]).astype(str)
    block_idx = np.asarray(z["block_idx"]).astype(np.int64)
    block_chrom = np.asarray(z["block_chrom"]).astype(str)
    block_pos_start = np.asarray(z["block_pos_start"]).astype(np.int64)
    block_pos_end = np.asarray(z["block_pos_end"]).astype(np.int64)
    block_n_haps = np.asarray(z["block_n_haps"]).astype(np.int32)
    block_n_kmers = np.asarray(z["block_n_kmers"]).astype(np.int32)
    block_kmer_offsets = np.asarray(z["block_kmer_offsets"]).astype(np.int64)
    block_hap_offsets = np.asarray(z["block_hap_offsets"]).astype(np.int64)
    kmer_strings = np.asarray(z["kmer_strings"]).astype(str)
    cn_indptr = np.asarray(z["cn_indptr"]).astype(np.int64)
    cn_indices = np.asarray(z["cn_indices"]).astype(np.int64)
    cn_data = np.asarray(z["cn_data"]).astype(np.uint8)
    founder_to_class = np.asarray(z["founder_to_class"]).astype(np.int32)
    k = int(z["k"])

    total_haps = int(block_hap_offsets[-1])
    total_kmers = int(block_kmer_offsets[-1])
    cn = csr_matrix((cn_data, cn_indices, cn_indptr),
                    shape=(total_haps, total_kmers))

    return dict(
        ecotypes=ecotypes,
        block_idx=block_idx,
        block_chrom=block_chrom,
        block_pos_start=block_pos_start,
        block_pos_end=block_pos_end,
        block_n_haps=block_n_haps,
        block_n_kmers=block_n_kmers,
        block_kmer_offsets=block_kmer_offsets,
        block_hap_offsets=block_hap_offsets,
        kmer_strings=kmer_strings,
        cn=cn,
        founder_to_class=founder_to_class,
        k=k,
    )


def smooth_h_across_blocks(
    h_founder_per_block: dict[int, np.ndarray],
    block_pos_start: np.ndarray,
    block_pos_end: np.ndarray,
    n_eco: int,
    recomb_rate: float = 4e-8,
    alpha: float = 0.7,
    n_passes: int = 3,
):
    """Li-Stephens-style HMM smoothing pass on per-block h_founder estimates.

    For each block b, blend its h with a recomb-distance-weighted average of
    its left and right valid neighbors:

        h_b' = alpha * h_b + (1 - alpha) * weighted_avg(h_neighbors)

    Apply n_passes times. Keeps blocks with sparse evidence (single-founder
    dominant) consistent with neighbors that share ancestry, addressing the
    "winner-takes-all" failure mode in heavy-recombination regimes.

    h_founder_per_block: dict {block_idx: (n_eco,) array}, only for blocks
        with successful EM (others are skipped).
    Returns: a new dict with the same keys, smoothed h vectors.
    """
    if not h_founder_per_block or n_passes <= 0 or alpha >= 1.0:
        return h_founder_per_block

    valid_blocks = sorted(h_founder_per_block.keys())
    n_blocks_total = len(block_pos_start)
    mid_pos = (block_pos_start + block_pos_end) / 2.0

    current = dict(h_founder_per_block)
    for _it in range(n_passes):
        new_h: dict[int, np.ndarray] = {}
        for b in valid_blocks:
            h_b = current[b]
            # Find left/right valid neighbors
            left_b = right_b = None
            for j in range(b - 1, -1, -1):
                if j in current:
                    left_b = j
                    break
            for j in range(b + 1, n_blocks_total):
                if j in current:
                    right_b = j
                    break
            neighbor_avg = None
            if left_b is not None and right_b is not None:
                gap_l = abs(mid_pos[b] - mid_pos[left_b])
                gap_r = abs(mid_pos[b] - mid_pos[right_b])
                w_l = float(np.exp(-recomb_rate * gap_l))
                w_r = float(np.exp(-recomb_rate * gap_r))
                norm = w_l + w_r
                if norm > 0:
                    neighbor_avg = (w_l * current[left_b] + w_r * current[right_b]) / norm
            elif left_b is not None:
                neighbor_avg = current[left_b]
            elif right_b is not None:
                neighbor_avg = current[right_b]
            if neighbor_avg is not None:
                h_smoothed = alpha * h_b + (1.0 - alpha) * neighbor_avg
                # Renormalize to simplex
                s = h_smoothed.sum()
                if s > 0:
                    h_smoothed = h_smoothed / s
                new_h[b] = h_smoothed.astype(np.float32)
            else:
                new_h[b] = h_b
        current = new_h
    return current


def project_h_class_to_founders(h_class: np.ndarray, f2c_b: np.ndarray) -> np.ndarray:
    """Distribute per-class freq equally among founders mapping to that class.

    Within a sequence-unique class, all founders have *identical* block
    sequence — they're statistically interchangeable, so we split h_class[C]
    evenly among them.

    h_class: (n_uniq,) — frequencies on the within-block class simplex.
    f2c_b:   (n_eco,) — founder-to-class mapping for this block.
    Returns: (n_eco,) — per-founder freq, sums to ~1.
    """
    n_uniq = h_class.shape[0]
    n_per_class = np.bincount(f2c_b, minlength=n_uniq).astype(np.float32)
    n_per_class[n_per_class == 0] = 1  # defensive
    h_founder = h_class[f2c_b] / n_per_class[f2c_b]
    return h_founder.astype(np.float32)


def run_chrom_bigld_haplotype(
    chrom: str,
    block_hap_npz_path: str,
    cn_var,                 # sparse (n_eco, n_records) or dense
    var_meta,
    reads_input,
    threads: int = 4,
    em_max_iter: int = 200,
    min_block_kmer_count: int = 30,
    inference_unit: str = "class",
    smooth_passes: int = 0,
    smooth_alpha: float = 0.7,
    smooth_recomb_rate: float = 4e-8,
    dirichlet_alpha: float = 1.0,        # 1.0 = no regularization
    fallback_to_global: bool = False,    # blend per-block h with chrom-global h
    fallback_eff_n_threshold: float = 0.0,  # if eff_n < threshold, blend toward global
):
    """Run bigld_haplotype mode on one chromosome.

    Returns (idx, freqs, h_hap_packs, elapsed)
      idx: (n_records_chrom,) record indices in cn_var
      freqs: (n_records_chrom,) projected alt-freqs
      h_hap_packs: dict with per-block diagnostics
      elapsed: wall time
    """
    t_chrom = time.time()
    print(f"[{chrom}] loading block-hap cn from {block_hap_npz_path}...", flush=True)
    bh = load_block_hap_cn(block_hap_npz_path)
    n_blocks = len(bh["block_idx"])
    n_eco = len(bh["ecotypes"])
    total_kmers = int(bh["block_kmer_offsets"][-1])
    total_haps = int(bh["block_hap_offsets"][-1])
    print(f"  blocks={n_blocks}, n_eco={n_eco}, total_kmers={total_kmers:,}, "
          f"total_haps={total_haps:,}", flush=True)

    # Filter to this chrom (the npz might already be chrom-specific, but be safe)
    keep = bh["block_chrom"] == chrom
    if not keep.all():
        print(f"  WARN: block-hap npz contains chroms {set(bh['block_chrom'])}; "
              f"filtering to {chrom}", flush=True)
        # Build a chrom-filtered view (re-index)
        # For this prototype, expect npz to be chrom-specific. Bail if mixed.
        if not keep.any():
            print(f"  no blocks for {chrom}, skipping", flush=True)
            return None, None, None, 0.0

    # Count k-mers in reads (single jellyfish pass on union)
    t = time.time()
    kmer_list = list(bh["kmer_strings"])
    if isinstance(reads_input, str) and reads_input.endswith(".bam"):
        cd = count_kmers_in_bam(reads_input, kmer_list, k=bh["k"],
                                threads=threads, hash_size="3G")
    else:
        cd = count_kmers_in_fasta(reads_input, kmer_list, k=bh["k"],
                                  threads=threads, hash_size="3G")
    counts = np.array([cd[km] for km in kmer_list], dtype=np.int64)
    nz = (counts > 0).sum()
    print(f"  counted {len(kmer_list):,} k-mers in {time.time()-t:.0f}s; "
          f"nonzero {nz:,} ({nz/max(1,len(kmer_list))*100:.1f}%)", flush=True)

    # Estimate global lambda for diagnostics (cn over all blocks gives an estimate)
    cn_all = bh["cn"]
    ac_total = float(cn_all.sum())
    cov_est = float(counts.sum()) * n_eco / max(1, ac_total) if ac_total > 0 else 0.0
    print(f"  cov estimate: {cov_est:.1f}×", flush=True)

    # Per-block EM
    # inference_unit="class": EM on n_uniq[b]-dim simplex (haplotype-class), then
    #                        scatter to founders by founder_to_class.
    # inference_unit="founder": scatter cn from class to founder rows BEFORE EM,
    #                          then EM on 231-founder simplex directly. Founders
    #                          in the same class share identical cn rows (so
    #                          the EM can swap mass between them freely; we do
    #                          NOT split — accept multi-modal solutions).
    t = time.time()
    h_per_block: list[np.ndarray | None] = [None] * n_blocks  # h_class or h_founder
    block_status = np.zeros(n_blocks, dtype=np.int8)
    n_skipped = 0
    n_converged = 0
    for b in range(n_blocks):
        n_haps_b = int(bh["block_n_haps"][b])
        ks = int(bh["block_kmer_offsets"][b])
        ke = int(bh["block_kmer_offsets"][b + 1])
        hs = int(bh["block_hap_offsets"][b])
        he = int(bh["block_hap_offsets"][b + 1])
        if ke == ks or n_haps_b <= 1:
            block_status[b] = 1
            n_skipped += 1
            continue
        cn_b_class = cn_all[hs:he, ks:ke].toarray().astype(np.float32)  # (n_uniq, n_kmers)
        c_b = counts[ks:ke].astype(np.float32)
        if c_b.sum() < min_block_kmer_count:
            block_status[b] = 1
            n_skipped += 1
            continue

        # Build the cn matrix the EM operates on
        if inference_unit == "class":
            cn_for_em = cn_b_class  # (n_uniq, n_kmers)
            n_dim = n_haps_b
        else:  # "founder"
            f2c_b = bh["founder_to_class"][b]  # (n_eco,)
            cn_for_em = cn_b_class[f2c_b]      # (n_eco, n_kmers) — scatter rows
            n_dim = n_eco

        # Filter to nonzero-count k-mers (EM only needs those)
        nzm = c_b > 0
        cn_em = np.ascontiguousarray(cn_for_em[:, nzm])
        c_em = c_b[nzm]
        if cn_em.shape[1] < 1:
            # No nonzero k-mer evidence; skip
            block_status[b] = 1
            n_skipped += 1
            continue
        h_b, info = solve_em(c_em, cn_em, cov_est, max_iter=em_max_iter, tol=1e-7,
                             dirichlet_alpha=dirichlet_alpha)
        h_per_block[b] = h_b.astype(np.float32)
        block_status[b] = 0
        n_converged += 1

    print(f"  per-block EM in {time.time()-t:.0f}s: {n_converged} converged, "
          f"{n_skipped} skipped (too few k-mers/counts/under-determined)", flush=True)

    # Build per-block h_founder for projection.
    # For skipped blocks, fall back to a uniform 1/N h_hap (i.e. uniform founder freq within block).
    # (Could also compute a chrom-level global h as fallback; for now keep simple.)
    rec_chrom_all = np.asarray(var_meta["chrom"]).astype(str)
    rec_pos_all = np.asarray(var_meta["pos"]).astype(np.int64)
    chrom_mask = rec_chrom_all == str(chrom)
    rec_idx = np.where(chrom_mask)[0]
    rec_pos = rec_pos_all[rec_idx]
    cn_var_chrom = cn_var[:, rec_idx]
    if hasattr(cn_var_chrom, "toarray"):
        # Keep sparse for memory; we'll do per-block dense projections
        pass

    # Assign each record to a block by position (binary search)
    block_starts = bh["block_pos_start"]
    block_ends = bh["block_pos_end"]
    # Sorted by block_pos_start (should be from BigLD partition; verify)
    order = np.argsort(block_starts)
    bs = block_starts[order]
    be = block_ends[order]
    inv_order = np.empty_like(order)
    inv_order[order] = np.arange(len(order))
    rec_block = np.full(len(rec_pos), -1, dtype=np.int64)
    # For each record, find first block with block_start <= pos <= block_end
    ins = np.searchsorted(bs, rec_pos, side="right") - 1
    in_range = (ins >= 0) & (rec_pos <= be[np.clip(ins, 0, len(be) - 1)])
    rec_block[in_range] = order[ins[in_range]]

    # Project: per-record AF
    freqs = np.full(len(rec_pos), np.nan, dtype=np.float32)
    n_fallback = 0
    # Pre-compute h_founder per block (only for converged blocks)
    h_founder_per_block: dict[int, np.ndarray] = {}
    for b in range(n_blocks):
        if block_status[b] != 0:
            continue
        h_b = h_per_block[b]
        if inference_unit == "class":
            f2c_b = bh["founder_to_class"][b]
            h_founder = project_h_class_to_founders(h_b, f2c_b)
        else:  # founder mode — h_b is already on the 231-founder simplex
            h_founder = h_b.astype(np.float32)
        h_founder_per_block[b] = h_founder

    # Optional confidence-weighted blend with chrom-global h.
    # Compute global h from ALL blocks' counts (sum over blocks → single EM).
    if fallback_to_global and len(h_founder_per_block) > 0:
        t = time.time()
        # Aggregate global counts per kmer (rather than re-EM): just use mean of
        # per-block h_founder weighted by block reads.
        # Simpler: use mean of converged-block h as global proxy.
        global_h = np.mean(np.stack(list(h_founder_per_block.values())), axis=0)
        global_h = global_h / max(global_h.sum(), 1e-9)
        # Blend each block's h toward global, weighted by 1/eff_n
        # eff_n = 1/Σh² (low value → highly concentrated → uncertain on the other founders)
        # Use: blend_w = exp(-eff_n / threshold) when eff_n < threshold, else 0
        n_blended = 0
        for b in list(h_founder_per_block.keys()):
            h_b = h_founder_per_block[b]
            eff_n = 1.0 / max(np.sum(h_b ** 2), 1e-9)
            if fallback_eff_n_threshold > 0 and eff_n < fallback_eff_n_threshold:
                # blend toward global: more blend if eff_n is small
                w_local = float(eff_n / fallback_eff_n_threshold)
                w_global = 1.0 - w_local
                h_blended = w_local * h_b + w_global * global_h
                s = h_blended.sum()
                if s > 0:
                    h_blended /= s
                h_founder_per_block[b] = h_blended.astype(np.float32)
                n_blended += 1
        print(f"  blended {n_blended}/{len(h_founder_per_block)} low-confidence blocks "
              f"(eff_n threshold={fallback_eff_n_threshold}) [{time.time()-t:.0f}s]",
              flush=True)

    # Optional Li-Stephens-style smoothing across blocks
    if smooth_passes > 0 and len(h_founder_per_block) > 1:
        t = time.time()
        h_founder_per_block = smooth_h_across_blocks(
            h_founder_per_block,
            np.asarray(bh["block_pos_start"]).astype(np.float64),
            np.asarray(bh["block_pos_end"]).astype(np.float64),
            n_eco,
            recomb_rate=smooth_recomb_rate,
            alpha=smooth_alpha,
            n_passes=smooth_passes,
        )
        print(f"  smoothed h across blocks: passes={smooth_passes}, "
              f"alpha={smooth_alpha}, recomb_rate={smooth_recomb_rate} "
              f"[{time.time()-t:.0f}s]", flush=True)

    cn_var_dense_chrom = (cn_var_chrom.toarray() if hasattr(cn_var_chrom, "toarray")
                          else np.asarray(cn_var_chrom)).astype(np.float32)

    # For records whose block has converged, project via that block's h_founder.
    # For records out of any block, or in a skipped block, fall back to uniform 1/N
    uniform_h = np.full(n_eco, 1.0 / n_eco, dtype=np.float32)
    for r_local in range(len(rec_pos)):
        b = int(rec_block[r_local])
        if b == -1 or block_status[b] != 0:
            n_fallback += 1
            freqs[r_local] = float(uniform_h @ cn_var_dense_chrom[:, r_local])
        else:
            h_f = h_founder_per_block[b]
            freqs[r_local] = float(h_f @ cn_var_dense_chrom[:, r_local])

    print(f"  projected {len(rec_pos):,} records "
          f"({n_fallback:,} fallback to uniform)", flush=True)

    elapsed = time.time() - t_chrom
    print(f"[{chrom}] {elapsed:.0f}s total — "
          f"{len(rec_pos):,} records projected", flush=True)

    # Build a list aligned to n_blocks for the post-processed founder-level h
    # (reflects fallback-blend + smoothing if applied; same as raw when neither).
    h_per_block_post: list[np.ndarray | None] = [None] * n_blocks
    for b, h in h_founder_per_block.items():
        h_per_block_post[b] = h

    return rec_idx, freqs, dict(
        block_status=block_status,
        h_per_block=h_per_block,
        h_per_block_post=h_per_block_post,
        block_chrom=bh["block_chrom"],
        block_pos_start=bh["block_pos_start"],
        block_pos_end=bh["block_pos_end"],
    ), elapsed
