"""
Per-chromosome kMate driver — founder-mixture (h) estimation + AF projection.

ONE estimator, selected by the estimation UNIT (`--unit`). Every unit is fit the
same way: unit → distinct-haplotype collapse (K_b ≤ F, --haploblock-eps 0 = exact,
a no-op when all founders are distinct) → local EM → equal-split back to members →
project (see ALGORITHM.md §4.4). Two internal fit helpers back the units — they
differ only because the whole-chromosome one-h estimand supports extra outputs:
  * `_fit_unit_chrom`   — --unit chrom: a single h over the whole chromosome
        (selfing / inbred / F0 pools, e.g. SEEDMIX). The only unit that supports
        --h-only and --emit-af-se. (Deprecated alias: --block-mode global.)
  * `_fit_unit_blocks`  — --unit {ld,bp,tsv}: an independent h per block, LOCAL-ONLY
        (no anchor prior, no fallback → thin/empty blocks give NaN AF, no smoothing),
        projected per block. `ld` = r²-LD CompleteLDPartition blocks (production);
        `bp` = fixed --window-bp windows; `tsv` = explicit --blocks-tsv. (Deprecated
        alias: --block-mode window → bp.) Pass --no-local-only for the legacy
        anchored+smoothed "star2" recipe.

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
from .em_solver import solve_em, haploblock_collapse_indices
from .kmer_count import count_kmers_in_bam, count_kmers_in_fasta, query_kmer_db
from .block_em import (define_windows, assign_kmers_to_blocks,
                      assign_records_to_blocks, solve_em_per_block,
                      project_blocks_to_records, BlockSpec)
from .block_haplotype_em import smooth_h_across_blocks
from .h_uncertainty import (fisher_information_h, _tangent_pinv_on_support,
                            af_se_from_cov)

# Calibrated AF identifiability floor c for the total SE  sqrt(SE_Fisher^2 + c^2).
# The per-record AF error is dominated by a coverage-INDEPENDENT non-identifiability
# bias (collinear founders) that no variance term sees, so the Fisher delta-method
# SE alone under-covers (a bare 95% Fisher interval covers truth only ~12-19%); c
# restores ~95% interval coverage. Panel-specific — recompute with
# benchmarks/h_uncertainty/af_calibrate_floor.py; override with --af-id-floor.
# Calibrated on the 231 arch3 panel (closed-loop g0, cov 10/30) to cover ALL
# defined records at 95% (c=0.0186). Note only ~33% of 231 records are well-called;
# the well-called subset alone calibrates to a smaller c≈0.0102, but using that would
# under-cover the 67% partially-called records (cactus-only SVs), so we take the
# all-defined value (also matches p80's 0.018 — the floor is a stable panel property).
# NOTE (2026-07-07): this c was calibrated under the legacy "global" normalization.
# After the per_founder (Kf_w) fix it should be RE-calibrated on per_founder outputs
# (benchmarks/h_uncertainty/af_calibrate_floor.py) as part of the post-rerun refresh;
# left at 0.0186 until that calibration is actually run. Override with --af-id-floor.
AF_ID_FLOOR_DEFAULT = 0.0186


def _resolvability_from_J(J, support):
    """eff_rank (exp spectral entropy) + condition number of the support Fisher info
    — the h-resolvability diagnostic (low eff_rank / high cond ⇒ collinear founders
    recoverable only in aggregate; report instead of a per-founder h SE)."""
    s = np.asarray(support, int)
    if s.size <= 1:
        return float(s.size), float("inf")
    ev = np.clip(np.linalg.eigvalsh(J[np.ix_(s, s)]), 0, None)
    pos = ev[ev > ev.max() * 1e-12] if ev.max() > 0 else ev[:0]
    if pos.size == 0:
        return 0.0, float("inf")
    p = pos / pos.sum()
    return float(np.exp(-(p * np.log(p)).sum())), float(pos.max() / pos.min())


def _af_se_total_chunked(h, Sigma, support, var_pa_chrom, var_called_chrom,
                         floor, chunk_r=200_000):
    """Per-record total AF SE = sqrt(Fisher delta-method SE^2 + floor^2), computed in
    record blocks so the dense F×R_chunk var slices stay bounded in memory."""
    R = var_pa_chrom.shape[1]
    se = np.empty(R, dtype=np.float32)
    for s0 in range(0, R, chunk_r):
        sl = slice(s0, min(s0 + chunk_r, R))
        Vc = var_pa_chrom[:, sl]
        V = Vc.toarray().astype(np.float64) if hasattr(Vc, "toarray") else np.asarray(Vc, np.float64)
        if var_called_chrom is not None:
            Uc = var_called_chrom[:, sl]
            U = Uc.toarray().astype(np.float64) if hasattr(Uc, "toarray") else np.asarray(Uc, np.float64)
        else:
            U = np.ones_like(V)
        _, se_f = af_se_from_cov(h, Sigma, V, U, support=support)
        se[sl] = np.sqrt(se_f * se_f + floor * floor).astype(np.float32)
    return se


def _count_and_load_kmer_pa_dense(chrom, kmer_pa_prefix, reads_input, threads,
                                  kmer_db=None, hash_size="3G",
                                  max_kmer_cov_mult=0.0):
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

    max_kmer_cov_mult: repeat-contamination guard. A panel k-mer is allele-unique
    and (by index construction) panel-unique, so a single-copy k-mer can be hit at
    most ~`cov` times (allele frequency <= 1). A k-mer observed > max_kmer_cov_mult
    × cov is therefore hitting extra genomic copies it doesn't tag (panel-unique
    != genome-unique) — its count is zeroed so the downstream `counts > 0` filter
    DROPS it from the EM (the allele is still tagged by its other unique k-mers).
    cov is re-estimated on the cleaned counts. 0.0 (default) = off / legacy.
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
    if max_kmer_cov_mult and max_kmer_cov_mult > 0:
        thr = max_kmer_cov_mult * cov
        repeat = counts > thr
        n_rep = int(repeat.sum())
        if n_rep:
            counts[repeat] = 0                       # -> dropped by the counts>0 filter
            cov = counts.sum() * F / max(1, ac.sum())  # re-estimate on cleaned counts
        print(f"  [{chrom}] repeat-guard: zeroed {n_rep:,} k-mers with count > "
              f"{max_kmer_cov_mult:g}×cov (={thr:.0f}); cov re-est {cov:.1f}×", flush=True)
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


def _fit_unit_chrom(chrom, kmer_pa_prefix, var_pa, var_meta, reads_input, threads,
                         em_max_iter=200, kmer_weight="uniform", kmer_db=None,
                         hash_size="3G", h_only=False, max_kmer_cov_mult=0.0,
                         emit_af_se=False, af_id_floor=0.0, normalize="per_founder",
                         haploblock_eps=0.0):
    """Fit --unit chrom: one founder mixture h over the whole chromosome + AF
    projection. (Also reached via the deprecated --block-mode global alias.)

    kmer_weight: "uniform" (ω_k=1, MLE) or "inv_mb" (ω_k=1/m_b per-bubble
    de-replication; m_b = #k-mers sharing the k-mer's bubble_id).
    kmer_db: optional prebuilt Jellyfish DB (count-once path; see
    _count_and_load_kmer_pa_dense).
    emit_af_se: also compute the calibrated per-record AF SE (Fisher delta-method
    + identifiability floor) and an h-resolvability diagnostic; returned in `extra`.

    Returns (idx, freqs, info, h, elapsed, extra) where extra is None unless
    emit_af_se, else dict(af_se, eff_rank, cond, support_size).
    """
    t_chrom = time.time()
    kmer_pa_dense, counts, meta, cov, F, K = _count_and_load_kmer_pa_dense(
        chrom, kmer_pa_prefix, reads_input, threads, kmer_db=kmer_db,
        hash_size=hash_size, max_kmer_cov_mult=max_kmer_cov_mult)
    if kmer_pa_dense is None:
        return None, None, None, None, 0.0, None

    # Per-k-mer weight ω_k = 1/m_b (or uniform) — computed over ALL k-mers first.
    omega_full = None
    if kmer_weight == "inv_mb":
        bid = np.asarray(meta["bubble_id"]).astype(np.int64)
        m_b = np.bincount(bid)[bid].astype(np.float32)        # #k-mers per bubble
        omega_full = (1.0 / m_b).astype(np.float32)
    # per_founder normalizer over the FULL panel (all k-mers, incl. c_k=0) — must be
    # computed BEFORE the nonzero filter so it is not conditioned on which k-mers got
    # reads this run (that survivorship bias reintroduces founder collapse). See em_solver.
    w_full = np.ones(kmer_pa_dense.shape[1], np.float32) if omega_full is None else omega_full
    kfw_full = (kmer_pa_dense @ w_full).astype(np.float32)

    # kMate design (block → haploblock → EM, here the "block" is the whole chromosome):
    # COMPUTE the distinct haplotypes the panel resolves (eps=0 exact k-mer identity)
    # instead of ASSUMING all F founders are separately identifiable. Fit the EM over
    # the K_b ≤ F distinct haplotypes, then split each class frequency equally to its
    # member founders. K_b == F (the usual whole-chromosome case) is an EXACT no-op —
    # identical numbers, but data-derived rather than assumed.
    F_founders = kmer_pa_dense.shape[0]
    lab, reps, csize, Kb = haploblock_collapse_indices(kmer_pa_dense, eps=haploblock_eps)
    print(f"  [{chrom}] haploblocks: {Kb} distinct of {F_founders} founders "
          f"({'no collapse (K_b=F)' if Kb == F_founders else f'{F_founders-Kb} merged'})",
          flush=True)

    # Filter to nonzero-count k-mers (the EM only needs those)
    nz = counts > 0
    kmer_pa_em = np.ascontiguousarray(kmer_pa_dense[:, nz])
    counts_em = counts[nz].astype(np.float32)
    del kmer_pa_dense
    gc.collect()

    omega = None if omega_full is None else omega_full[nz]
    if omega is not None:
        print(f"  [{chrom}] ω_k=1/m_b weighting: m_b median={np.median(1.0/omega):.0f} "
              f"max={(1.0/omega).max():.0f}", flush=True)

    t = time.time()
    h_c = None
    if Kb == F_founders:
        # fast no-op path: no distinct-haplotype merging, fit founders directly
        # (byte-identical to the collapse path by permutation-equivariance).
        h, info = solve_em(counts_em, kmer_pa_em, cov, max_iter=em_max_iter, tol=1e-7,
                           omega=omega, normalize=normalize, kfw=kfw_full)
    else:
        h_c, info = solve_em(counts_em, kmer_pa_em[reps], cov, max_iter=em_max_iter, tol=1e-7,
                             omega=omega, normalize=normalize, kfw=kfw_full[reps])
        h = (h_c / csize)[lab].astype(np.float32)     # split haplotype freq to member founders
    print(f"  [{chrom}] EM solved in {info['iterations']} iters [{time.time()-t:.0f}s]; "
          f"eff_n_founders = {1/np.sum(h**2):.1f}", flush=True)

    # AF-certainty: observed Fisher info at ĥ (while the k-mer matrix is still alive).
    Sigma = support = af_diag = None
    if emit_af_se and not h_only:
        t_se = time.time()
        chunk = 1_000_000 if F > 120 else None
        if Kb == F_founders:
            # no collapse: uncertainty on the founders directly (unchanged path).
            J = fisher_information_h(h, kmer_pa_em, counts_em.astype(np.float64),
                                     omega=omega, chunk=chunk)
            support = np.flatnonzero(h > 1e-3)
            Sigma = _tangent_pinv_on_support(J, support, rcond=1e-2)
            eff_rank, cond = _resolvability_from_J(J, support)
            support_size = int(support.size)
        else:
            # COLLAPSED estimator: the point estimate is h = A·h_c with A[f,c] =
            # 1[lab[f]==c]/csize[c], so the correct sampling covariance is on the K_b
            # haplotype classes (h_c), mapped to founder space by Cov(h)=A·Σ_c·Aᵀ.
            # Computing the Fisher info on the full F founders instead would be singular
            # across k-mer-identical class members (the flat ridge the collapse removes).
            Jc = fisher_information_h(h_c, kmer_pa_em[reps], counts_em.astype(np.float64),
                                      omega=omega, chunk=chunk)
            support_c = np.flatnonzero(h_c > 1e-3)
            Sigma_c = _tangent_pinv_on_support(Jc, support_c, rcond=1e-2)
            eff_rank, cond = _resolvability_from_J(Jc, support_c)   # over the K_b classes
            support_size = int(support_c.size)
            A = np.zeros((F_founders, Kb), dtype=np.float64)
            A[np.arange(F_founders), lab] = 1.0 / csize[lab]
            Sigma = A @ Sigma_c @ A.T
            support = np.flatnonzero(h > 1e-3)                     # founder-space support
            J = Jc
        af_diag = dict(eff_rank=eff_rank, cond=cond, support_size=support_size)
        print(f"  [{chrom}] resolvability: eff_rank={eff_rank:.1f} of {support_size} "
              f"support {'founders' if Kb == F_founders else 'haploblocks'}, "
              f"cond={cond:.0f} [{time.time()-t_se:.0f}s]", flush=True)
        del J

    del kmer_pa_em, counts_em
    gc.collect()

    if h_only:
        # MOI / founder-mixture use: skip the AF projection (no var_pa load, no
        # per-record output) — we only need h.
        print(f"  [{chrom}] {time.time()-t_chrom:.0f}s total — h-only (projection skipped)",
              flush=True)
        return None, None, None, h, time.time() - t_chrom, None

    rec_chrom = np.asarray(var_meta["chrom"]).astype(str)
    idx = np.where(rec_chrom == str(chrom))[0]
    var_pa_chrom = var_pa[:, idx]
    freqs, info = _project_with_called_mask(var_pa_chrom, h, idx)
    info = info.astype(np.float32)

    extra = None
    if emit_af_se and Sigma is not None:
        var_called_chrom = _VAR_CALLED[:, idx] if _VAR_CALLED is not None else None
        af_se = _af_se_total_chunked(h, Sigma, support, var_pa_chrom,
                                     var_called_chrom, af_id_floor)
        extra = dict(af_se=af_se, **af_diag)

    elapsed = time.time() - t_chrom
    print(f"  [{chrom}] {elapsed:.0f}s total — {len(idx):,} records projected", flush=True)
    return idx, freqs, info, h, elapsed, extra


def load_blocks_tsv(path, chrom):
    """BlockSpec list for one chrom from an LD-block TSV
    (cols: chrom start_pos end_pos n_variants; TAIR10 1-based)."""
    import pandas as pd
    bt = pd.read_csv(path, sep="\t")
    bt = bt[bt["chrom"].astype(str) == str(chrom)]
    return [BlockSpec(chrom=str(chrom), start=int(s), end=int(e))
            for s, e in zip(bt["start_pos"], bt["end_pos"])]


def _fit_unit_blocks(chrom, kmer_pa_prefix, var_pa, var_meta, reads_input, threads,
                         window_bp=10_000, em_max_iter=200,
                         global_anchor_weight=0.3,
                         hmm_smooth_passes=5,
                         hmm_smooth_alpha=0.5,
                         hmm_smooth_recomb_rate=4e-8,
                         kmer_weight="uniform", kmer_db=None, hash_size="3G",
                         blocks_tsv=None, min_kmers_per_block=200,
                         local_only=True, max_kmer_cov_mult=0.0, normalize="per_founder",
                         haploblock_eps=0.0):
    """Fit --unit {ld,bp,tsv}: an independent h per block + per-block AF projection.
    (Also reached via the deprecated --block-mode window alias → bp.)

    Each window is fit independently: unit → haploblock collapse → per-window EM,
    then projected through var_pa with the missing-aware (called-mask) normalization.

    local_only (default True): the PRODUCTION window recipe — GLOBAL-FREE. No global
    anchor prior, no fallback to the chrom-wide h (thin/empty windows → NaN AF), and
    NO cross-window smoothing. A recombinant pool carries only a handful of
    haplotypes per window; fitting each window purely on its own k-mers (over the K_b
    haplotypes present) is what the data supports, and the chrom-wide mixture is the
    wrong prior for a local window. When local_only, global_anchor_weight and
    hmm_smooth_passes are FORCED to 0 regardless of the values passed.

    Set local_only=False to restore the legacy "star2" recipe (per-window EM anchored
    toward the chrom-wide h_global by global_anchor_weight, then Li-Stephens-style
    smoothed across windows by hmm_smooth_*). NOTE that under a genuine within-window
    haploblock collapse the anchor is split equally across k-mer-indistinguishable
    class members (see block_em._fit_one), so the anchored recipe cannot preserve
    intra-class prior asymmetry — another reason local-only is the default.
    """
    t_chrom = time.time()
    if local_only:
        if global_anchor_weight != 0.0:
            print(f"  [{chrom}] local-only: forcing global_anchor_weight "
                  f"{global_anchor_weight} → 0 (no global prior)", flush=True)
            global_anchor_weight = 0.0
        if hmm_smooth_passes != 0:
            print(f"  [{chrom}] local-only: forcing hmm_smooth_passes "
                  f"{hmm_smooth_passes} → 0 (no cross-window smoothing)", flush=True)
            hmm_smooth_passes = 0
    kmer_pa_dense, counts, meta, cov, F, K = _count_and_load_kmer_pa_dense(
        chrom, kmer_pa_prefix, reads_input, threads, kmer_db=kmer_db,
        hash_size=hash_size, max_kmer_cov_mult=max_kmer_cov_mult)
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
        omega=omega, local_only=local_only, normalize=normalize,
        haploblock_eps=haploblock_eps,
    )
    _low_label = "NaN'd (local-only)" if local_only else "fallbacks"
    print(f"  [{chrom}] block-EM {time.time()-t:.0f}s "
          f"({(status==0).sum()}/{n_blocks} local fits, "
          f"{(status>=1).sum()} {_low_label})", flush=True)
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
    # Global-free: records with no window (rec_block == -1) must NOT borrow
    # global_h. Feed a NaN vector as the projection fallback so they → NaN AF.
    # NaN h_blocks (low/empty windows) already propagate to NaN through h@var_pa.
    proj_fallback_h = np.full_like(global_h, np.nan) if local_only else global_h
    freqs, info = project_blocks_to_records(
        h_blocks, proj_fallback_h, var_pa_chrom, rec_block,
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
    ap.add_argument("--normalize", default="per_founder", choices=["per_founder", "global"],
                    help="EM M-step normalization. 'per_founder' (DEFAULT) divides each founder's "
                         "update by its own observed k-mer content Kf_w (RNA-seq effective-length "
                         "correction) — removes the completeness bias that otherwise collapses "
                         "k-mer-poor founders to ~0. 'global' is the LEGACY multinomial "
                         "normalization by the global count total (kept for reproducing old runs; "
                         "under-calls founder frequencies). Applied in both global and window modes.")
    ap.add_argument("--unit", default=None,
                    choices=["chrom", "ld", "bp", "tsv"],
                    help="Estimation UNIT (one estimator; the mode IS the unit). "
                         "'chrom' (DEFAULT): one h per chromosome — the production "
                         "estimator for selfing / inbred / F0 pools (e.g. GrENE-Net); "
                         "robust on both uniform and sparse panels; supports --h-only "
                         "and --emit-af-se. 'ld': r²-LD blocks from var_pa "
                         "(CompleteLDPartition, --ld-r2) — per-block resolution, but "
                         "collapses in low-diversity blocks (centromere), so it is WRONG "
                         "for selfing pools (see docs/EM_UNIT_CHOICE_AND_NONIDENTIFIABILITY.md). "
                         "'bp': fixed --window-bp windows. 'tsv': explicit --blocks-tsv. "
                         "Each unit is fit locally: haploblock-collapse (--haploblock-eps) "
                         "→ EM → project (no anchor / no smoothing / no fallback).")
    ap.add_argument("--ld-r2", type=float, default=0.1,
                    help="r² cutoff for --unit ld CompleteLDPartition (default 0.1).")
    ap.add_argument("--ld-blocks", default=None,
                    help="Precomputed LD-blocks TSV for --unit ld (chrom start_pos end_pos "
                         "n_variants). If omitted, blocks are computed from --var-pa at "
                         "--ld-r2 and cached next to it (LD blocks are a panel property).")
    ap.add_argument("--ld-window", type=int, default=100,
                    help="CompleteLDPartition window (markers) for --unit ld (default 100).")
    ap.add_argument("--block-mode", default=None, choices=["global", "window"],
                    help="DEPRECATED alias for --unit (kept for back-compat): "
                         "global→--unit chrom, window→--unit bp. Prefer --unit.")
    ap.add_argument("--window-bp", type=int, default=10_000,
                    help="fixed-window size for --unit bp (production window: 10 kb)")
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
                         "h_global. ONLY applies with --no-local-only (the legacy "
                         "'star2' recipe); local-only (the default) forces it to 0.")
    ap.add_argument("--local-only", action=argparse.BooleanOptionalAction, default=True,
                    help="GLOBAL-FREE window mode (DEFAULT, production recipe for "
                         "recombinant pools): no global anchor prior, no global "
                         "fallback, and NO cross-window smoothing. Forces "
                         "--global-anchor-weight and --hmm-smooth-passes to 0; windows "
                         "below --min-kmers-per-block (and empty windows / unassigned "
                         "records) get NaN AF instead of the chrom-wide h. global_h is "
                         "still stored for diagnostics only. Pass --no-local-only to "
                         "restore the legacy anchored+smoothed 'star2' recipe. "
                         "(Ignored in --block-mode global.)")
    ap.add_argument("--hmm-smooth-passes", type=int, default=5,
                    help="Post-EM Li-Stephens-style smoothing passes on h_blocks. "
                         "ONLY applies with --no-local-only; local-only (the default) "
                         "forces it to 0. Legacy 'star2': 5.")
    ap.add_argument("--hmm-smooth-alpha", type=float, default=0.5,
                    help="α for HMM smoothing: h_b' = α·h_b + (1-α)·neighbor_avg. "
                         "Smaller = more smoothing. Production: 0.5.")
    ap.add_argument("--hmm-smooth-recomb-rate", type=float, default=4e-8,
                    help="Recomb rate per bp for distance-weighted neighbor averaging "
                         "(default 4e-8 = 4 cM/Mb, A. thaliana).")
    ap.add_argument("--max-kmer-cov-mult", type=float, default=5.0,
                    help="Repeat-contamination guard (ON by default). Drop any k-mer "
                         "observed > this multiple of the sample's estimated coverage: a "
                         "single-copy, allele-unique k-mer caps at ~cov (allele freq ≤1), "
                         "so a far-higher count means it recurs in a genomic repeat the "
                         "panel didn't model (panel-unique ≠ genome-unique). Such k-mers "
                         "are excluded from the EM and cov is re-estimated. 0 = off "
                         "(legacy / byte-identical to pre-guard runs).")
    ap.add_argument("--emit-af-se", action="store_true",
                    help="GLOBAL mode only: replace the legacy binomial `se` column with "
                         "a calibrated AF SE = sqrt(Fisher delta-method SE^2 + c^2), where "
                         "c (--af-id-floor) is the panel identifiability floor — the per-AF "
                         "error is bias-dominated, so the Fisher SE alone under-covers. Also "
                         "stores an h-resolvability diagnostic (eff_rank/cond) in the h npz. "
                         "Adds an observed-Fisher-info compute per chrom.")
    ap.add_argument("--af-id-floor", type=float, default=AF_ID_FLOOR_DEFAULT,
                    help=f"Identifiability floor c for --emit-af-se (default "
                         f"{AF_ID_FLOOR_DEFAULT}, calibrated on the 231 arch3 panel; "
                         f"recompute per panel with af_calibrate_floor.py).")
    ap.add_argument("--haploblock-eps", type=float, default=0.0,
                    help="Haplotype-merge tolerance for the block→haploblock→EM collapse, "
                         "as a FRACTION of a unit's k-mers. 0.0 (default) = exact "
                         "byte-identical k-mer presence patterns are one haplotype "
                         "(lossless; K_b==F is an exact no-op). >0 = also merge founders "
                         "whose presence differs in ≤ eps·(unit k-mers) — APPROXIMATE, "
                         "merges near-indistinguishable founders (both modes).")
    args = ap.parse_args()

    # Resolve the estimation UNIT. --unit is authoritative; --block-mode is a
    # deprecated alias (global→chrom, window→bp); default is chrom — the production
    # estimator for selfing/inbred pools (ld collapses low-diversity blocks; see
    # docs/EM_UNIT_CHOICE_AND_NONIDENTIFIABILITY.md).
    if args.unit is not None:
        if args.block_mode is not None:
            print(f"  NOTE: both --unit and --block-mode given; --unit '{args.unit}' wins.",
                  flush=True)
        unit = args.unit
    elif args.block_mode is not None:
        unit = {"global": "chrom", "window": "bp"}[args.block_mode]
        print(f"  NOTE: --block-mode {args.block_mode} is DEPRECATED → --unit {unit}.",
              flush=True)
    else:
        unit = "chrom"
    args.unit = unit

    if args.emit_af_se and unit != "chrom":
        sys.exit("ERROR: --emit-af-se is implemented for --unit chrom only")
    if unit == "tsv" and not args.blocks_tsv:
        sys.exit("ERROR: --unit tsv requires --blocks-tsv <path>")

    # --local-only (default True) forces no anchor / no smoothing / no fallback in the
    # per-unit (non-chrom) path; _fit_unit_blocks does that forcing itself.

    print(f"=== {args.sample} (per-chrom unit={unit}"
          f"{', h-only' if args.h_only else ''}) ===", flush=True)
    print(f"  reads: {args.reads}", flush=True)
    t0 = time.time()

    global _VAR_CALLED
    if args.h_only:
        if unit != "chrom":
            sys.exit("ERROR: --h-only is supported in --unit chrom only")
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
        # Calibrated AF SE (Fisher + identifiability floor), filled per chrom when
        # --emit-af-se; replaces the legacy binomial `se` column at write time.
        se_af_global = np.full(n_records, np.nan, dtype=np.float32) if args.emit_af_se else None
        # Panel-level per-record called counts (h-independent QC metric). When
        # var_called is unavailable, assume the panel size F (all called).
        F_total = var_pa.shape[0]
        if _VAR_CALLED is not None:
            n_called_per_rec = np.asarray(_VAR_CALLED.sum(axis=0)).flatten().astype(np.int32)
        else:
            n_called_per_rec = np.full(n_records, F_total, dtype=np.int32)

    for chrom in args.chroms:
        if unit == "chrom":
            idx, freqs, info, h, _, extra = _fit_unit_chrom(
                chrom, args.kmer_pa_prefix, var_pa, var_meta,
                reads_input, args.threads, kmer_weight=args.kmer_weight,
                kmer_db=kmer_db, hash_size=args.hash_size, h_only=args.h_only,
                max_kmer_cov_mult=args.max_kmer_cov_mult,
                emit_af_se=args.emit_af_se, af_id_floor=args.af_id_floor,
                normalize=args.normalize, haploblock_eps=args.haploblock_eps)
            if h is None:
                continue
            h_save[chrom] = h
            if extra is not None:
                se_af_global[idx] = extra["af_se"]
                h_save[f"{chrom}_eff_rank"] = np.float32(extra["eff_rank"])
                h_save[f"{chrom}_cond"] = np.float32(extra["cond"])
                h_save[f"{chrom}_support_size"] = np.int32(extra["support_size"])
            if args.h_only:
                gc.collect()
                continue
        else:  # per-unit local fit: bp windows | tsv blocks | ld blocks
            # Resolve this chrom's block TSV. ld: compute-or-load (panel property);
            # tsv: the given --blocks-tsv; bp: None (fixed --window-bp windows).
            blocks_tsv = None
            if unit == "tsv":
                blocks_tsv = args.blocks_tsv
            elif unit == "ld":
                if args.ld_blocks:
                    blocks_tsv = args.ld_blocks
                else:
                    from .ld_partition import ld_blocks_tsv as _ld_tsv
                    vp_prefix = args.var_pa[:-len(".var_pa.npz")] if \
                        args.var_pa.endswith(".var_pa.npz") else args.var_pa
                    cache = f"{vp_prefix}.ld_blocks_r2_{args.ld_r2:.2f}_{chrom}.tsv"
                    print(f"  [{chrom}] --unit ld: r²={args.ld_r2} blocks "
                          f"({'cached' if os.path.exists(cache) else 'computing'}) → {cache}",
                          flush=True)
                    blocks_tsv = _ld_tsv(vp_prefix, chrom, args.ld_r2, cache,
                                         window=args.ld_window)
            idx, freqs, info, pack, _ = _fit_unit_blocks(
                chrom, args.kmer_pa_prefix, var_pa, var_meta,
                reads_input, args.threads,
                window_bp=args.window_bp,
                global_anchor_weight=args.global_anchor_weight,
                hmm_smooth_passes=args.hmm_smooth_passes,
                hmm_smooth_alpha=args.hmm_smooth_alpha,
                hmm_smooth_recomb_rate=args.hmm_smooth_recomb_rate,
                kmer_weight=args.kmer_weight,
                kmer_db=kmer_db, hash_size=args.hash_size,
                blocks_tsv=blocks_tsv,
                min_kmers_per_block=args.min_kmers_per_block,
                local_only=args.local_only,
                max_kmer_cov_mult=args.max_kmer_cov_mult,
                normalize=args.normalize,
                haploblock_eps=args.haploblock_eps,
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
        if args.emit_af_se:
            # Calibrated AF SE = sqrt(Fisher delta-method SE^2 + identifiability floor^2),
            # filled per chrom above. The principled per-record uncertainty.
            se_global = se_af_global
        else:
            # Legacy SE (Wald, using n_called as effective N). NaN if p is NaN or
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
        # ATOMIC write: build the full per-chrom TSV in a .tmp then os.replace() it
        # into place. Under preemptible (lowprio + --requeue) array runs, a task
        # killed mid-write must NOT leave a partial-but-non-empty ${SAMPLE}_${CHR}.tsv
        # — the runner's resume guard is a non-empty check, so a truncated file would
        # be skipped as "done" and silently concatenated into a chrom-truncated
        # genome-wide TSV. os.replace is atomic on the same filesystem.
        tmp_out = args.out + ".tmp"
        with open(tmp_out, "w") as f:
            f.write("chrom\tpos\tref_len\talt_len\talt_freq\tinfo\tn_called\tse\n")
            for i in range(n_records):
                af = freqs_global[i]
                af_str = f"{af:.5f}" if np.isfinite(af) else "NaN"
                inf_str = f"{info_global[i]:.5f}" if np.isfinite(info_global[i]) else "NaN"
                nc_str = f"{int(n_called_per_rec[i])}"
                se_str = f"{se_global[i]:.5f}" if np.isfinite(se_global[i]) else "NaN"
                f.write(f"{chrom_arr[i]}\t{pos_arr[i]}\t{ref_arr[i]}\t{alt_arr[i]}\t"
                        f"{af_str}\t{inf_str}\t{nc_str}\t{se_str}\n")
        os.replace(tmp_out, args.out)
        print(f"  wrote {n_records:,} records to {args.out}", flush=True)

    suffix = ".h_per_chrom.npz" if unit == "chrom" else ".h_blocks_per_chrom.npz"
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
