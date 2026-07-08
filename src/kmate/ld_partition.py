"""r²-LD block partitions from kMate's OWN founder panel (var_pa) — no external
SNP panel. Markers are the panel's var_pa records (founder × variant), filtered to
a common, well-called set; blocks are the hapFIRE CompleteLDPartition algorithm
(r² cutoff + fixed window), tiled to genomic midpoints so every downstream k-mer
lands in exactly one block.

Used by the `--unit ld` estimator mode (per_sample_per_chrom): the LD blocks are a
PANEL property (sample-independent), so they are computed once per (panel, chrom, r²)
and cached to a TSV that the block-EM path consumes via `--blocks-tsv`.

`find_ld` / `complete_ld_partition` are copied faithfully from hapFIRE
(external/HapFIRE/haplotype_generation.py); they are pure-numpy (no cvxpy).
"""
from __future__ import annotations
import os
import numpy as np
import scipy.sparse as sp

# TAIR10 chromosome lengths (bp) — so the final block tiles to the true chrom end.
TAIR10_CHR_LEN = {
    "Chr1": 30427671, "Chr2": 19698289, "Chr3": 23459830,
    "Chr4": 18585056, "Chr5": 26975502,
}


def find_ld(i, snps, cutoff, window_size):
    n_inds, n_snps = snps.shape
    left = max(i - window_size, 0)
    right = min(i + window_size, n_snps)
    left_cor = np.matmul(np.transpose(snps[:, left:(i + 1)]), snps[:, i]) / n_inds
    left_cor_rev = np.flip(left_cor)
    right_cor = np.matmul(np.transpose(snps[:, i:(right + 1)]), snps[:, i]) / n_inds
    left_list_ = np.where(left_cor_rev ** 2 > cutoff)[0]
    for j in range(len(left_list_) - 1):
        if left_list_[j + 1] - left_list_[j] > 25:
            left_list_ = left_list_[:j + 1]
            break
    left_list_ = np.flip(left_list_) * -1
    right_list_ = np.where(right_cor ** 2 > cutoff)[0]
    for j in range(len(right_list_) - 1):
        if right_list_[j + 1] - right_list_[j] > 25:
            right_list_ = right_list_[:j + 1]
            break
    return np.unique(np.concatenate((left_list_, right_list_))) + i


def complete_ld_partition(standardized_genotype_matrix, cutoff, window_size):
    """Return (boundary, alone_SNPs) — hapFIRE's independent-LD-block partition.
    `boundary` is a list of [left_idx, right_idx] (or [idx]) in marker-index space."""
    n_inds, n_snps = standardized_genotype_matrix.shape
    snp_list = {}
    max_list = []
    boundary = []
    alone = []
    for i in range(n_snps):
        snp_list[i] = find_ld(i, standardized_genotype_matrix, 0.1, 50)
        if len(snp_list[i]) == 1:
            alone.append(i)
    for i in range(n_snps):
        snp_list[i] = find_ld(i, standardized_genotype_matrix, cutoff, window_size)
    for i in range(len(snp_list)):
        max_list.append(int(np.max(snp_list[i])) if len(snp_list[i]) > 0 else i)
    cummax = [max_list[0]]
    for i in range(1, len(max_list)):
        cummax.append(max(max_list[i], cummax[i - 1]))
    idx = np.where(np.array(cummax) - np.arange(n_snps) == 0)[0]
    b_ = np.concatenate(([-1], idx))
    for i in range(len(b_) - 1):
        l = b_[i] + 1
        r = b_[i + 1]
        boundary.append([l, r] if r - l > 0 else [l])
    j = 0
    while j < len(boundary):
        if len(boundary[j]) == 1:
            if j == 0:
                boundary[j + 1][0] = boundary[j][0]
                del boundary[j]
            else:
                boundary[j - 1][1] = boundary[j][0]
                del boundary[j]
        else:
            j += 1
    return boundary, alone


def _standardize(G):
    """Column-standardize (mean 0, unit variance); constant columns → 0.
    Matches sklearn.preprocessing.scale (ddof=0) without the sklearn dependency."""
    mu = G.mean(axis=0)
    sd = G.std(axis=0)
    sd_safe = np.where(sd == 0, 1.0, sd)
    return (G - mu) / sd_safe


def compute_ld_blocks(var_pa_prefix, chrom, r2, window=100, maf=0.05,
                      callrate_min=0.9, chrom_len=None):
    """Compute r²-LD blocks for one chromosome from a var_pa panel.

    var_pa_prefix: path prefix such that `{prefix}.var_pa.npz`, `.var_called.npz`,
        `.meta.npz` exist (founder × variant panel + call mask + pos/chrom meta).
    Returns a list of (chrom, start_pos, end_pos, n_variants) rows tiling the whole
    chromosome (contiguous, so every k-mer position lands in exactly one block).
    Missingness-aware: markers filtered to callrate≥callrate_min and called-MAF in
    (maf, 1-maf); uncalled genotypes mean-imputed to the marker's called freq.
    """
    vp = sp.load_npz(var_pa_prefix + ".var_pa.npz").tocsr()
    vc = sp.load_npz(var_pa_prefix + ".var_called.npz").tocsr()
    meta = np.load(var_pa_prefix + ".meta.npz", allow_pickle=True)
    F = vp.shape[0]
    pos = np.asarray(meta["pos"]).astype(np.int64)
    rec_chrom = np.asarray(meta["chrom"]).astype(str) if "chrom" in meta.files else None

    on_chrom = (rec_chrom == str(chrom)) if rec_chrom is not None else np.ones(len(pos), bool)
    called = np.asarray(vc.sum(axis=0)).ravel().astype(np.float64)
    alt = np.asarray(vp.sum(axis=0)).ravel().astype(np.float64)
    callrate = called / F
    with np.errstate(invalid="ignore", divide="ignore"):
        freqc = np.where(called > 0, alt / called, 0.0)
    common = np.flatnonzero(on_chrom & (callrate >= callrate_min)
                            & (freqc > maf) & (freqc < 1 - maf))
    common = common[np.argsort(pos[common], kind="stable")]     # position-sorted
    if len(common) == 0:
        raise ValueError(f"compute_ld_blocks: no common markers on {chrom}")
    cpos = pos[common]

    G = np.asarray(vp[:, common].todense(), dtype=np.float64)   # F × n_common
    Cs = np.asarray(vc[:, common].todense(), dtype=bool)
    fc = freqc[common]
    rows, cols = np.where(~Cs)                                  # uncalled → marker called-freq
    G[rows, cols] = fc[cols]
    del Cs
    Gs = _standardize(G)
    del G

    boundary, _ = complete_ld_partition(Gs, cutoff=r2, window_size=window)
    edges = sorted((int(b[0]), int(b[-1] if len(b) > 1 else b[0])) for b in boundary)

    if chrom_len is None:
        chrom_len = TAIR10_CHR_LEN.get(str(chrom))
        if chrom_len is None:
            chrom_len = int(cpos.max())                          # fallback: last marker
    # tile block boundaries to genomic midpoints between adjacent LD blocks
    starts = [1]
    for k in range(len(edges) - 1):
        mid = (int(cpos[edges[k][1]]) + int(cpos[edges[k + 1][0]])) // 2
        starts.append(mid + 1)
    out = []
    for k, (l, r) in enumerate(edges):
        s = starts[k]
        e = (starts[k + 1] - 1) if k + 1 < len(starts) else int(chrom_len)
        out.append((str(chrom), s, e, r - l + 1))
    return out


def ld_blocks_tsv(var_pa_prefix, chrom, r2, cache_path, window=100, maf=0.05,
                  callrate_min=0.9, chrom_len=None, recompute=False):
    """Return `cache_path`, computing + writing the r²-LD blocks TSV if absent.

    The TSV (chrom, start_pos, end_pos, n_variants) is what the block-EM path loads
    via `--blocks-tsv`. LD blocks are a panel property, so the cache is keyed by
    (panel, chrom, r²) and reused across samples.
    """
    if os.path.exists(cache_path) and not recompute:
        return cache_path
    rows = compute_ld_blocks(var_pa_prefix, chrom, r2, window=window, maf=maf,
                             callrate_min=callrate_min, chrom_len=chrom_len)
    tmp = cache_path + ".tmp"
    with open(tmp, "w") as fh:
        fh.write("chrom\tstart_pos\tend_pos\tn_variants\n")
        for c, s, e, n in rows:
            fh.write(f"{c}\t{s}\t{e}\t{n}\n")
    os.replace(tmp, cache_path)
    return cache_path
