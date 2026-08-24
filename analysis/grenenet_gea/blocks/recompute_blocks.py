#!/usr/bin/env python
"""Recompute LD haploblocks on OUR panel (merged_231, all classes) using HapFM's
exact partition method: coarse CompleteLDPartition (custom windowed-r2) + fine
BigLD (gpart BigLD.R, CLQmode="density").

This is faithful reuse of HapFM/bin/block_partition.py + BigLD.R — the four
partition functions below are copied verbatim from HapFM (find_ld,
CompleteLDPartition, BigLD_partition, keep_largest_overlaps). The only new code is
the input side: we build the 231-founder genotype matrix from var_pa (all-class:
SNP+indel+SV) instead of a VCF, impute founder missingness to the per-variant
mean, MAF-filter, and define blocks on UNIQUE positions (our biallelic
decomposition has co-located records sharing a position).

Output: <out>_blocks_clq{CLQ}.tsv  cols: chrom start_pos end_pos n_variants
Run inside the `kmate` env (numpy/scipy); needs `--rscript` pointing at the
`bigld` env's Rscript and `--bigld-dir` at HapFM/bin (where BigLD.R lives).
"""
import argparse, os, subprocess, sys
import numpy as np
import scipy.sparse as sp

# ----------------------------------------------------------------------------
# HapFM partition functions — COPIED VERBATIM from HapFM/bin/block_partition.py
# and utility_functions.keep_largest_overlaps (do not "improve"; faithful reuse).
# ----------------------------------------------------------------------------
def find_ld(i, snps, cutoff, window_size):
    n_inds, n_snps = snps.shape
    left = max(i - window_size, 0)
    right = min(i + window_size, n_snps)
    left_snps_window = snps[:, left:(i + 1)]
    right_snps_window = snps[:, i:(right + 1)]
    left_cor = np.matmul(np.transpose(left_snps_window), snps[:, i]) / n_inds
    left_cor_rev = np.flip(left_cor)
    right_cor = np.matmul(np.transpose(right_snps_window), snps[:, i]) / n_inds
    left_list_ = np.where(left_cor_rev ** 2 > cutoff)[0]
    for j in range(len(left_list_) - 1):
        if left_list_[j + 1] - left_list_[j] > 10:
            left_list_ = left_list_[:j + 1]
            break
    left_list_ = np.flip(left_list_) * -1
    right_list_ = np.where(right_cor ** 2 > cutoff)[0]
    for j in range(len(right_list_) - 1):
        if right_list_[j + 1] - right_list_[j] > 10:
            right_list_ = right_list_[:j + 1]
            break
    SNPinLD_index = np.unique(np.concatenate((left_list_, right_list_))) + i
    return SNPinLD_index


def CompleteLDPartition(standardized_genotype_matrix, cutoff, window_size):
    n_inds, n_snps = standardized_genotype_matrix.shape
    snp_list = {}
    cummax_list = []
    max_list = []
    boundary = []
    alone_SNPs_index = []
    for i in range(n_snps):
        snp_list[i] = find_ld(i, snps=standardized_genotype_matrix, cutoff=0.1, window_size=50)
        if len(snp_list[i]) == 1:
            alone_SNPs_index.append(i)
    print("QUALITY CHECK: %d snps not in LD (r2<0.1) with 50 neighbours." % len(alone_SNPs_index))
    print("window_size %d, correlation cutoff %f" % (window_size, cutoff))
    for i in range(n_snps):
        snp_list[i] = find_ld(i, snps=standardized_genotype_matrix, cutoff=cutoff, window_size=window_size)
    for i in range(len(snp_list)):
        max_list.append(np.max(snp_list[i]) if len(snp_list[i]) > 0 else i)
    cummax_list.append(max_list[0])
    for i in range(1, len(max_list)):
        cummax_list.append(max(max_list[i], cummax_list[i - 1]))
    idx = np.where(cummax_list - np.array(range(n_snps)) == 0)[0]
    boundary_ = np.concatenate(([-1], np.array(idx)))
    for i in range(len(boundary_) - 1):
        left = boundary_[i] + 1
        right = boundary_[i + 1]
        boundary.append([left, right] if right - left > 0 else [left])
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
    return boundary, alone_SNPs_index


def keep_largest_overlaps(blocks):
    ivals = [(int(s), int(e)) for s, e in blocks]
    n = len(ivals)
    if n <= 1:
        return blocks
    adj = [[] for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            a0, a1 = ivals[i]
            b0, b1 = ivals[j]
            if not (a1 < b0 or b1 < a0):
                adj[i].append(j)
                adj[j].append(i)
    visited = [False] * n
    result = []
    for i in range(n):
        if visited[i]:
            continue
        stack = [i]
        visited[i] = True
        cluster = [i]
        while stack:
            u = stack.pop()
            for v in adj[u]:
                if not visited[v]:
                    visited[v] = True
                    stack.append(v)
                    cluster.append(v)
        best = max(cluster, key=lambda k: (ivals[k][1] - ivals[k][0], -ivals[k][0], ivals[k][1]))
        result.append(list(ivals[best]))
    return sorted(result, key=lambda x: x[0])


def BigLD_partition(DIR, IndepLD_breakpoints_index, geno_matrix, variant_names,
                    variant_positions, CLQcut, prefix, ch, maf, rscript):
    fine_breakpoints_ch = []
    for I in range(len(IndepLD_breakpoints_index)):
        left = IndepLD_breakpoints_index[I][0]
        right = IndepLD_breakpoints_index[I][1]
        tmp_names = variant_names[left:right + 1]
        tmp_positions = variant_positions[left:right + 1]
        tmp_matrix = geno_matrix[:, left:right + 1]
        import pandas as pd
        pd.DataFrame(tmp_matrix, columns=tmp_names).to_csv(
            f"{prefix}_{I}_geno_matrix.btmp", sep="\t", header=True, index=False)
        with open(f"{prefix}_{I}_snpINFO.btmp", "w") as INFO:
            INFO.write("chrN\trsID\tbp\n")
            for j in range(len(tmp_positions)):
                INFO.write(f"tmp\t{tmp_names[j]}\t{tmp_positions[j]}\n")
        cmd = (f"{rscript} {DIR}/BigLD.R -g {prefix}_{I}_geno_matrix.btmp "
               f"-s {prefix}_{I}_snpINFO.btmp -c {CLQcut} -m {maf} -o {prefix}_{I}")
        try:
            blocks = []
            subprocess.check_call(cmd, shell=True)
            with open(f"{prefix}_{I}_res_btmp.txt") as INPUT:
                INPUT.readline()
                for line in INPUT:
                    items = line.split("\t")
                    blocks.append([variant_positions.index(int(items[5])),
                                   variant_positions.index(int(items[6]))])
            blocks[-1][1] = right
            blocks = keep_largest_overlaps(blocks)
            fine_breakpoints_ch.extend(blocks)
        except subprocess.CalledProcessError:
            print("BigLD could not partition region %i-%i; keeping coarse." % (left, right))
            fine_breakpoints_ch.append(IndepLD_breakpoints_index[I])
        subprocess.call(f"rm -f {prefix}_{I}_*btmp*", shell=True)
    return fine_breakpoints_ch


# ----------------------------------------------------------------------------
# Input side (new): build the founder genotype matrix from var_pa, all classes.
# ----------------------------------------------------------------------------
def build_common_matrix(var_pa_prefix, maf, min_called_frac):
    """Return (standardized matrix [F x Ncommon], raw 0/1 matrix, positions list)
    on the MAF-filtered, unique-position common set (one variant per position)."""
    vp = sp.load_npz(f"{var_pa_prefix}.var_pa.npz").toarray().astype(np.float64)   # F x R, 0/1
    vc = sp.load_npz(f"{var_pa_prefix}.var_called.npz").toarray().astype(np.float64)
    meta = np.load(f"{var_pa_prefix}.meta.npz", allow_pickle=True)
    pos = meta["pos"].astype(np.int64)
    F = vp.shape[0]
    n_called = vc.sum(0)
    n_alt = vp.sum(0)
    af = np.divide(n_alt, n_called, out=np.full_like(n_alt, np.nan), where=n_called > 0)
    keep = (af > maf) & (af < 1 - maf) & (n_called >= min_called_frac * F)
    # unique positions: keep the FIRST kept variant at each position (blocks are intervals)
    order = np.where(keep)[0]
    seen = set(); uniq = []
    for c in order:
        p = int(pos[c])
        if p not in seen:
            seen.add(p); uniq.append(c)
    uniq = np.array(uniq)
    print(f"  variants kept after MAF/called filter: {keep.sum()}; unique positions: {len(uniq)}")
    raw = vp[:, uniq].copy()           # F x Nuniq, 0/1
    miss = vc[:, uniq] == 0
    afc = af[uniq]
    # impute missing -> per-variant mean (the standard LD missing-handling)
    raw[miss] = np.repeat(afc[None, :], F, axis=0)[miss]
    # standardize per column (mean 0, sd 1) so matmul/n = Pearson r (HapFM uses preprocessing.scale)
    mu = raw.mean(0); sd = raw.std(0); sd[sd == 0] = 1.0
    std = (raw - mu) / sd
    return std, raw, [int(p) for p in pos[uniq]]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--var-pa-prefix", required=True,
                    help="e.g. panel/arch3/chr1/var_pa_231_arch3_chr1")
    ap.add_argument("--chrom", required=True)
    ap.add_argument("--clqcut", type=float, default=0.5, help="BigLD fine r2 cutoff")
    ap.add_argument("--corr", type=float, default=0.5, help="coarse CompleteLDPartition r2")
    ap.add_argument("--window", type=int, default=50)
    ap.add_argument("--maf", type=float, default=0.05)
    ap.add_argument("--min-called-frac", type=float, default=0.5)
    ap.add_argument("--out", required=True, help="output prefix")
    ap.add_argument("--no-bigld", action="store_true",
                    help="skip BigLD fine-split; output the coarse CompleteLDPartition "
                         "blocks directly (pure numpy, no R/gpart). Block tightness is "
                         "controlled by --corr.")
    ap.add_argument("--rscript", default=None, help="path to bigld env Rscript (BigLD mode)")
    ap.add_argument("--bigld-dir", default=None, help="dir with BigLD.R (BigLD mode)")
    a = ap.parse_args()

    print(f"[{a.chrom}] building common matrix (maf>{a.maf}) ...")
    std, raw, positions = build_common_matrix(a.var_pa_prefix, a.maf, a.min_called_frac)
    print(f"[{a.chrom}] coarse CompleteLDPartition (corr={a.corr}) ...")
    coarse, _ = CompleteLDPartition(std, cutoff=a.corr, window_size=a.window)
    print(f"[{a.chrom}] {len(coarse)} coarse independent-LD regions")
    if a.no_bigld:
        fine = [b if len(b) == 2 else [b[0], b[0]] for b in coarse]
        tag = f"corr{a.corr}"
    else:
        names = [str(p) for p in positions]
        print(f"[{a.chrom}] fine BigLD (CLQcut={a.clqcut}) ...")
        fine = BigLD_partition(a.bigld_dir, coarse, raw, names, positions,
                               a.clqcut, a.out, a.chrom, a.maf, a.rscript)
        tag = f"clq{a.clqcut}"
    outf = f"{a.out}_blocks_{tag}.tsv"
    with open(outf, "w") as o:
        o.write("chrom\tstart_pos\tend_pos\tn_variants\n")
        for b in fine:
            l, r = b[0], b[1]
            o.write(f"{a.chrom}\t{positions[l]}\t{positions[r]}\t{r - l + 1}\n")
    print(f"[{a.chrom}] wrote {len(fine)} blocks -> {outf}")


if __name__ == "__main__":
    main()
