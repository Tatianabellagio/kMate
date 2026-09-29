"""Genome-wide verification: r2=0.1 LD blocks + eps=0 haploblock K_b on Chr1-5.

For each chromosome:
  1. define r2=0.1 LD blocks from arch3 var_pa (SAME filtering + CompleteLDPartition
     as gen_ld_partitions.py -- copied verbatim), tiled to genomic midpoints;
  2. compute K_b per block (exact k-mer identity) on the arch3_filt2inv k-mer panel,
     with the CORRECT per-k-mer -> block assignment via bubble_id.
Writes ld_blocks_r2_0.10_<Chr>.tsv per chrom + a combined genome TSV, and a
K_b summary per chromosome. Confirms whether all chromosomes land at K_b~231.
"""
import time, json
import numpy as np, scipy.sparse as sp
from sklearn import preprocessing

ROOT = "/global/scratch/users/tbellg/kmate"
OUT = f"{ROOT}/analysis/grenenet_selection/blocks/hap_blocks"
KMDIR = f"{ROOT}/data/kmer_pa_231_arch3_filt2inv"
CUTOFF = 0.1; WINDOW = 100; MAF = 0.05; CALLRATE_MIN = 0.9; F = 231
KMER_FLOOR = 10_000
CHR_LEN = {"Chr1": 30427671, "Chr2": 19698289, "Chr3": 23459830,
           "Chr4": 18585056, "Chr5": 26975502}


# --- CompleteLDPartition copied verbatim from gen_ld_partitions.py / hapFIRE ---
def find_ld(i, snps, cutoff, window_size):
    n_inds, n_snps = snps.shape
    left = max(i - window_size, 0); right = min(i + window_size, n_snps)
    left_cor = np.matmul(np.transpose(snps[:, left:(i + 1)]), snps[:, i]) / n_inds
    left_cor_rev = np.flip(left_cor)
    right_cor = np.matmul(np.transpose(snps[:, i:(right + 1)]), snps[:, i]) / n_inds
    left_list_ = np.where(left_cor_rev ** 2 > cutoff)[0]
    for j in range(len(left_list_) - 1):
        if left_list_[j + 1] - left_list_[j] > 25:
            left_list_ = left_list_[:j + 1]; break
    left_list_ = np.flip(left_list_) * -1
    right_list_ = np.where(right_cor ** 2 > cutoff)[0]
    for j in range(len(right_list_) - 1):
        if right_list_[j + 1] - right_list_[j] > 25:
            right_list_ = right_list_[:j + 1]; break
    return np.unique(np.concatenate((left_list_, right_list_))) + i


def CompleteLDPartition(G, cutoff, window_size):
    n_inds, n_snps = G.shape
    snp_list = {}; max_list = []; boundary = []; alone = []
    for i in range(n_snps):
        snp_list[i] = find_ld(i, G, 0.1, 50)
        if len(snp_list[i]) == 1:
            alone.append(i)
    for i in range(n_snps):
        snp_list[i] = find_ld(i, G, cutoff, window_size)
    for i in range(len(snp_list)):
        max_list.append(np.max(snp_list[i]) if len(snp_list[i]) > 0 else i)
    cummax = [max_list[0]]
    for i in range(1, len(max_list)):
        cummax.append(max(max_list[i], cummax[i - 1]))
    idx = np.where(cummax - np.array(range(n_snps)) == 0)[0]
    b_ = np.concatenate(([-1], np.array(idx)))
    for i in range(len(b_) - 1):
        l = b_[i] + 1; r = b_[i + 1]
        boundary.append([l, r] if r - l > 0 else [l])
    j = 0
    while j < len(boundary):
        if len(boundary[j]) == 1:
            if j == 0:
                boundary[j + 1][0] = boundary[j][0]; del boundary[j]
            else:
                boundary[j - 1][1] = boundary[j][0]; del boundary[j]
        else:
            j += 1
    return boundary, alone


def blocks_for_chrom(chrom):
    cl = chrom.lower()
    VP = f"{ROOT}/panel/arch3/{cl}/var_pa_231_arch3_{cl}"
    vp = sp.load_npz(VP + ".var_pa.npz").tocsr()
    vc = sp.load_npz(VP + ".var_called.npz").tocsr()
    pos = np.load(VP + ".meta.npz", allow_pickle=True)["pos"].astype(np.int64)
    called = np.asarray(vc.sum(0)).ravel().astype(float)
    alt = np.asarray(vp.sum(0)).ravel().astype(float)
    with np.errstate(invalid="ignore", divide="ignore"):
        freqc = np.where(called > 0, alt / called, 0.0)
    common = np.flatnonzero((called / F >= CALLRATE_MIN) & (freqc > MAF) & (freqc < 1 - MAF))
    common = common[np.argsort(pos[common], kind="stable")]
    cpos = pos[common]
    G = np.asarray(vp[:, common].todense(), dtype=np.float64)
    Cs = np.asarray(vc[:, common].todense(), dtype=bool)
    fc = freqc[common]
    miss = ~Cs; rr, cc = np.where(miss); G[rr, cc] = fc[cc]
    Gs = preprocessing.scale(G)
    boundary, alone = CompleteLDPartition(Gs, CUTOFF, WINDOW)
    edges = sorted((int(b[0]), int(b[-1])) for b in boundary)
    starts = [1]
    for k in range(len(edges) - 1):
        mid = (int(cpos[edges[k][1]]) + int(cpos[edges[k + 1][0]])) // 2
        starts.append(mid + 1)
    rows = []
    for k, (l, r) in enumerate(edges):
        s = starts[k]; e = (starts[k + 1] - 1) if k + 1 < len(starts) else CHR_LEN[chrom]
        rows.append((chrom, s, e, r - l + 1))
    return rows, len(common)


def kb_for_blocks(chrom, blocks):
    pre = f"{KMDIR}/kmer_pa_{chrom}"
    m = np.load(pre + ".meta.npz", allow_pickle=True)
    bub = m["bubble_id"].astype(np.int64)
    bs = m["bubble_start"].astype(np.int64); be = m["bubble_end"].astype(np.int64)
    kpos = (bs[bub] + be[bub]) // 2
    K = sp.load_npz(pre + ".kmer_pa.npz").astype(bool).tocsc()
    order = np.argsort(kpos); spos = kpos[order]
    Kb = []; nk = []
    for (_, s, e, _) in blocks:
        lo = np.searchsorted(spos, s, "left"); hi = np.searchsorted(spos, e, "right")
        cols = order[lo:hi]
        if len(cols) == 0:
            continue
        sub = K[:, cols].tocsr()
        sig = [hash(sub.indices[sub.indptr[i]:sub.indptr[i + 1]].tobytes()) for i in range(F)]
        Kb.append(len(np.unique(sig))); nk.append(len(cols))
    return np.array(Kb), np.array(nk)


def main():
    t0 = time.time()
    combined = []
    summary = {}
    print(f"{'chrom':6s} {'blocks':>6s} {'Kb_min':>6s} {'Kb_med':>6s} {'Kb_max':>6s} "
          f"{'nk_med':>9s} {'Kb<200&krich':>12s} {'Kb<100':>6s}")
    for chrom in ["Chr1", "Chr2", "Chr3", "Chr4", "Chr5"]:
        blocks, ncommon = blocks_for_chrom(chrom)
        with open(f"{OUT}/ld_blocks_r2_0.10_{chrom}.tsv", "w") as fh:
            fh.write("chrom\tstart_pos\tend_pos\tn_variants\n")
            for c, s, e, n in blocks:
                fh.write(f"{c}\t{s}\t{e}\t{n}\n")
        combined.extend(blocks)
        Kb, nk = kb_for_blocks(chrom, blocks)
        real = int(((Kb < 200) & (nk > KMER_FLOOR)).sum())
        deep = int((Kb < 100).sum())
        summary[chrom] = dict(n_blocks=len(blocks), n_common_markers=ncommon,
                              Kb_min=int(Kb.min()), Kb_median=int(np.median(Kb)), Kb_max=int(Kb.max()),
                              nk_median=int(np.median(nk)),
                              n_Kb_lt200_krich=real, n_Kb_lt100=deep)
        print(f"{chrom:6s} {len(blocks):6d} {Kb.min():6d} {int(np.median(Kb)):6d} {Kb.max():6d} "
              f"{int(np.median(nk)):9,d} {real:12d} {deep:6d}  [+{time.time()-t0:.0f}s]", flush=True)
    with open(f"{OUT}/ld_blocks_r2_0.10_genome.tsv", "w") as fh:
        fh.write("chrom\tstart_pos\tend_pos\tn_variants\n")
        for c, s, e, n in combined:
            fh.write(f"{c}\t{s}\t{e}\t{n}\n")
    json.dump(summary, open(f"{OUT}/kb_genomewide_summary.json", "w"), indent=2)
    print(f"\nwrote per-chrom + genome block TSVs + kb_genomewide_summary.json ({time.time()-t0:.0f}s)")
    print("CONCLUSION: if Kb_median ~231 and n_Kb_lt200_krich ~0 on every chrom, all "
          "chromosomes agree -> haploblocks reduce to 231 at r2=0.1 genome-wide.")


if __name__ == "__main__":
    main()
