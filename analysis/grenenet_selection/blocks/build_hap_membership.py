#!/usr/bin/env python
"""STEP 1 of haploblock-frequency projection: build the persisted founder->haplotype
membership (the projection matrix M_b) on the dynld-K500 units.

For each dynld-K500 unit (analysis/grenenet_selection/blocks_mcf90/chr{N}_units_dynld_K500.tsv),
partition the 231 founders into haplotypes using the FAITHFUL HapFM-xmeans rule
(copied verbatim from block_cluster_pc1ve.cluster_founders, + seeded for reproducibility):
  <2 unique haplotypes  -> 1 cluster (the whole block segregates as one unit)
  2..6 unique           -> each unique haplotype is its own cluster (deterministic)
  >=7 unique            -> HapFM xmeans (k-means++ init, up to 30 clusters); seeded

Founder genotype matrix + unit slicing are built IDENTICALLY to dynamic_ld_blocks.py
(build_common_matrix MAF=0.05, MINCF=0.9, geno = raw>=0.5), so the [lo,hi) variant slice
of each unit reproduces the dynld map exactly. Founder ORDER == window-run h founder order
(var_pa meta 'founders'; verified equal to *.h_blocks_per_chrom.npz 'founders').

Output per chrom: analysis/grenenet_selection/blocks_mcf90/hap_membership/{chrlc}_hapmemb_K500.npz
  founders     (231,)        founder ids (projection-aligned to h_blocks)
  unit_chrom/start/end       (U,)  per-unit interval (== dynld unit map)
  unit_nvar/kmers/covered    (U,)  per-unit n_variants, panel_kmers, covered flag
  unit_d_uniq                (U,)  # unique founder haplotypes in the unit
  unit_k                     (U,)  # haplotype clusters kept
  unit_n_eff                 (U,)  1/sum(freq^2), founder-level effective # haplotypes
  labels        (U, 231) int16    founder -> cluster-id-within-unit  (the partition)
  hap_offset    (U+1,) int64       prefix-sum of unit_k -> global hap-id ranges
                                    (global hap id of founder f in unit u = hap_offset[u]+labels[u,f])

Also writes a flat haplotype registry CSV (one row per (unit,haplotype)) for inspection.

Run in the kmate env (needs pyclustering for xmeans).
  python build_hap_membership.py --chrom Chr1
"""
import argparse, os, random, sys, warnings
import numpy as np
import pandas as pd

# pyclustering 0.10.x calls numpy.warnings (removed in numpy>=1.24) -> shim before import
if not hasattr(np, "warnings"):
    np.warnings = warnings

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from recompute_blocks import build_common_matrix

MAF, MINCF = 0.05, 0.9          # IDENTICAL to dynamic_ld_blocks.py


def common_is_sv(var_pa_prefix, maf, min_called_frac):
    """is_sv aligned 1:1 to build_common_matrix's returned positions (SAME MAF/called filter
    and SAME first-kept-variant-per-unique-position selection), True where the chosen variant
    is an SV (|alt_len-ref_len| > 50). Used by --exclude-sv to remove SV columns from the
    founder clustering matrix so the clq90 haplotypes are NOT co-defined by the SVs later
    tested for tagging them (circularity control)."""
    import scipy.sparse as sp
    vp = sp.load_npz(f"{var_pa_prefix}.var_pa.npz").tocsc()
    vc = sp.load_npz(f"{var_pa_prefix}.var_called.npz").tocsc()
    meta = np.load(f"{var_pa_prefix}.meta.npz", allow_pickle=True)
    pos = meta["pos"].astype(np.int64)
    dl = np.abs(meta["alt_len"].astype(np.int64) - meta["ref_len"].astype(np.int64))
    F = vp.shape[0]
    na = np.asarray(vp.sum(0)).ravel(); nc = np.asarray(vc.sum(0)).ravel()
    af = np.divide(na, nc, out=np.full(len(na), np.nan), where=nc > 0)
    keep = (af > maf) & (af < 1 - maf) & (nc >= min_called_frac * F)
    seen = set(); is_sv = []
    for c in np.where(keep)[0]:
        p = int(pos[c])
        if p not in seen:
            seen.add(p); is_sv.append(bool(dl[c] > 50))
    return np.array(is_sv, bool)
BR = "analysis/grenenet_selection/blocks_mcf90"
HAPFM_KMAX = 7                  # HapFM: cluster (xmeans) only when #unique haplotypes >= this
SEED = 0


def xmeans_clustering(array):
    """HapFM utility_functions.xmeans_clustering (k-means++ init, X-Means up to 30
    clusters, BIC). array: list of binary haplotype vectors. REPRODUCIBLE: uses
    random_state=SEED and ccore=False (the pyclustering C++ backend ignores the seed,
    so ccore=True is nondeterministic; verified ccore=False is bit-stable run-to-run).
    Only deviation from HapFM's call is the optimizer backend (Python vs C++); same
    criterion + kmax."""
    from pyclustering.cluster.xmeans import xmeans
    from pyclustering.cluster.center_initializer import kmeans_plusplus_initializer
    random.seed(SEED); np.random.seed(SEED)
    initial_centers = kmeans_plusplus_initializer(array, 2, random_state=SEED).initialize()
    inst = xmeans(array, initial_centers, 30, random_state=SEED, ccore=False)
    inst.process()
    clusters_ = inst.get_clusters()
    labels = [0] * len(array)
    for i, grp in enumerate(clusters_):
        for j in grp:
            labels[j] = i
    return np.asarray(labels)


def cluster_founders(sub):
    """sub: F x m int 0/1 founder block matrix -> (per-founder cluster label, #unique haps).
    Dedup haplotypes, cluster the unique ones (HapFM rule), map founders back."""
    F = sub.shape[0]
    uniq, inv = np.unique(sub, axis=0, return_inverse=True)   # inv: founder -> uniq idx
    d = uniq.shape[0]
    if d < 2:
        return np.zeros(F, dtype=int), d
    if d < HAPFM_KMAX:
        uniq_labels = np.arange(d)                            # each unique hap its own cluster
    else:
        uniq_labels = xmeans_clustering([row.astype(float).tolist() for row in uniq])
    # compress labels to contiguous 0..k-1
    _, lab = np.unique(uniq_labels[inv], return_inverse=True)
    return lab, d


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--chrom", default="Chr1")
    ap.add_argument("--K", type=int, default=500)
    ap.add_argument("--units", default=None,
                    help="unit-map TSV (chrom,start_pos,end_pos,n_variants[,panel_kmers,covered]); "
                         "default = the dynld-K{K} map")
    ap.add_argument("--tag", default=None,
                    help="output suffix -> {chrlc}_hapmemb_{tag}.npz; default = K{K}")
    ap.add_argument("--exclude-sv", action="store_true",
                    help="drop SV columns (|alt_len-ref_len|>50) from the clustering matrix so "
                         "haplotypes are defined by SNP+indel only (SV self-tagging control)")
    a = ap.parse_args()
    chrlc = a.chrom.lower()
    units_path = a.units or f"{BR}/{chrlc}_units_dynld_K{a.K}.tsv"
    tag = a.tag or f"K{a.K}"

    print(f"[{a.chrom}] founder matrix (MAF={MAF}, MINCF={MINCF}) ...", flush=True)
    _, raw, positions = build_common_matrix(f"panel/arch3/{chrlc}/var_pa_231_arch3_{chrlc}", MAF, MINCF)
    geno = (raw >= 0.5).astype(np.int8)
    positions = np.asarray(positions)
    if a.exclude_sv:
        is_sv = common_is_sv(f"panel/arch3/{chrlc}/var_pa_231_arch3_{chrlc}", MAF, MINCF)
        assert len(is_sv) == geno.shape[1] == len(positions), (len(is_sv), geno.shape, len(positions))
        keepcol = ~is_sv
        geno = geno[:, keepcol]
        positions = positions[keepcol]
        print(f"[{a.chrom}] --exclude-sv: dropped {int(is_sv.sum()):,} SV cols -> "
              f"{int(keepcol.sum()):,} SNP/indel cols for clustering", flush=True)
    founders = np.load(f"panel/arch3/{chrlc}/var_pa_231_arch3_{chrlc}.meta.npz",
                       allow_pickle=True)["founders"].astype(str)
    F = geno.shape[0]
    assert F == len(founders) == 231, (F, len(founders))

    U = pd.read_csv(units_path, sep="\t")
    n = len(U)
    # k-mer coverage columns are dynld-only (irrelevant to the GWAS, which reads only
    # labels/offset/founders/start/end); default them when the unit map lacks them (e.g. clq0.9).
    kmers_col = U["panel_kmers"].to_numpy(np.int64) if "panel_kmers" in U else np.zeros(n, np.int64)
    covered_col = U["covered"].to_numpy(bool) if "covered" in U else np.ones(n, bool)
    print(f"[{a.chrom}] {n:,} units ({tag}) to partition from {units_path}", flush=True)

    labels = np.zeros((n, F), dtype=np.int16)
    d_uniq = np.zeros(n, dtype=np.int32)
    k_clust = np.zeros(n, dtype=np.int32)
    n_eff = np.zeros(n, dtype=np.float64)
    reg = []   # flat (unit, haplotype) registry

    for ui in range(n):
        s, e = int(U.start_pos[ui]), int(U.end_pos[ui])
        lo = int(np.searchsorted(positions, s))
        hi = int(np.searchsorted(positions, e, side="right"))
        sub = geno[:, lo:hi]
        if sub.shape[1] == 0:                      # no kept variant in interval -> single block unit
            lab = np.zeros(F, dtype=int); d = 1
        else:
            lab, d = cluster_founders(sub)
        labels[ui] = lab
        d_uniq[ui] = d
        k = int(lab.max()) + 1
        k_clust[ui] = k
        freq = np.bincount(lab, minlength=k) / F
        n_eff[ui] = 1.0 / np.sum(freq ** 2)
        for c in range(k):
            reg.append(dict(chrom=a.chrom, unit_idx=ui, unit_start=s, unit_end=e,
                            cluster=c, n_founders=int((lab == c).sum()),
                            cluster_freq=round(float(freq[c]), 4),
                            n_variants=int(sub.shape[1]),
                            panel_kmers=int(kmers_col[ui]), covered=bool(covered_col[ui])))
        if (ui + 1) % 2000 == 0:
            print(f"  ...{ui+1}/{n} units", flush=True)

    hap_offset = np.concatenate([[0], np.cumsum(k_clust)]).astype(np.int64)
    total_haps = int(hap_offset[-1])

    out = f"{BR}/hap_membership/{chrlc}_hapmemb_{tag}.npz"
    np.savez_compressed(
        out,
        founders=founders,
        unit_chrom=U.chrom.values.astype(str),
        unit_start=U.start_pos.values.astype(np.int64),
        unit_end=U.end_pos.values.astype(np.int64),
        unit_nvar=U.n_variants.values.astype(np.int64),
        unit_kmers=kmers_col,
        unit_covered=covered_col,
        unit_d_uniq=d_uniq, unit_k=k_clust, unit_n_eff=n_eff,
        labels=labels, hap_offset=hap_offset,
    )
    reg_df = pd.DataFrame(reg)
    reg_csv = f"{BR}/hap_membership/{chrlc}_hapmemb_{tag}_registry.csv"
    reg_df.to_csv(reg_csv, index=False)

    nx = int((d_uniq >= HAPFM_KMAX).sum())
    print(f"\n=== {a.chrom} haplotype membership ({tag}) ===")
    print(f"units: {n:,}   total haplotypes: {total_haps:,}   "
          f"(median {int(np.median(k_clust))} hap/unit, max {int(k_clust.max())})")
    print(f"unique-haplotype distribution: median d_uniq {int(np.median(d_uniq))}, "
          f"{nx:,} units ({100*nx/n:.1f}%) used xmeans (d_uniq>={HAPFM_KMAX})")
    print(f"founder-level n_eff: median {np.median(n_eff):.2f}, "
          f"%units n_eff>2 {100*(n_eff>2).mean():.0f}%")
    print(f"covered units: {int(covered_col.sum()):,} ({100*covered_col.mean():.0f}%)")
    print(f"-> {out}\n-> {reg_csv}")


if __name__ == "__main__":
    main()
