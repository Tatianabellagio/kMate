#!/usr/bin/env python
"""PROTOTYPE (Chr1): does HapFM haplotype-clustering turn low-VE multi-haplotype
blocks into high-VE single-haplotype units?

For each multi-haplotype Chr1 block (gated on EFFECTIVE diversity n_eff>GATE):
  1. cluster the 231 founders' haplotypes with HapFM's default xmeans
     (faithful reuse of utility_functions.xmeans_clustering; <7 unique haps ->
     each unique haplotype is its own cluster, exactly as HapFM does).
  2. for each cluster, take its SIGNATURE variants (alt in >=50% of the cluster's
     founders AND >=0.5 higher than outside) = the SNPs/indels/SVs that travel
     with that haplotype.
  3. compute PC1 variance-explained of those signature variants' AF across the
     gen9 POOLS (out-of-sample: pools never touched clustering).

Compare per-cluster PC1-VE to the whole-block PC1-VE. If clustering is real, a
block that was low-VE (mixed) should decompose into clusters each with HIGH VE.

Run in kmate env (needs pyclustering for xmeans). Usage:
  block_cluster_pc1ve.py [--gate 2.0] [--chrom Chr1] [--max-blocks N]
"""
import argparse, os, sys, warnings
import numpy as np
import pandas as pd
# pyclustering 0.10.x calls numpy.warnings, removed in numpy>=1.24 -> shim it
if not hasattr(np, "warnings"):
    np.warnings = warnings
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from recompute_blocks import build_common_matrix
from eval_block_coherence import load_allclass_af, pc1_var_explained

BR = "analysis/grenenet_selection/blocks/results/blocks_recompute"
MAF, MINCF = 0.05, 0.5
HAPFM_KMAX = 7          # HapFM clusters only when #unique haplotypes >= this
SIG_PRESENT = 0.5       # variant alt-freq in cluster >= this
SIG_DIFF = 0.5          # and (in - out) >= this  -> a haplotype-defining variant


def xmeans_clustering(array):
    """FAITHFUL copy of HapFM utility_functions.xmeans_clustering (k-means++ init,
    X-Means up to 30 clusters). array: list of binary haplotype vectors."""
    from pyclustering.cluster.xmeans import xmeans
    from pyclustering.cluster.center_initializer import kmeans_plusplus_initializer
    initial_centers = kmeans_plusplus_initializer(array, 2).initialize()
    inst = xmeans(array, initial_centers, 30)
    inst.process()
    clusters_ = inst.get_clusters()
    labels = [0] * len(array)
    for i, grp in enumerate(clusters_):
        for j in grp:
            labels[j] = i
    return np.asarray(labels)


def cluster_founders(sub):
    """sub: F x m int 0/1 founder block matrix -> per-founder cluster label.
    Dedup haplotypes, cluster the unique ones (HapFM rule), map founders back."""
    F = sub.shape[0]
    uniq, inv = np.unique(sub, axis=0, return_inverse=True)   # inv: founder->uniq idx
    d = uniq.shape[0]
    if d < 2:
        return np.zeros(F, dtype=int), d
    if d < HAPFM_KMAX:
        uniq_labels = np.arange(d)                            # each unique hap its own cluster
    else:
        uniq_labels = xmeans_clustering([row.astype(float).tolist() for row in uniq])
    return uniq_labels[inv], d


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--chrom", default="Chr1")
    ap.add_argument("--gate", type=float, default=2.0, help="cluster blocks with n_eff > gate")
    ap.add_argument("--max-blocks", type=int, default=0, help="0 = all gated blocks")
    a = ap.parse_args()
    chrlc = a.chrom.lower()

    print(f"[{a.chrom}] founder matrix + gen9 pools ...", flush=True)
    _, raw, positions = build_common_matrix(f"panel/arch3/{chrlc}/var_pa_231_arch3_{chrlc}", MAF, MINCF)
    positions = np.asarray(positions)
    geno = (raw >= 0.5).astype(np.int8)
    ppos, AF = load_allclass_af(a.chrom, "gen9")      # ppos sorted asc, AF pools x records

    counts = pd.read_csv(f"{BR}/block_haplotype_counts_clq0.9.csv")
    counts = counts[counts.chrom == a.chrom]
    gated = counts[counts.n_eff > a.gate].reset_index(drop=True)
    if a.max_blocks:
        gated = gated.iloc[:a.max_blocks]
    print(f"[{a.chrom}] {len(gated)} blocks with n_eff>{a.gate} to cluster", flush=True)

    def pool_pc1ve(sig_positions):
        """PC1-VE on gen9 pools for the gen9 records at these founder positions."""
        if len(sig_positions) < 2:
            return np.nan, 0
        ix = np.searchsorted(ppos, sig_positions)
        ix = ix[(ix < len(ppos))]
        ix = ix[ppos[ix] == sig_positions[:len(ix)]] if len(ix) else ix
        # robust exact-match gather
        keep = []
        for p in sig_positions:
            j = np.searchsorted(ppos, p)
            if j < len(ppos) and ppos[j] == p:
                keep.append(j)
        if len(keep) < 2:
            return np.nan, len(keep)
        return pc1_var_explained(AF[:, keep]), len(keep)

    rows = []
    for bi in range(len(gated)):
        b = gated.iloc[bi]
        lo = int(np.searchsorted(positions, b.start))
        hi = int(np.searchsorted(positions, b.end, side="right"))
        sub = geno[:, lo:hi]
        bpos = positions[lo:hi]
        if sub.shape[1] < 2:
            continue
        # whole-block PC1-VE on pools (the baseline we want to beat)
        block_ve, block_nmatch = pool_pc1ve(bpos)
        try:
            labels, d_uniq = cluster_founders(sub)
        except Exception as e:
            n_fail = locals().get("n_fail", 0) + 1
            if n_fail <= 5:
                print(f"  block {bi} cluster fail: {e}", flush=True)
            continue
        k = len(np.unique(labels))
        for c in np.unique(labels):
            inmask = labels == c
            n_in = int(inmask.sum())
            mean_in = sub[inmask].mean(0)
            mean_out = sub[~inmask].mean(0) if n_in < sub.shape[0] else np.zeros(sub.shape[1])
            sig = (mean_in >= SIG_PRESENT) & (mean_in - mean_out >= SIG_DIFF)
            sigpos = bpos[sig]
            cve, nmatch = pool_pc1ve(sigpos)
            rows.append(dict(block=bi, start=int(b.start), end=int(b.end),
                             n_eff=round(float(b.n_eff), 2), nvar=int(sub.shape[1]),
                             k=k, d_uniq=int(d_uniq), cluster=int(c),
                             cluster_freq=round(n_in / sub.shape[0], 3),
                             n_sig=int(sig.sum()), n_sig_matched=nmatch,
                             block_ve=round(block_ve, 3) if np.isfinite(block_ve) else np.nan,
                             cluster_ve=round(cve, 3) if np.isfinite(cve) else np.nan))
        if (bi + 1) % 2000 == 0:
            print(f"  ...{bi+1}/{len(gated)} blocks", flush=True)

    R = pd.DataFrame(rows)
    out = f"{BR}/{chrlc}_cluster_pc1ve_gate{a.gate}.csv"
    R.to_csv(out, index=False)

    # ---- summary: does per-cluster VE beat block VE? ----
    valid = R[np.isfinite(R.cluster_ve)]
    real = valid[(valid.n_sig_matched >= 2) & (valid.cluster_freq >= 2/231)]  # drop private singletons
    print(f"\n=== {a.chrom} per-cluster PC1-VE (gate n_eff>{a.gate}) ===")
    print(f"clustered blocks: {R.block.nunique()};  clusters scored: {len(valid)};  "
          f"non-singleton clusters w/ >=2 matched sig vars: {len(real)}")
    bl = R.groupby('block').block_ve.first().dropna()
    print(f"\nWHOLE-BLOCK PC1-VE (baseline):  median {bl.median():.3f}  "
          f"%>=0.7 {100*(bl>=0.7).mean():.0f}%")
    print(f"PER-CLUSTER PC1-VE (signature): median {real.cluster_ve.median():.3f}  "
          f"%>=0.7 {100*(real.cluster_ve>=0.7).mean():.0f}%")
    print("\n  per-cluster median VE by block effective-diversity bin:")
    for lo, hi, lab in [(2, 4, "mild  2-4"), (4, 7, "mod   4-7"),
                        (7, 12, "high  7-12"), (12, 1e9, "vhigh >12")]:
        m = real[(real.n_eff > lo) & (real.n_eff <= hi)]
        if len(m):
            print(f"    {lab}: clusters {len(m):5d}  median cluster_ve {m.cluster_ve.median():.3f}  "
                  f"%>=0.7 {100*(m.cluster_ve>=0.7).mean():.0f}%")
    print(f"\nwrote -> {out}")


if __name__ == "__main__":
    main()
