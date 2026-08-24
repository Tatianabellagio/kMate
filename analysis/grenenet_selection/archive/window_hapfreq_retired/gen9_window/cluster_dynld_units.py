#!/usr/bin/env python
"""Apply the HapFM xmeans haplotype-clustering (yesterday's technique, faithful reuse
from block_cluster_pc1ve.py) to the PRODUCTION dynld K500 units — the units we ran
window-mode kMate on. For each unit, cluster the 231 founders by their block-haplotype;
each cluster = a distinct founder-haplotype whose per-sample frequency = sum of the
window h over its founders. Founder order matches the h npz (verified), so labels align.

Outputs (results/grenenet_gea/gen9_window/clusters/):
  {chrlc}_clusters.csv   one row per unit x cluster: chrom,unit_idx,start,end,nvar,
                         n_eff,n_uniq,cluster,n_founders,cluster_freq(founding),n_sig
  {chrlc}_labels.npz     labels [n_units x 231] int8 (founder->cluster per unit),
                         start[n_units], end[n_units], n_eff[n_units]

Run in kmate env (pyclustering). Usage: cluster_dynld_units.py Chr1
"""
import os, sys, warnings
import numpy as np, pandas as pd
if not hasattr(np, "warnings"):
    np.warnings = warnings                      # pyclustering 0.10 / numpy>=1.24 shim
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from recompute_blocks import build_common_matrix

MAF, MINCF = 0.05, 0.9
HAPFM_KMAX = 7
SIG_PRESENT, SIG_DIFF = 0.5, 0.5
MD = "results/grenenet_gea/blocks_mcf90"
OUT = "results/grenenet_gea/gen9_window/clusters"

def xmeans_clustering(array):
    from pyclustering.cluster.xmeans import xmeans
    from pyclustering.cluster.center_initializer import kmeans_plusplus_initializer
    ic = kmeans_plusplus_initializer(array, 2).initialize()
    inst = xmeans(array, ic, 30); inst.process()
    labels = [0] * len(array)
    for i, grp in enumerate(inst.get_clusters()):
        for j in grp:
            labels[j] = i
    return np.asarray(labels)

def cluster_founders(sub):
    """sub: 231 x m int 0/1 -> per-founder cluster label, n_unique_haplotypes."""
    F = sub.shape[0]
    uniq, inv = np.unique(sub, axis=0, return_inverse=True)
    d = uniq.shape[0]
    if d < 2:
        return np.zeros(F, int), d
    ul = np.arange(d) if d < HAPFM_KMAX else xmeans_clustering([r.astype(float).tolist() for r in uniq])
    return ul[inv], d

def main():
    chrom = sys.argv[1]; chrlc = chrom.lower()
    _, raw, positions = build_common_matrix(f"panel/arch3/{chrlc}/var_pa_231_arch3_{chrlc}", MAF, MINCF)
    positions = np.asarray(positions); geno = (raw >= 0.5).astype(np.int8)
    units = pd.read_csv(f"{MD}/{chrlc}_units_dynld_K500.tsv", sep="\t")
    os.makedirs(OUT, exist_ok=True)
    rows = []; labels_all = np.full((len(units), 231), -1, np.int8)
    neff_all = np.zeros(len(units)); st = units.start_pos.to_numpy(); en = units.end_pos.to_numpy()
    for ui in range(len(units)):
        lo = int(np.searchsorted(positions, st[ui]))
        hi = int(np.searchsorted(positions, en[ui], side="right"))
        sub = geno[:, lo:hi]
        if sub.shape[1] < 1:
            continue
        uniq, inv = np.unique(sub, axis=0, return_inverse=True)
        cnt = np.bincount(inv); fr = cnt / 231.0
        n_eff = float(1.0 / np.sum(fr ** 2)); neff_all[ui] = n_eff
        labels, d = cluster_founders(sub)
        labels_all[ui] = labels
        for c in np.unique(labels):
            inm = labels == c; n_in = int(inm.sum())
            mean_in = sub[inm].mean(0)
            mean_out = sub[~inm].mean(0) if n_in < sub.shape[0] else np.zeros(sub.shape[1])
            n_sig = int(((mean_in >= SIG_PRESENT) & (mean_in - mean_out >= SIG_DIFF)).sum())
            rows.append((chrom, ui, int(st[ui]), int(en[ui]), int(sub.shape[1]),
                         round(n_eff, 3), int(d), int(c), n_in, round(n_in / 231.0, 4), n_sig))
        if (ui + 1) % 1000 == 0:
            print(f"  {chrom}: {ui+1}/{len(units)}", flush=True)
    df = pd.DataFrame(rows, columns=["chrom", "unit_idx", "start", "end", "nvar", "n_eff",
                                     "n_uniq", "cluster", "n_founders", "cluster_freq", "n_sig"])
    df.to_csv(f"{OUT}/{chrlc}_clusters.csv", index=False)
    np.savez(f"{OUT}/{chrlc}_labels.npz", labels=labels_all, start=st, end=en, n_eff=neff_all)
    print(f"{chrom}: {len(units)} units, {len(df)} clusters; "
          f"n_eff>2 units {int((neff_all>2).sum())} -> {OUT}/{chrlc}_clusters.csv", flush=True)

if __name__ == "__main__":
    main()
