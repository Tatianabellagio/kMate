#!/usr/bin/env python
"""Two-decision frontier data: for a given CLQcut block map (Chr1), emit per-BLOCK
and per-CLUSTER PC1-VE so the notebook can synthesize the (CLQcut x clustering-gate)
tradeoff between coherence (per-unit PC1-VE) and #test-units.

Per block: n_eff (effective # founder haplotypes) + whole-block PC1-VE on gen9 pools.
Per cluster (only for blocks with n_eff>CLUSTER_GATE): HapFM xmeans cluster ->
signature variants -> PC1-VE on gen9 pools.

A "gate G" config (synthesized in the notebook) = block-level unit if n_eff<=G, else
that block's cluster units. Run for clqcut 0.5/0.7/0.9 to overlay frontiers.

Run in kmate env. Usage: block_unit_frontier.py --clqcut 0.9 [--chrom Chr1] [--cluster-gate 2]
"""
import argparse, os, sys
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from recompute_blocks import build_common_matrix
from eval_block_coherence import load_allclass_af, pc1_var_explained
from block_haplotype_counts import hap_counts
from block_cluster_pc1ve import cluster_founders, SIG_PRESENT, SIG_DIFF   # shim applied on import

BR = os.environ.get("BLOCKS_DIR", "analysis/grenenet_gea/blocks_recompute")
MAF = 0.05
MINCF = float(os.environ.get("MINCF", "0.5"))   # must match the block-build filter


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--chrom", default="Chr1")
    ap.add_argument("--clqcut", required=True)
    ap.add_argument("--cluster-gate", type=float, default=2.0)
    a = ap.parse_args()
    chrlc = a.chrom.lower()

    _, raw, positions = build_common_matrix(f"panel/arch3/{chrlc}/var_pa_231_arch3_{chrlc}", MAF, MINCF)
    positions = np.asarray(positions)
    geno = (raw >= 0.5).astype(np.int8)
    ppos, AF = load_allclass_af(a.chrom, "gen9")
    bl = pd.read_csv(f"{BR}/{chrlc}_clq{a.clqcut}_blocks_clq{a.clqcut}.tsv", sep="\t")
    print(f"[{a.chrom} clq{a.clqcut}] {len(bl)} blocks", flush=True)

    def pool_pc1ve(sig_positions):
        keep = []
        for p in sig_positions:
            j = np.searchsorted(ppos, p)
            if j < len(ppos) and ppos[j] == p:
                keep.append(j)
        if len(keep) < 2:
            return np.nan, len(keep)
        return pc1_var_explained(AF[:, keep]), len(keep)

    brows, crows = [], []
    for bi, (s, e, nv) in enumerate(zip(bl.start_pos, bl.end_pos, bl.n_variants)):
        lo = int(np.searchsorted(positions, s)); hi = int(np.searchsorted(positions, e, side="right"))
        sub = geno[:, lo:hi]; bpos = positions[lo:hi]
        if sub.shape[1] < 2:
            continue
        nd, neff, ng = hap_counts(sub)
        bve, bm = pool_pc1ve(bpos)
        brows.append((a.clqcut, int(s), int(e), int(sub.shape[1]), nd, round(neff, 3),
                      round(bve, 4) if np.isfinite(bve) else np.nan, bm))
        if neff > a.cluster_gate:
            try:
                labels, d = cluster_founders(sub)
            except Exception:
                continue
            k = len(np.unique(labels))
            for c in np.unique(labels):
                inmask = labels == c; n_in = int(inmask.sum())
                mean_in = sub[inmask].mean(0)
                mean_out = sub[~inmask].mean(0) if n_in < sub.shape[0] else np.zeros(sub.shape[1])
                sig = (mean_in >= SIG_PRESENT) & (mean_in - mean_out >= SIG_DIFF)
                cve, nm = pool_pc1ve(bpos[sig])
                crows.append((a.clqcut, int(s), int(e), round(neff, 3), k, int(c),
                              round(n_in / sub.shape[0], 3), int(sig.sum()), nm,
                              round(cve, 4) if np.isfinite(cve) else np.nan))
        if (bi + 1) % 3000 == 0:
            print(f"  ...{bi+1}/{len(bl)}", flush=True)

    B = pd.DataFrame(brows, columns=["clqcut", "start", "end", "nvar", "n_distinct",
                                     "n_eff", "block_ve", "block_nmatch"])
    C = pd.DataFrame(crows, columns=["clqcut", "start", "end", "n_eff", "k", "cluster",
                                     "cluster_freq", "n_sig", "n_sig_matched", "cluster_ve"])
    B.to_csv(f"{BR}/{chrlc}_frontier_blocks_clq{a.clqcut}.csv", index=False)
    C.to_csv(f"{BR}/{chrlc}_frontier_clusters_clq{a.clqcut}.csv", index=False)
    print(f"[{a.chrom} clq{a.clqcut}] blocks={len(B)} clusters={len(C)} -> frontier CSVs", flush=True)


if __name__ == "__main__":
    main()
