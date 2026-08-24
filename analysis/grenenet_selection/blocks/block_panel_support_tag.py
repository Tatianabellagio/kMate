#!/usr/bin/env python
"""Per-variant haplotype tagging + panel support (Chr1, locked CLQ0.9 + gate2).

For each block variant (the unique-position set that defines the blocks), record:
  - founder call-rate (panel support, from var_called)
  - the haplotype-CLUSTER it tags = argmax founder-enrichment over the block's
    HapFM-xmeans clusters (n_eff>2 blocks) / the single block haplotype (n_eff<=2).
This lets the notebook ask: if we drop low-panel-support variants at threshold X,
do they still have a good-support variant tagging the SAME haplotype in the block?

Run in kmate env. Output: chr1_panel_support_tag.csv
"""
import os, sys, argparse
import numpy as np
import pandas as pd
import scipy.sparse as sp
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from block_haplotype_counts import hap_counts
from block_cluster_pc1ve import cluster_founders   # numpy.warnings shim applied on import

BR = "analysis/grenenet_selection/blocks_recompute"
MAF, MINCF = 0.05, 0.5
GATE = 2.0


def load(pref):
    """Replicate recompute_blocks.build_common_matrix selection, also returning the
    per-variant founder call-rate aligned to the kept unique positions."""
    vp = sp.load_npz(f"{pref}.var_pa.npz").toarray().astype(np.float64)
    vc = sp.load_npz(f"{pref}.var_called.npz").toarray().astype(np.float64)
    pos = np.load(f"{pref}.meta.npz", allow_pickle=True)["pos"].astype(np.int64)
    F = vp.shape[0]
    n_called = vc.sum(0); n_alt = vp.sum(0)
    af = np.divide(n_alt, n_called, out=np.full_like(n_alt, np.nan), where=n_called > 0)
    keep = (af > MAF) & (af < 1 - MAF) & (n_called >= MINCF * F)
    order = np.where(keep)[0]
    seen = set(); uniq = []
    for c in order:
        p = int(pos[c])
        if p not in seen:
            seen.add(p); uniq.append(c)
    uniq = np.array(uniq)
    raw = vp[:, uniq].copy(); miss = vc[:, uniq] == 0; afc = af[uniq]
    raw[miss] = np.repeat(afc[None, :], F, axis=0)[miss]
    geno = (raw >= 0.5).astype(np.int8)
    callrate = n_called[uniq] / F
    return geno, callrate, np.array([int(p) for p in pos[uniq]])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--chrom", default="Chr1")
    a = ap.parse_args()
    chrlc = a.chrom.lower()
    pref = f"panel/arch3/{chrlc}/var_pa_231_arch3_{chrlc}"
    print(f"[{a.chrom}] loading founder matrix + call-rate ...", flush=True)
    geno, callrate, positions = load(pref)
    bl = pd.read_csv(f"{BR}/{chrlc}_clq0.9_blocks_clq0.9.tsv", sep="\t")
    print(f"{len(bl)} blocks; {len(positions)} block variants", flush=True)

    rows = []
    for bi, (s, e) in enumerate(zip(bl.start_pos, bl.end_pos)):
        lo = int(np.searchsorted(positions, s)); hi = int(np.searchsorted(positions, e, side="right"))
        if hi - lo < 1:
            continue
        sub = geno[:, lo:hi]; bpos = positions[lo:hi]; cr = callrate[lo:hi]
        m = sub.shape[1]
        nd, neff, ng = hap_counts(sub) if m >= 2 else (1, 1.0, 1)
        if m >= 2 and neff > GATE:
            try:
                labels, _ = cluster_founders(sub)
            except Exception:
                labels = np.zeros(sub.shape[0], dtype=int)
            cls = np.unique(labels)
            # per-cluster founder-mean per variant -> tag = argmax enrichment
            means = np.stack([sub[labels == c].mean(0) for c in cls], axis=0)  # ncl x m
            best = means.argmax(0)
            tag = cls[best]
            srt = np.sort(means, axis=0)
            strength = srt[-1] - (srt[-2] if means.shape[0] > 1 else srt[-1] * 0)  # top - 2nd
        else:
            tag = np.zeros(m, dtype=int)
            strength = np.full(m, np.nan)
        for j in range(m):
            rows.append((int(bpos[j]), bi, round(float(neff), 3), int(tag[j]),
                         round(float(strength[j]), 3) if np.isfinite(strength[j]) else np.nan,
                         round(float(cr[j]), 4)))
        if (bi + 1) % 5000 == 0:
            print(f"  ...{bi+1}/{len(bl)}", flush=True)

    df = pd.DataFrame(rows, columns=["pos", "block", "n_eff", "tag_cluster",
                                     "tag_strength", "call_rate"])
    out = f"{BR}/{chrlc}_panel_support_tag.csv"
    df.to_csv(out, index=False)
    print(f"wrote {len(df)} variant rows -> {out}", flush=True)
    print(f"  call_rate: median {df.call_rate.median():.3f}; "
          f"distinct (block,tag) groups: {df.groupby(['block','tag_cluster']).ngroups}")


if __name__ == "__main__":
    main()
