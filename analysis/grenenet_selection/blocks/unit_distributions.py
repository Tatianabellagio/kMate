#!/usr/bin/env python
"""Per-unit metrics for the final dynamic-LD unit map: n_variants, panel_kmers (from the map),
mean within-unit founder LD (r2), and PC1-VE on the gen9 pools. One row per unit.
Output: analysis/grenenet_selection/blocks_mcf90/final_units_dynld_K500_metrics.csv
"""
import os, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from recompute_blocks import build_common_matrix
from eval_block_coherence import load_allclass_af, pc1_var_explained

MAF, MINCF = 0.05, 0.9
BR = "analysis/grenenet_selection/blocks_mcf90"
CHROMS = ["Chr1", "Chr2", "Chr3", "Chr4", "Chr5"]


def mean_within_r2(std_sub):
    m = std_sub.shape[1]
    if m < 2:
        return np.nan
    C = (std_sub.T @ std_sub) / std_sub.shape[0]    # m x m correlation (std cols)
    r2 = C ** 2
    iu = np.triu_indices(m, k=1)
    return float(np.nanmean(r2[iu]))


def main():
    rows = []
    for chrom in CHROMS:
        chrlc = chrom.lower()
        print(f"[{chrom}] founder matrix + gen9 pools ...", flush=True)
        std, _, positions = build_common_matrix(f"panel/arch3/{chrlc}/var_pa_231_arch3_{chrlc}", MAF, MINCF)
        positions = np.asarray(positions)
        ppos, AF = load_allclass_af(chrom, "gen9")
        U = pd.read_csv(f"{BR}/{chrlc}_units_dynld_K500.tsv", sep="\t")
        for _, r in U.iterrows():
            lo = int(np.searchsorted(positions, r.start_pos)); hi = int(np.searchsorted(positions, r.end_pos, side="right"))
            mr2 = mean_within_r2(std[:, lo:hi]) if hi - lo >= 2 else np.nan
            pl = np.searchsorted(ppos, r.start_pos); ph = np.searchsorted(ppos, r.end_pos, side="right")
            ve = pc1_var_explained(AF[:, pl:ph]) if ph - pl >= 2 else np.nan
            rows.append((chrom, int(r.start_pos), int(r.end_pos), int(r.n_variants),
                         int(r.panel_kmers), int(ph - pl), bool(r.covered),
                         round(mr2, 4) if np.isfinite(mr2) else np.nan,
                         round(ve, 4) if np.isfinite(ve) else np.nan))
        print(f"[{chrom}] {len(U)} units done", flush=True)
    D = pd.DataFrame(rows, columns=["chrom", "start", "end", "n_variants", "panel_kmers",
                                    "n_gen9", "covered", "mean_r2", "pc1_ve"])
    out = f"{BR}/final_units_dynld_K500_metrics.csv"
    D.to_csv(out, index=False)
    print(f"\nwrote {len(D)} units -> {out}")
    for col in ["n_variants", "panel_kmers", "mean_r2", "pc1_ve"]:
        v = D[col].dropna()
        print(f"  {col:>12}: median {v.median():.3f}  q25 {v.quantile(.25):.3f}  q75 {v.quantile(.75):.3f}")


if __name__ == "__main__":
    main()
