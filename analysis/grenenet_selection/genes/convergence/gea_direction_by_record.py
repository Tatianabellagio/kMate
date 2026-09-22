#!/usr/bin/env python
"""Warm or cold, by region: the direction table behind notebooks/gea_hit_landscape.ipynb §7.

For every tested non-SNP record (gea_hit_landscape.csv.gz, split records collapsed as in the
notebook) this takes the GEA's own signed result (lfmm_signed/, the identical LFMM with B and z
kept; p reproduced exactly) on the nine temperature axes, picks the axis where the record is
MOST significant, and records the sign there: z > 0 = the ALT allele (for a deletion, the
deletion) rose more in warmer gardens. One sign per record, so a hit on several correlated
temperature axes is counted once. `temp_hit` = a GEA hit on at least one temperature axis.

Temperature axes: bio1, bio5, bio6, bio8, bio9, bio10, bio11, pc1 (-> +bio1), pc3 (-> +bio10).
Writes results/gea_direction_by_record.csv.gz. env: kmate. Compute node.
"""
import os, sys
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import gea_direction as GD                                       # noqa: E402

OUT = f"{HERE}/results"
K4 = ["chrom", "pos", "ref_len", "alt_len"]


def main():
    U = pd.read_csv(f"{OUT}/gea_hit_landscape.csv.gz")
    U = U.sort_values(["hit", "min_p"], ascending=[False, True]).drop_duplicates(K4)
    U = U[U.kind.isin(["deletion", "insertion"])].copy()
    U["mafr"] = U.MAF.round(6)
    best_p = pd.Series(np.inf, index=U.index); best_z = pd.Series(np.nan, index=U.index)
    best_ax = pd.Series("", index=U.index)
    for ax in GD.TEMP:
        S = pd.read_csv(f"{HERE}/lfmm_signed/lfmm_nonsnp_gen9_{ax}_signed.csv.gz")
        S["mafr"] = S.MAF.round(6)
        m = U[K4 + ["mafr"]].reset_index().merge(S[K4 + ["mafr", "pval", "z"]].drop_duplicates(K4 + ["mafr"]),
                                                 on=K4 + ["mafr"], how="left").set_index("index")
        better = m.pval < best_p
        best_p[better], best_z[better], best_ax[better] = m.pval[better], m.z[better], ax
        print(f"  {ax}", flush=True)
    U["temp_axis"], U["temp_p"], U["temp_z"] = best_ax, best_p, best_z
    U["warm"] = U.temp_z > 0
    U["temp_hit"] = U.hit_axes.fillna("").str.split(",").map(lambda s: any(a in GD.TEMP for a in s))
    cols = K4 + ["cls", "kind", "size", "MAF", "block", "tier_1kb", "hit", "temp_hit", "temp_axis",
                 "temp_p", "temp_z", "warm"]
    U[cols].to_csv(f"{OUT}/gea_direction_by_record.csv.gz", index=False)
    print(f"{len(U):,} records; {int(U.temp_hit.sum())} hits on a temperature axis; "
          f"warm share among them {U[U.temp_hit].warm.mean():.3f}, among all tested {U.warm.mean():.3f}")


if __name__ == "__main__":
    main()
