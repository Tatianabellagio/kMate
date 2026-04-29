#!/usr/bin/env python3
"""
Use the GrENE-net seed mix recipe (per-ecotype seed proportions) to compute a
"panel-restricted truth" SV frequency:

    truth_panel(v) = Σ_{e ∈ 80 panel ∩ recipe} seed_prop(e) · G_SV(e, v)

This is the lower bound of the true seedmix SV freq (we only count carriers
whose SV genotype we know). hapFIRE projection should track this closely if
hapFIRE recovers the correct ecotype frequencies on the 80 founders. freqk
tries to estimate the true freq, which includes mass on the unobserved 151.

Outputs a 3-way comparison: truth_panel vs hapfire_proj vs freqk.
"""
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path("/home/tbellagio/scratch/hapfire_sv")
SAMPLE_IDS = [f"SEEDMIX_S{i}" for i in range(1, 9)]


def main():
    t0 = time.time()
    # Load founder × SV matrix
    G = pd.read_parquet(ROOT / "data" / "founder_sv_matrix.parquet")
    meta = pd.read_parquet(ROOT / "data" / "founder_sv_meta.parquet")
    ecotype_cols = [c for c in G.columns if c not in ("chrom", "pos", "alt_idx")]
    df_e = G[ecotype_cols]
    G_collapsed = df_e.T.groupby(df_e.columns, axis=0).max().T
    G_mat = G_collapsed.to_numpy(dtype=np.float32)
    ecotype_cols = G_collapsed.columns.tolist()
    print(f"founder×SV matrix: {G_mat.shape}, {len(ecotype_cols)} unique ecotypes")

    # Load seed mix recipe → seed_prop for each ecotype
    recipe = pd.read_csv(ROOT / "data" / "seedmix_recipe.tsv", sep="\t", dtype={"ID": str})
    seed_prop = dict(zip(recipe["ID"], recipe["seed_prop"]))
    panel_props = np.array([seed_prop.get(e, 0.0) for e in ecotype_cols], dtype=np.float64)
    print(f"recipe coverage on panel: total mass = {panel_props.sum():.4f}")
    # The 80 panel ecotypes account for `panel_props.sum()` of the seedmix.

    # truth_panel = G_SV @ seed_prop  (vector of length n_sv)
    truth_panel = G_mat @ panel_props
    print(f"truth_panel quantiles: {np.quantile(truth_panel, [0.5, 0.9, 0.99, 1.0]).round(4)}")

    # Load freqk and hapfire projection
    freqk = pd.read_parquet("/home/tbellagio/scratch/freqk_gr/results/seedmix_p_wide.parquet")
    proj  = pd.read_parquet(ROOT / "results" / "hapfire_proj_svs_seedmix.parquet")
    G_keys = G[["chrom", "pos", "alt_idx"]].copy()
    G_keys["chrom"] = G_keys["chrom"].astype(str)
    G_keys["pos"]   = G_keys["pos"].astype(np.int32)
    G_keys["alt_idx"] = G_keys["alt_idx"].astype(np.int8)
    G_keys["truth_panel"] = truth_panel
    G_keys = G_keys.merge(meta[["chrom","pos","alt_idx","var_type","sv_size","ac_in_panel"]],
                          on=["chrom","pos","alt_idx"])

    # Take just S1 sample first for plotting; build a tall df with all comparisons
    rows = []
    for sid in SAMPLE_IDS:
        sub = G_keys.copy()
        sub["sample_id"] = sid
        # hapfire_proj
        proj_sub = proj[["chrom","pos","alt_idx",sid]].rename(columns={sid:"hapfire_proj_af"})
        proj_sub["chrom"] = proj_sub["chrom"].astype(str)
        sub = sub.merge(proj_sub, on=["chrom","pos","alt_idx"], how="left")
        # freqk
        freqk_sub = freqk[freqk["alt_idx"].isin([1,2,3,4,5,6,7,8,9])][["chrom","pos","alt_idx",sid]].rename(columns={sid:"freqk_af"})
        freqk_sub["chrom"] = freqk_sub["chrom"].astype(str)
        sub = sub.merge(freqk_sub, on=["chrom","pos","alt_idx"], how="left")
        rows.append(sub)
    full = pd.concat(rows, ignore_index=True)
    full.to_csv(ROOT / "results" / "seedmix_3way_comparison.tsv.gz",
                sep="\t", index=False, compression="gzip", float_format="%.4f")
    print(f"wrote: {ROOT / 'results' / 'seedmix_3way_comparison.tsv.gz'} ({len(full):,} rows)")

    # Per-sample correlations: hapfire_proj vs truth_panel and freqk vs truth_panel and hapfire vs freqk
    print("\nPer-sample correlations (truth_panel uses recipe; ascertainment = panel coverage):")
    print("                     n_with_freqk  hapf_vs_truth  freqk_vs_truth  hapf_vs_freqk")
    for sid in SAMPLE_IDS:
        sub = full[full["sample_id"] == sid]
        # All SVs in panel: have truth_panel + hapfire_proj for all
        has_freqk = sub["freqk_af"].notna()
        x_truth = sub["truth_panel"].to_numpy(dtype=float)
        y_hapf  = sub["hapfire_proj_af"].to_numpy(dtype=float)
        # hapf vs truth (over ALL panel SVs, not just freqk-overlap)
        ok_ht = np.isfinite(x_truth) & np.isfinite(y_hapf)
        r_ht = np.corrcoef(x_truth[ok_ht], y_hapf[ok_ht])[0,1] if ok_ht.sum()>1 else float("nan")

        # restrict to freqk-overlap for the other two
        sub2 = sub[has_freqk]
        x = sub2["freqk_af"].to_numpy(dtype=float)
        y = sub2["hapfire_proj_af"].to_numpy(dtype=float)
        z = sub2["truth_panel"].to_numpy(dtype=float)
        ok_fz = np.isfinite(x) & np.isfinite(z)
        r_fz = np.corrcoef(x[ok_fz], z[ok_fz])[0,1] if ok_fz.sum()>1 else float("nan")
        ok_hf = np.isfinite(x) & np.isfinite(y)
        r_hf = np.corrcoef(x[ok_hf], y[ok_hf])[0,1] if ok_hf.sum()>1 else float("nan")
        print(f"  {sid}:  n={int(has_freqk.sum())}  ht={r_ht:.3f}  ft={r_fz:.3f}  hf={r_hf:.3f}")

    # Diagnostic: is hapfire_proj ≈ truth_panel? If yes, hapFIRE recovered ecotype freq accurately
    sid = "SEEDMIX_S1"
    sub = full[full["sample_id"] == sid].copy()
    sub["abs_err_hapfire_vs_truth"] = (sub["hapfire_proj_af"] - sub["truth_panel"]).abs()
    sub["abs_err_freqk_vs_truth"]   = (sub["freqk_af"] - sub["truth_panel"]).abs()
    print(f"\nFor S1: comparing each method vs truth_panel (which is the LOWER bound — only 80-panel mass)")
    print(f"  hapfire_proj vs truth_panel: MAE = {sub['abs_err_hapfire_vs_truth'].mean():.4f}, "
          f"R = {sub[['hapfire_proj_af','truth_panel']].corr().iloc[0,1]:.3f}")
    print(f"  freqk vs truth_panel       : MAE = {sub['abs_err_freqk_vs_truth'].mean():.4f}, "
          f"R = {sub[['freqk_af','truth_panel']].corr().iloc[0,1]:.3f}")
    print(f"  (truth_panel mean = {sub['truth_panel'].mean():.4f}, "
          f"hapf mean = {sub['hapfire_proj_af'].mean():.4f}, "
          f"freqk mean = {sub['freqk_af'].mean():.4f})")
    print()
    print(f"Expected behavior:")
    print(f"  - hapfire_proj_af should match truth_panel almost perfectly (~no error)")
    print(f"    if hapFIRE recovered uniform-ish ecotype freq for the 80 panel founders.")
    print(f"  - freqk recovers the TRUE freq (which includes mass on unobserved 151).")
    print(f"  - So freqk > hapfire_proj typically, by a factor reflecting unobserved-founder mass.")


if __name__ == "__main__":
    main()
