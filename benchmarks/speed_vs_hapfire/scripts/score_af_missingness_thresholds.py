#!/usr/bin/env python3
"""kMate SNP AF R2/RMSE at three reference-panel completeness thresholds
(all records / n_called==231 exact zero-missing / n_called>=208 <=10% missing),
alongside hapFIRE's fixed zero-missingness RMSE/R2 from hapfire_vs_kmate_af_table.tsv.

This is the detailed missingness-breakdown behind the >=90%-called filter that
score_snp_vs_nonsnp.py applies by default: it quantifies how much the fair-
comparison filter changes kMate's AF accuracy vs hapFIRE across the grid
(archive/plot_consolidation_2026-07-21/scripts/plot_h_and_af_90pctcomplete.py
plotted this breakdown; superseded by hapfire_vs_kmate_af_accuracy.png /
hapfire_vs_kmate_h_accuracy.png, which already apply the filter). Reuses
already-computed kMate tsv outputs + recomb_truth_raw.tsv.gz + the hapFIRE af
table -- no new compute.

Run with the `basic` env:
    /global/home/users/tbellg/miniforge3/envs/basic/bin/python \
        benchmarks/speed_vs_hapfire/scripts/score_af_missingness_thresholds.py
"""
import os
import numpy as np
import pandas as pd

ROOT = "/global/scratch/users/tbellg/kmate"
P231 = f"{ROOT}/benchmarks/p231"
OUT = f"{ROOT}/benchmarks/speed_vs_hapfire/results"
os.makedirs(OUT, exist_ok=True)

POOL_SIZES = [2, 5, 20, 50, 150, 231]
DEPTHS = [1, 10]
SEEDS = [42, 43, 44, 45, 46]
REL_EPS = 1e-6
FULL_PANEL = 231
THRESH_90 = int(np.ceil(0.9 * FULL_PANEL))  # 208


def r2_rmse(est, truth):
    d = est - truth
    ss_tot = np.sum((truth - truth.mean()) ** 2)
    degenerate = ss_tot < (REL_EPS * truth.mean()) ** 2 * len(truth)
    r2 = np.nan if degenerate else 1 - np.sum(d ** 2) / ss_tot
    rmse = float(np.sqrt(np.mean(d ** 2)))
    return float(r2), rmse


rows = []
for n in POOL_SIZES:
    for cov in DEPTHS:
        for seed in SEEDS:
            sim = f"{P231}/sims/cov{cov}_n{n}_g0_s{seed}_hotspots_p231_chr1"
            kmate_dir = f"{P231}/results/kmate_chrom_poolsize_depth/n{n}_cov{cov}_s{seed}"
            kmate_tsv = f"{kmate_dir}/p231_chrom_psd_n{n}_cov{cov}_s{seed}.tsv"
            truth_path = f"{sim}/recomb_truth_raw.tsv.gz"
            if not (os.path.exists(kmate_tsv) and os.path.exists(truth_path)):
                print(f"  [skip] n={n} cov={cov} s={seed}: missing output")
                continue

            KEYS = ["chrom", "pos", "ref_len", "alt_len"]
            est = pd.read_csv(kmate_tsv, sep="\t")
            truth = pd.read_csv(truth_path, sep="\t").dropna(subset=["truth_af"])
            truth = truth.copy(); truth["occ"] = truth.groupby(KEYS).cumcount()
            est = est.copy(); est["occ"] = est.groupby(KEYS).cumcount()
            m = truth.merge(est[KEYS + ["occ", "alt_freq", "n_called"]], on=KEYS + ["occ"], how="inner")

            is_snp = (m.ref_len == 1) & (m.alt_len == 1)
            snp = m.loc[is_snp]

            # (a) unrestricted
            a_r2, a_rmse = r2_rmse(snp["alt_freq"].values, snp["truth_af"].values)
            a_n = len(snp)
            # (b) exact zero missingness
            b_mask = snp["n_called"] == FULL_PANEL
            b_r2, b_rmse = r2_rmse(snp.loc[b_mask, "alt_freq"].values, snp.loc[b_mask, "truth_af"].values)
            b_n = int(b_mask.sum())
            # (c) <=10% missing (n_called >= 208)
            c_mask = snp["n_called"] >= THRESH_90
            c_r2, c_rmse = r2_rmse(snp.loc[c_mask, "alt_freq"].values, snp.loc[c_mask, "truth_af"].values)
            c_n = int(c_mask.sum())

            rows.append(dict(N=n, depth=cov, seed=seed,
                             all_n=a_n, all_R2=a_r2, all_RMSE=a_rmse,
                             zero_n=b_n, zero_frac=b_n / a_n if a_n else np.nan,
                             zero_R2=b_r2, zero_RMSE=b_rmse,
                             p90_n=c_n, p90_frac=c_n / a_n if a_n else np.nan,
                             p90_R2=c_r2, p90_RMSE=c_rmse))
            print(f"n={n:>3} cov={cov:>2}x s={seed}: all n={a_n:,} R2={a_r2:.4f} RMSE={a_rmse:.4f} | "
                  f"zero-miss n={b_n:,} ({b_n/a_n:.1%}) R2={b_r2:.4f} RMSE={b_rmse:.4f} | "
                  f"<=10%miss n={c_n:,} ({c_n/a_n:.1%}) R2={c_r2:.4f} RMSE={c_rmse:.4f}")

df = pd.DataFrame(rows)
pc_path = f"{OUT}/kmate_af_missingness_percondition.tsv"
df.to_csv(pc_path, sep="\t", index=False)
print(f"\nwrote {pc_path} ({len(df)} rows)")

# summary: mean over seeds, per N per depth, vs hapfire reference
hf = pd.read_csv(f"{OUT}/hapfire_vs_kmate_af_table.tsv", sep="\t")
summ = df.groupby(["N", "depth"]).agg(
    kmate_RMSE_all=("all_RMSE", "mean"),
    kmate_RMSE_zeromissing=("zero_RMSE", "mean"),
    kmate_RMSE_90pctcomplete=("p90_RMSE", "mean"),
    kmate_R2_all=("all_R2", "mean"),
    kmate_R2_zeromissing=("zero_R2", "mean"),
    kmate_R2_90pctcomplete=("p90_R2", "mean"),
    zero_frac_mean=("zero_frac", "mean"),
    p90_frac_mean=("p90_frac", "mean"),
).reset_index()
hf_summ = hf.groupby(["N", "depth"]).agg(
    hapfire_RMSE=("hapfire_af_RMSE", "mean"),
    hapfire_R2=("hapfire_af_R2", "mean"),
).reset_index()
final = summ.merge(hf_summ, on=["N", "depth"], how="left").sort_values(["depth", "N"])
sm_path = f"{OUT}/kmate_af_missingness_summary.tsv"
final.to_csv(sm_path, sep="\t", index=False)
print(f"wrote {sm_path}")
pd.set_option("display.width", 220, "display.max_columns", 30)
print("\n=== SUMMARY (mean over seeds) ===")
print(final.to_string(index=False))
