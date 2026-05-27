#!/usr/bin/env python
"""
Direct contamination measurement at non-231-private SNPs.

For each SEEDMIX sample:
  - Load the per-site AD output from 02_measure_af.sh
  - Join with diagnostic-site annotations (AF in non-231 panel, carrier list)
  - Compute observed AF = ad_alt / (ad_ref + ad_alt) per site
  - Apply quality filters (min DP, min AF in panel, etc.)
  - Estimate global contamination fraction:
        Under "uniform contamination from non-231 founders at fraction c",
        and at sites where panel-AF in non-231 is p, the expected seedmix-BAM
        AF is approximately c * p   (plus sequencing error e).
        So: c_hat = (mean obs_AF - e) / mean(p)
  - Estimate per-accession contamination via private-SNP analysis (sites
    where exactly ONE non-231 ecotype carries the alt; their AF in the
    seedmix BAM tells us about that specific accession's contribution).

Output:
  results/contamination_summary.tsv     per-sample headline numbers
  results/per_accession.tsv             per-non-231-accession contamination estimate

Reference values:
  - Sequencing error rate (Illumina): ~0.3-0.5%
  - 1/231 = 0.43%   (per-founder uniform recipe expectation)
"""

import argparse
import sys
from pathlib import Path
import numpy as np
import pandas as pd


ROOT = Path("/home/tbellagio/scratch/hapfire_sv/contamination_test/option1_af")
DEFAULT_DIAG_TSVS = sorted(ROOT.glob("sites/diag_chr*.tsv"))
DEFAULT_AF_TSVS = sorted(ROOT.glob("af/seedmix_S*.af.tsv"))
RESULTS = ROOT / "results"


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--min-dp", type=int, default=8,
                   help="Minimum DP per site to include (default 8 ~ 1x coverage at the site)")
    p.add_argument("--min-panel-af", type=float, default=0.05,
                   help="Minimum AF in non-231 panel per site (excludes very rare sites)")
    p.add_argument("--seq-error", type=float, default=0.005,
                   help="Per-base sequencing error rate to subtract (default 0.5%)")
    p.add_argument("--diag-tsvs", nargs="*", default=[str(p) for p in DEFAULT_DIAG_TSVS])
    p.add_argument("--af-tsvs", nargs="*", default=[str(p) for p in DEFAULT_AF_TSVS])
    p.add_argument("--out-prefix", default=str(RESULTS))
    return p.parse_args()


def load_diagnostic_sites(tsvs):
    dfs = [pd.read_csv(t, sep="\t", dtype={"chrom": str}) for t in tsvs]
    df = pd.concat(dfs, ignore_index=True)
    print(f"loaded {len(df):,} diagnostic sites across {len(tsvs)} chrom files", file=sys.stderr)
    return df


def load_af(tsvs):
    """Returns dict: sample -> df with chrom, pos, ad_ref, ad_alt, dp, af"""
    out = {}
    for t in tsvs:
        path = Path(t)
        # filename pattern: seedmix_S<N>.af.tsv
        sample = path.stem.replace(".af", "")
        df = pd.read_csv(t, sep="\t", dtype={"chrom": str})
        out[sample] = df
        print(f"  {sample}: {len(df):,} sites measured", file=sys.stderr)
    return out


def estimate_contamination_fraction(obs_af, panel_af, seq_error):
    """
    obs_af = c * panel_af + seq_error  (under uniform-contamination model)
    => c_hat = mean((obs_af - seq_error) / panel_af) over sites where panel_af > 0
    """
    per_site = (obs_af - seq_error) / panel_af
    return per_site.mean(), np.median(per_site)


def estimate_via_regression(obs_af, panel_af, dp):
    """
    Fit obs_af = c * panel_af + intercept (= seq error)
    DP-weighted least squares - sites with more reads carry more weight.
    Returns (c_hat, intercept_hat, r2)
    """
    # weights = depth (sites with more reads have lower variance)
    w = np.asarray(dp, dtype=float)
    x = np.asarray(panel_af, dtype=float)
    y = np.asarray(obs_af, dtype=float)
    # weighted least squares via numpy
    W = np.sqrt(w)
    A = np.column_stack([x * W, W])  # [c, intercept]
    b = y * W
    sol, *_ = np.linalg.lstsq(A, b, rcond=None)
    c_hat, intercept_hat = sol[0], sol[1]
    # r2 (un-weighted, on residuals)
    y_pred = c_hat * x + intercept_hat
    ss_res = np.sum((y - y_pred) ** 2)
    ss_tot = np.sum((y - y.mean()) ** 2)
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 0
    return c_hat, intercept_hat, r2


def main():
    args = parse_args()
    Path(args.out_prefix).mkdir(parents=True, exist_ok=True)

    sites = load_diagnostic_sites(args.diag_tsvs)
    sites["chrom_pos"] = sites["chrom"].astype(str) + ":" + sites["pos"].astype(str)
    sites = sites.set_index("chrom_pos")

    af_per_sample = load_af(args.af_tsvs)

    # Apply panel AF filter (exclude very rare sites where signal would be too low)
    n_before = len(sites)
    sites_used = sites[sites["af_non231"] >= args.min_panel_af].copy()
    print(f"\nKept {len(sites_used):,}/{n_before:,} sites with panel AF >= {args.min_panel_af}", file=sys.stderr)
    print(f"  panel AF range: [{sites_used['af_non231'].min():.3f}, {sites_used['af_non231'].max():.3f}]", file=sys.stderr)
    print(f"  panel AF median: {sites_used['af_non231'].median():.3f}", file=sys.stderr)

    summary = []
    for sample, af_df in af_per_sample.items():
        af_df = af_df.copy()
        af_df["chrom_pos"] = af_df["chrom"].astype(str) + ":" + af_df["pos"].astype(str)
        af_df = af_df.set_index("chrom_pos")

        # join with diagnostic annotations
        joined = af_df.join(sites_used[["af_non231", "ac_non231", "carriers"]], how="inner")

        # filter by min DP
        n_total = len(joined)
        joined = joined[joined["dp"] >= args.min_dp]
        n_kept = len(joined)

        if n_kept == 0:
            print(f"WARN: {sample} has 0 sites after DP filter", file=sys.stderr)
            continue

        # diagnostic AFs and depth stats
        mean_af = joined["af"].mean()
        median_af = joined["af"].median()
        mean_dp = joined["dp"].mean()
        mean_panel_af = joined["af_non231"].mean()

        # contamination fraction estimate (uniform-from-non-231 model)
        c_mean, c_median = estimate_contamination_fraction(
            joined["af"].values, joined["af_non231"].values, args.seq_error
        )

        # naive: just observed AF minus error
        c_naive = mean_af - args.seq_error

        # WLS regression of obs_af on panel_af; slope = c_hat, intercept = seq error
        c_reg, eps_reg, r2 = estimate_via_regression(
            joined["af"].values, joined["af_non231"].values, joined["dp"].values
        )

        # bin by panel AF to see if relationship is linear (sanity check)
        joined["af_bin"] = pd.cut(joined["af_non231"], bins=[0, 0.1, 0.2, 0.3, 0.5, 1.0])
        af_by_bin = joined.groupby("af_bin", observed=True)["af"].mean()

        summary.append({
            "sample": sample,
            "n_sites_total": n_total,
            "n_sites_used": n_kept,
            "mean_dp": mean_dp,
            "mean_obs_af": mean_af,
            "median_obs_af": median_af,
            "mean_panel_af": mean_panel_af,
            "c_naive_mean_af_minus_error": c_naive,
            "c_uniform_contam_mean": c_mean,
            "c_uniform_contam_median": c_median,
            "c_regression_slope": c_reg,
            "regression_intercept": eps_reg,
            "regression_r2": r2,
            "obs_af_by_panel_af_0.05_0.10": af_by_bin.get(pd.Interval(0.0, 0.1), np.nan),
            "obs_af_by_panel_af_0.10_0.20": af_by_bin.get(pd.Interval(0.1, 0.2), np.nan),
            "obs_af_by_panel_af_0.20_0.30": af_by_bin.get(pd.Interval(0.2, 0.3), np.nan),
            "obs_af_by_panel_af_0.30_0.50": af_by_bin.get(pd.Interval(0.3, 0.5), np.nan),
            "obs_af_by_panel_af_0.50_1.00": af_by_bin.get(pd.Interval(0.5, 1.0), np.nan),
        })

    summary = pd.DataFrame(summary)
    out = Path(args.out_prefix) / "contamination_summary.tsv"
    summary.to_csv(out, sep="\t", index=False, float_format="%.5f")
    print(f"\nwrote {out}", file=sys.stderr)
    print()
    print("=== CONTAMINATION SUMMARY ===")
    cols = ["sample", "n_sites_used", "mean_dp", "mean_obs_af",
            "mean_panel_af", "c_naive_mean_af_minus_error",
            "c_uniform_contam_mean", "c_regression_slope",
            "regression_intercept", "regression_r2"]
    print(summary[cols].to_string(index=False, float_format="%.4f"))

    print()
    print("=== AF VS PANEL-AF (sanity check: should be ~linear if uniform contamination) ===")
    bin_cols = [c for c in summary.columns if c.startswith("obs_af_by_panel_af")]
    print(summary[["sample"] + bin_cols].to_string(index=False, float_format="%.4f"))

    # Per-accession analysis: sites where exactly ONE non-231 carrier
    print()
    print("=== PER-ACCESSION CONTAMINATION (private SNPs only) ===")
    private_sites = sites[sites["ac_non231"] == 1].copy()
    print(f"private sites (AC=1 in non-231): {len(private_sites):,}", file=sys.stderr)
    if len(private_sites) > 0:
        private_sites["chrom_pos"] = private_sites.index
        private_sites["accession"] = private_sites["carriers"]

        rows = []
        for sample, af_df in af_per_sample.items():
            af_df = af_df.copy()
            af_df["chrom_pos"] = af_df["chrom"].astype(str) + ":" + af_df["pos"].astype(str)
            af_df = af_df.set_index("chrom_pos")
            joined = af_df.join(private_sites[["accession", "af_non231"]], how="inner")
            joined = joined[joined["dp"] >= args.min_dp]
            # for AC=1 sites in panel, panel_af = 1/(2*910) (1 carrier homozygous in inbred Arabidopsis)
            # so under uniform contamination at fraction c: expected obs_af = c * 1/910
            # but if contamination is FROM that specific accession only, expected obs_af = c_specific * 1.0
            # so obs_af at AC=1 sites tells us about specific-accession contribution
            grouped = joined.groupby("accession").agg(
                n_sites=("af", "count"),
                mean_obs_af=("af", "mean"),
                median_obs_af=("af", "median"),
                mean_dp=("dp", "mean"),
            ).reset_index()
            grouped["sample"] = sample
            rows.append(grouped)

        all_per_acc = pd.concat(rows, ignore_index=True)
        # aggregate across samples
        per_acc = all_per_acc.groupby("accession").agg(
            n_sites=("n_sites", "first"),
            mean_obs_af=("mean_obs_af", "mean"),
            max_obs_af=("mean_obs_af", "max"),
            mean_dp=("mean_dp", "mean"),
        ).reset_index().sort_values("mean_obs_af", ascending=False)

        out2 = Path(args.out_prefix) / "per_accession.tsv"
        per_acc.to_csv(out2, sep="\t", index=False, float_format="%.5f")
        print(f"wrote {out2}", file=sys.stderr)

        # filter to accessions with sufficient evidence (≥10 private sites, mean DP ≥ 5)
        well_supported = per_acc[(per_acc["n_sites"] >= 10) & (per_acc["mean_dp"] >= 5)]
        print(f"\nTop 20 non-231 accessions by mean obs AF at their private SNPs (≥10 sites, mean DP ≥ 5):")
        print(well_supported.head(20).to_string(index=False, float_format="%.4f"))


if __name__ == "__main__":
    main()
