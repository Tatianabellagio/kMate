"""Append a kallisto-EM run's per-record AFs to the session_summary parquets so
the RECOMB_SWEEP_RESULTS.ipynb plots render the new method side-by-side with
the existing methods.

Specifically:
- Load /tmp/<run>.tsv (chrom, pos, ref_len, alt_len, alt_freq)
- Inner-join against af_long.parquet's rows for the target regime (which already
  carry the canonical filter set + var_type + ld_tier per record)
- Use af_long's first method (e.g. window_10kb) on that regime as the
  template for the record set; copy chrom/pos/ref_len/alt_len/truth_af/var_type/
  rec_idx/max_r2_xwu/sv_size from it.
- Build new rows with method=<method-name>, est_af=<kal alt_freq>, regime=<regime>
- Append to af_long.parquet (write back to disk)
- Also compute per-stratification R²/MAE/RMSE/slope/intercept and append to
  af_summary.parquet
"""
from __future__ import annotations
import argparse
from pathlib import Path

import numpy as np
import pandas as pd

RESULTS_DIR = Path("/carnegie/nobackup/scratch/tbellagio/hapfire_sv/results/session_summary")


def linreg_stats(t: np.ndarray, p: np.ndarray) -> dict:
    """R², MAE, RMSE, slope, intercept on finite pairs."""
    f = np.isfinite(t) & np.isfinite(p)
    if f.sum() < 5:
        return dict(n=int(f.sum()), R2=np.nan, MAE=np.nan, RMSE=np.nan,
                    slope=np.nan, intercept=np.nan)
    t = t[f]; p = p[f]
    r = float(np.corrcoef(t, p)[0, 1])
    r2 = r * r
    mae = float(np.mean(np.abs(p - t)))
    rmse = float(np.sqrt(np.mean((p - t) ** 2)))
    # slope, intercept of p = a*t + b
    if np.std(t) > 0:
        slope, intercept = np.polyfit(t, p, 1)
    else:
        slope, intercept = np.nan, np.nan
    return dict(n=int(f.sum()), R2=float(r2), MAE=float(mae), RMSE=float(rmse),
                slope=float(slope), intercept=float(intercept))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tsv", required=True,
                    help="kallisto-EM TSV with chrom/pos/ref_len/alt_len/alt_freq")
    ap.add_argument("--method-name", required=True,
                    help="method label to use in the parquet (e.g. kallisto_em)")
    ap.add_argument("--regime", required=True,
                    help="regime label (e.g. n50_g3)")
    ap.add_argument("--template-method", default="window_10kb",
                    help="existing method whose record set we mirror (default window_10kb)")
    ap.add_argument("--af-long",   default=str(RESULTS_DIR / "af_long.parquet"))
    ap.add_argument("--af-summary",default=str(RESULTS_DIR / "af_summary.parquet"))
    ap.add_argument("--dry-run", action="store_true",
                    help="don't write parquet, just print what would be appended")
    args = ap.parse_args()

    print(f"loading {args.af_long}...")
    af_long = pd.read_parquet(args.af_long)
    af_summary = pd.read_parquet(args.af_summary)

    # If method already exists for this regime, drop those rows first
    drop_long = (af_long.method == args.method_name) & (af_long.regime == args.regime)
    if drop_long.any():
        print(f"  removing {drop_long.sum():,} existing rows for "
              f"{args.method_name} @ {args.regime}")
        af_long = af_long[~drop_long].reset_index(drop=True)
    drop_sum = ((af_summary.method == args.method_name)
                & (af_summary.regime == args.regime))
    if drop_sum.any():
        af_summary = af_summary[~drop_sum].reset_index(drop=True)

    # Build the template subset
    tmpl = af_long[(af_long.method == args.template_method)
                   & (af_long.regime == args.regime)].copy()
    if len(tmpl) == 0:
        raise SystemExit(f"no template rows for "
                         f"method={args.template_method}, regime={args.regime}")
    print(f"template rows: {len(tmpl):,}")

    # Load kal TSV. The TSV has one row per cn_var record in the same order
    # as the cn_var meta file. rec_idx in af_long corresponds to position in
    # that meta. So we align positionally via rec_idx — NO key-merge (which
    # would inflate rows due to multi-allelic duplicate (chrom,pos,ref,alt)).
    print(f"loading {args.tsv}...")
    kal = pd.read_csv(args.tsv, sep="\t")
    print(f"  kal rows: {len(kal):,} (full cn_var record set)")
    rec_idx = tmpl.rec_idx.values
    if rec_idx.max() >= len(kal):
        raise SystemExit(f"rec_idx max {rec_idx.max()} >= kal rows {len(kal)} — "
                         f"kal TSV is shorter than expected")
    # Sanity: chrom/pos in template should match kal[rec_idx] chrom/pos
    sample_n = min(20, len(tmpl))
    sample_idxs = np.linspace(0, len(tmpl) - 1, sample_n).astype(int)
    mismatches = 0
    for i in sample_idxs:
        tr = tmpl.iloc[i]; kr = kal.iloc[int(tr.rec_idx)]
        if tr.chrom != kr.chrom or tr.pos != kr.pos:
            mismatches += 1
    if mismatches:
        raise SystemExit(f"{mismatches}/{sample_n} sanity-check mismatches — "
                         f"kal TSV ordering does NOT match cn_var meta order")
    print(f"  positional sanity check passed ({sample_n} samples)")

    new_rows = tmpl.copy()
    new_rows["est_af"] = kal["alt_freq"].values[rec_idx].astype(np.float32)
    new_rows["method"] = args.method_name
    new_rows = new_rows[af_long.columns]
    new_rows.method = new_rows.method.astype(af_long.method.dtype)

    # Compute summary stats per (var_type, ld_tier)
    sum_rows = []
    for vt in new_rows.var_type.unique():
        for tier in new_rows.ld_tier.unique():
            sub = new_rows[(new_rows.var_type == vt) & (new_rows.ld_tier == tier)]
            if len(sub) == 0:
                continue
            t = sub.truth_af.values.astype(float)
            p = sub.est_af.values.astype(float)
            s = linreg_stats(t, p)
            sum_rows.append(dict(regime=args.regime, method=args.method_name,
                                 var_type=vt, ld_tier=tier, **s))
        # Also add an "all" tier
        sub_all = new_rows[new_rows.var_type == vt]
        if len(sub_all) == 0:
            continue
        t = sub_all.truth_af.values.astype(float)
        p = sub_all.est_af.values.astype(float)
        s = linreg_stats(t, p)
        sum_rows.append(dict(regime=args.regime, method=args.method_name,
                             var_type=vt, ld_tier="all", **s))

    new_summary = pd.DataFrame(sum_rows)
    print(f"new summary rows: {len(new_summary)}")
    print("Summary preview:")
    show = new_summary[(new_summary.ld_tier == "all")
                       | (new_summary.var_type == "SNP")]
    print(show[["var_type", "ld_tier", "n", "R2", "MAE"]]
          .sort_values(["var_type", "ld_tier"]).to_string(index=False))

    if args.dry_run:
        print("\n[dry-run] not writing parquet")
        return

    af_long_new = pd.concat([af_long, new_rows], ignore_index=True)
    af_summary_new = pd.concat([af_summary, new_summary], ignore_index=True)
    print(f"writing {args.af_long} ({len(af_long_new):,} rows)...")
    af_long_new.to_parquet(args.af_long, index=False)
    print(f"writing {args.af_summary} ({len(af_summary_new):,} rows)...")
    af_summary_new.to_parquet(args.af_summary, index=False)
    print("done.")


if __name__ == "__main__":
    main()
