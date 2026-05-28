"""
Calibrate per-record alt_freq predictions from per_sample_driver output.

Background:
    The 231-founder cn (with 151 imputed founders) produces per-record
    alt_freq predictions that are systematically over-predicted by a
    constant slope ≈ 1.42, with near-zero intercept. The bias comes from
    rank-deficiency at the founder level: EM puts ~1.4× too much mass on
    imputed founders relative to recipe truth, but the SHAPE of the
    per-record alt_freq distribution is preserved (Pearson r > 0.99).

    Validated on SEEDMIX_S1, S2, S3 (231-founder cn): slope = 1.418,
    1.414, 1.423; raw R² = 0.74-0.76; calibrated R² = 0.989-0.991.

Usage:
    # 1. Compute calibration from a known-truth pool (e.g. SEEDMIX vs recipe):
    python calibrate_alt_freqs.py compute \\
        --predicted-tsv results/seedmix_231/SEEDMIX_S1.tsv \\
        --recipe ../data/seedmix_recipe_normalized.tsv \\
        --cn-var data/cn_var_231.cn_var.npz \\
        --cn-var-meta data/cn_var_231.meta.npz \\
        --out calibration_seedmix_s1.json

    # 2. Apply to evolved sample TSVs:
    python calibrate_alt_freqs.py apply \\
        --calibration calibration_seedmix_s1.json \\
        --in-tsv results/per_sample_231/MLFH010120180409.tsv \\
        --out-tsv results/per_sample_231_calibrated/MLFH010120180409.tsv
"""
from __future__ import annotations
import argparse, json, os
import numpy as np, pandas as pd
from scipy.sparse import load_npz


def cmd_compute(args):
    cn_var = load_npz(args.cn_var)
    meta = np.load(args.cn_var_meta, allow_pickle=True)
    founders = list(meta["founders"])
    recipe = pd.read_csv(args.recipe, sep="\t")
    rd = dict(zip(recipe.ID.astype(str), recipe.seed_prop))
    h_truth = np.array([rd.get(str(f), 0.0) for f in founders])
    h_truth = h_truth / h_truth.sum()
    af_truth = (h_truth.astype(np.float32) @ cn_var.toarray()).astype(np.float64)

    df = pd.read_csv(args.predicted_tsv, sep="\t")
    if len(df) != len(af_truth):
        raise SystemExit(f"len(predicted)={len(df)} ≠ len(truth)={len(af_truth)}")
    af_pred = df["alt_freq"].to_numpy()

    # Linear regression
    slope, intercept = np.polyfit(af_truth, af_pred, 1)
    # Inverse: predicted = slope*truth + intercept → calibrated = (predicted - intercept) / slope
    r = float(np.corrcoef(af_pred, af_truth)[0, 1])
    rmse_raw = float(np.sqrt(np.mean((af_pred - af_truth)**2)))
    af_calib = (af_pred - intercept) / slope
    rmse_cal = float(np.sqrt(np.mean((af_calib - af_truth)**2)))
    r2_raw = 1 - ((af_pred - af_truth)**2).sum() / ((af_truth - af_truth.mean())**2).sum()
    r2_cal = 1 - ((af_calib - af_truth)**2).sum() / ((af_truth - af_truth.mean())**2).sum()

    print(f"  slope:     {slope:.4f}")
    print(f"  intercept: {intercept:.6f}")
    print(f"  Pearson r: {r:.4f}")
    print(f"  RMSE (raw):        {rmse_raw:.4f}  R²(raw): {r2_raw:.4f}")
    print(f"  RMSE (calibrated): {rmse_cal:.4f}  R²(calibrated): {r2_cal:.4f}")

    out = {
        "slope": float(slope),
        "intercept": float(intercept),
        "source_predicted": os.path.abspath(args.predicted_tsv),
        "source_recipe": os.path.abspath(args.recipe),
        "source_cn_var": os.path.abspath(args.cn_var),
        "n_records": int(len(af_truth)),
        "raw_R2": float(r2_raw),
        "calibrated_R2": float(r2_cal),
        "raw_RMSE": rmse_raw,
        "calibrated_RMSE": rmse_cal,
        "pearson_r": r,
    }
    with open(args.out, "w") as fh:
        json.dump(out, fh, indent=2)
    print(f"\nWrote {args.out}")


def cmd_apply(args):
    with open(args.calibration) as fh:
        cal = json.load(fh)
    slope = cal["slope"]
    intercept = cal["intercept"]
    print(f"  applying slope={slope:.4f}, intercept={intercept:.6f}")
    df = pd.read_csv(args.in_tsv, sep="\t")
    df["alt_freq_calibrated"] = ((df["alt_freq"] - intercept) / slope).clip(0, 1)
    df["alt_freq_raw"] = df["alt_freq"]
    df["alt_freq"] = df["alt_freq_calibrated"]   # canonical column = calibrated
    os.makedirs(os.path.dirname(args.out_tsv) or ".", exist_ok=True)
    df.to_csv(args.out_tsv, sep="\t", index=False)
    print(f"  wrote {args.out_tsv}  ({len(df):,} records)")


def main():
    ap = argparse.ArgumentParser()
    sp = ap.add_subparsers(dest="cmd", required=True)

    a = sp.add_parser("compute", help="Compute calibration from a known-truth pool")
    a.add_argument("--predicted-tsv", required=True)
    a.add_argument("--recipe", required=True)
    a.add_argument("--cn-var", required=True)
    a.add_argument("--cn-var-meta", required=True)
    a.add_argument("--out", required=True)
    a.set_defaults(func=cmd_compute)

    b = sp.add_parser("apply", help="Apply calibration to a per-sample TSV")
    b.add_argument("--calibration", required=True)
    b.add_argument("--in-tsv", required=True)
    b.add_argument("--out-tsv", required=True)
    b.set_defaults(func=cmd_apply)

    args = ap.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
