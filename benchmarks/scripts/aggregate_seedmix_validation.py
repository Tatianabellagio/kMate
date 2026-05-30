"""
Aggregate validate_seedmix_recipe metrics across 8 SEEDMIX replicates and
optionally across two panels (82-founder vs 231-founder), producing a single
summary TSV.

For each (panel, replicate, ac_bin) cell we report n, R², RMSE, Pearson r.

Usage:
    python aggregate_seedmix_validation.py \\
        --label-82 "82-founder" \\
        --dir-82 results/seedmix_82 \\
        --var-pa-82 data/var_pa_82.var_pa.npz \\
        --meta-82 data/var_pa_82.meta.npz \\
        --panel-map ../data/sv_panel_to_accession_id.tsv \\
        --label-231 "231-founder" \\
        --dir-231 results/seedmix_231 \\
        --var-pa-231 data/var_pa_231.var_pa.npz \\
        --meta-231 data/var_pa_231.meta.npz \\
        --recipe ../data/seedmix_recipe_normalized.tsv \\
        --out results/seedmix_aggregate.tsv

The 231-founder panel uses 1001G IDs as founder names directly (no
panel-map needed).
"""
from __future__ import annotations
import argparse, glob, os
import numpy as np
import pandas as pd
from scipy.sparse import load_npz


AC_BINS = [
    ("AC=1",     1, 1),
    ("AC 2-4",   2, 4),
    ("AC 5-10",  5, 10),
    ("AC>10",   11, 999_999),
    ("ALL",      0, 999_999),
]


def metrics(y_pred, y_true):
    y_pred = np.asarray(y_pred, dtype=np.float64)
    y_true = np.asarray(y_true, dtype=np.float64)
    rmse = float(np.sqrt(np.mean((y_pred - y_true) ** 2)))
    r2 = float("nan")
    if y_true.std() > 1e-12:
        ss_res = ((y_true - y_pred) ** 2).sum()
        ss_tot = ((y_true - y_true.mean()) ** 2).sum()
        r2 = float(1 - ss_res / ss_tot)
    r = float("nan")
    if y_pred.std() > 1e-12 and y_true.std() > 1e-12:
        r = float(np.corrcoef(y_pred, y_true)[0, 1])
    return r2, rmse, r


def project_recipe(var_pa, founders, panel_map_path, recipe_path):
    """Return (expected_alt_freq vector, panel_total_mass, h_norm)."""
    if panel_map_path:
        pm = pd.read_csv(panel_map_path, sep="\t")
        asm_to_1001g = dict(zip(pm.Assembly_ID.astype(str),
                                pm.Accession_ID.astype(str)))
        founder_1001g = [asm_to_1001g.get(str(f), None) for f in founders]
    else:
        founder_1001g = [str(f) for f in founders]

    recipe = pd.read_csv(recipe_path, sep="\t")
    rd = dict(zip(recipe.ID.astype(str), recipe.seed_prop))
    h = np.array([rd.get(fid, 0.0) if fid else 0.0 for fid in founder_1001g])
    panel_total = float(h.sum())
    h_norm = h / panel_total if panel_total > 0 else h
    expected = np.asarray((h_norm @ var_pa)).flatten()
    return expected, panel_total, h_norm


def evaluate_panel(label, output_dir, var_pa_path, meta_path,
                   panel_map_path, recipe_path):
    print(f"\n=== {label} ===")
    print(f"  output dir: {output_dir}")
    var_pa = load_npz(var_pa_path)
    meta = np.load(meta_path, allow_pickle=True)
    founders = list(meta["founders"])
    F, N = var_pa.shape
    print(f"  var_pa: {var_pa.shape}, nnz={var_pa.nnz:,}")

    ac = np.asarray(var_pa.sum(axis=0)).flatten()
    expected, panel_total, h_norm = project_recipe(
        var_pa, founders, panel_map_path, recipe_path)
    print(f"  recipe panel mass: {panel_total*100:.1f}%, "
          f"effective n founders (1/Σh²): {1/np.sum(h_norm**2):.1f}")

    rows = []
    tsvs = sorted(glob.glob(os.path.join(output_dir, "*.tsv")))
    print(f"  found {len(tsvs)} per-sample TSVs")
    for path in tsvs:
        sid = os.path.splitext(os.path.basename(path))[0]
        df = pd.read_csv(path, sep="\t")
        if len(df) != N:
            print(f"  SKIP {sid}: {len(df):,} records ≠ {N:,} var_pa records")
            continue
        pred = df["alt_freq"].to_numpy()
        for bin_label, lo, hi in AC_BINS:
            mask = (ac >= lo) & (ac <= hi)
            if mask.sum() < 100:
                continue
            r2, rmse, r = metrics(pred[mask], expected[mask])
            rows.append({
                "panel": label, "sample": sid, "ac_bin": bin_label,
                "n": int(mask.sum()), "r2": r2, "rmse": rmse, "pearson_r": r,
            })
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--label-82", default="82-founder")
    ap.add_argument("--dir-82", required=True)
    ap.add_argument("--var-pa-82", required=True)
    ap.add_argument("--meta-82", required=True)
    ap.add_argument("--panel-map", default=None,
                    help="Assembly_ID → Accession_ID mapping (only needed if "
                         "var_pa founder names are Assembly IDs, e.g. for the "
                         "82-founder var_pa). Omit for 231-founder var_pa "
                         "where founders are already 1001G IDs.")
    ap.add_argument("--label-231", default="231-founder")
    ap.add_argument("--dir-231", default=None,
                    help="omit if 231 results aren't ready yet — only 82 reported")
    ap.add_argument("--var-pa-231", default=None)
    ap.add_argument("--meta-231", default=None)
    ap.add_argument("--recipe", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    rows = evaluate_panel(
        args.label_82, args.dir_82, args.var_pa_82, args.meta_82,
        args.panel_map, args.recipe,
    )
    if args.dir_231 and os.path.isdir(args.dir_231):
        rows += evaluate_panel(
            args.label_231, args.dir_231, args.var_pa_231, args.meta_231,
            None, args.recipe,
        )

    if not rows:
        raise RuntimeError("No metrics computed — no per-sample TSVs found?")

    out = pd.DataFrame(rows)
    out.to_csv(args.out, sep="\t", index=False)

    # Pretty print: pivot per (panel, ac_bin) summary
    print(f"\n=== Summary ===")
    summary = out.groupby(["panel", "ac_bin"]).agg(
        n_samples=("sample", "nunique"),
        r2_mean=("r2", "mean"),
        r2_min=("r2", "min"),
        r2_max=("r2", "max"),
        rmse_mean=("rmse", "mean"),
        pearson_r_mean=("pearson_r", "mean"),
    ).reset_index()
    print(summary.to_string(index=False))
    print(f"\nfull table → {args.out}")


if __name__ == "__main__":
    main()
