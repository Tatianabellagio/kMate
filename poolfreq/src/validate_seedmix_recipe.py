"""
Validate per_sample_driver output against the SEEDMIX recipe at the
per-VCF-record level.

For each biallelic record r in cn_var:
    expected_alt_freq[r] = sum_f (recipe[f] · cn_var[f, r])

where recipe is normalized so that sum_f recipe[f] = 1 over the panel founders.

For the 82-founder panel: only the 80 GrENE-overlap founders carry recipe
mass (the recipe restricted to the panel sums to ~35% of total — the
"missing-mass" problem). For the 231-founder panel: all founders carry
recipe mass and sum to 1.

The driver's output (one TSV per SEEDMIX replicate) is compared at the
record level: R² and RMSE of predicted alt_freq vs expected_alt_freq.

Usage:
    python validate_seedmix_recipe.py \\
        --cn-var data/cn_var_82.cn_var.npz \\
        --cn-var-meta data/cn_var_82.meta.npz \\
        --panel-map ../data/sv_panel_to_accession_id.tsv \\
        --recipe ../data/seedmix_recipe_normalized.tsv \\
        --driver-output results/smoke/SEEDMIX_S1.tsv \\
        --label "SEEDMIX_S1 (82-founder panel)"

For 231-founder panel where founder names ARE the recipe IDs, omit
--panel-map.
"""
from __future__ import annotations
import argparse
import numpy as np
import pandas as pd
from scipy.sparse import load_npz


def metrics(y_pred, y_true):
    y_pred = np.asarray(y_pred, dtype=np.float64)
    y_true = np.asarray(y_true, dtype=np.float64)
    rmse = float(np.sqrt(np.mean((y_pred - y_true) ** 2)))
    if y_true.std() > 1e-12:
        ss_res = ((y_true - y_pred) ** 2).sum()
        ss_tot = ((y_true - y_true.mean()) ** 2).sum()
        r2 = float(1 - ss_res / ss_tot)
    else:
        r2 = float("nan")
    # Also report Pearson r, since R² can be very negative when prediction is
    # systematically biased — Pearson tells us shape agreement separately
    if y_pred.std() > 1e-12 and y_true.std() > 1e-12:
        r = float(np.corrcoef(y_pred, y_true)[0, 1])
    else:
        r = float("nan")
    return r2, rmse, r


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cn-var", required=True)
    ap.add_argument("--cn-var-meta", required=True)
    ap.add_argument("--panel-map", default=None,
                    help="Assembly_ID → Accession_ID TSV. Required if cn_var "
                         "founder names are Assembly_IDs (82-founder panel). "
                         "Omit if founders are already 1001G IDs (231-founder).")
    ap.add_argument("--recipe", required=True,
                    help="TSV with columns ID, seed_prop")
    ap.add_argument("--driver-output", required=True,
                    help="per_sample_driver TSV with columns chrom, pos, "
                         "ref_len, alt_len, alt_freq")
    ap.add_argument("--label", default="seedmix")
    ap.add_argument("--out-tsv", default=None,
                    help="optional: write per-record (expected, predicted) "
                         "to this path for downstream plotting")
    args = ap.parse_args()

    print(f"=== {args.label} ===")
    print(f"\n[1] Loading cn_var")
    cn_var = load_npz(args.cn_var)
    meta = np.load(args.cn_var_meta, allow_pickle=True)
    founders = list(meta["founders"])
    F, N = cn_var.shape
    print(f"  cn_var shape: {cn_var.shape}, nnz={cn_var.nnz:,}")
    print(f"  first 3 founder names: {founders[:3]}")

    # 2) Map cn_var founder names → 1001G accession IDs (if panel_map given)
    if args.panel_map:
        pm = pd.read_csv(args.panel_map, sep="\t")
        asm_to_1001g = dict(zip(pm.Assembly_ID.astype(str),
                                pm.Accession_ID.astype(str)))
        founder_1001g = [asm_to_1001g.get(str(f), None) for f in founders]
    else:
        founder_1001g = [str(f) for f in founders]
    print(f"  founders mapped to 1001G IDs: "
          f"{sum(x is not None for x in founder_1001g)}/{F}")

    # 3) Load recipe → assign per-founder weight
    recipe = pd.read_csv(args.recipe, sep="\t")
    recipe_dict = dict(zip(recipe.ID.astype(str), recipe.seed_prop))
    h = np.array([recipe_dict.get(fid, 0.0) if fid else 0.0
                  for fid in founder_1001g])
    panel_mass = h.sum()
    print(f"\n[2] Recipe panel mass: {panel_mass:.4f} ({panel_mass*100:.1f}% of "
          f"the seedmix is on this panel)")
    if panel_mass > 0:
        h_norm = h / panel_mass
    else:
        raise RuntimeError("Recipe has no overlap with cn_var founders")
    print(f"  effective n founders (1/Σh²): {1/np.sum(h_norm**2):.1f}")
    top_idx = np.argsort(h_norm)[::-1][:5]
    print(f"  top 5 founders by recipe weight:")
    for i in top_idx:
        print(f"    {founder_1001g[i]:>10}  {h_norm[i]:.4f}")

    # 4) Project recipe → expected per-record alt_freq
    # Use h_norm (sums to 1 over panel). For 231-founder panel, panel_mass≈1
    # so h_norm ≈ h. For 82-founder panel, h_norm rescales recipe to the
    # in-panel mass — this is the right comparison for h estimated under
    # simplex normalization.
    print(f"\n[3] Projecting recipe → expected per-record alt_freq")
    # cn_var is sparse F × N. h is F-vector. Want N-vector.
    expected = (h_norm @ cn_var.toarray() if F < 500
                else h_norm @ cn_var)  # h @ sparse works directly
    if hasattr(expected, "toarray"):
        expected = np.asarray(expected).flatten()
    expected = np.asarray(expected).flatten()
    print(f"  expected alt_freq: min={expected.min():.4f}, "
          f"max={expected.max():.4f}, mean={expected.mean():.4f}, "
          f"median={float(np.median(expected)):.4f}")

    # 5) Load driver output and align to cn_var record order
    print(f"\n[4] Loading driver output: {args.driver_output}")
    df = pd.read_csv(args.driver_output, sep="\t")
    print(f"  records: {len(df):,}")
    if len(df) != N:
        raise RuntimeError(
            f"driver output has {len(df):,} records but cn_var has {N:,}. "
            f"They should match if both come from the same biallelic VCF.")
    predicted = df["alt_freq"].to_numpy()

    # 6) Metrics
    r2, rmse, r = metrics(predicted, expected)
    print(f"\n[5] Comparison (n={N:,} records)")
    print(f"  R²:        {r2:>8.4f}")
    print(f"  RMSE:      {rmse:>8.4f}")
    print(f"  Pearson r: {r:>8.4f}")

    # 7) Stratify by AC (informative records have AC > 1)
    ac = np.asarray(cn_var.sum(axis=0)).flatten()
    bins = [(1, 1, "AC=1 (singletons)"),
            (2, 4, "AC 2–4 (rare)"),
            (5, 10, "AC 5–10 (low-freq)"),
            (11, F, "AC>10 (common)")]
    print(f"\n  Stratified by panel allele count:")
    print(f"  {'bin':<25} {'n':>10} {'R²':>8} {'RMSE':>8} {'r':>8}")
    print(f"  " + "-" * 65)
    for lo, hi, label in bins:
        mask = (ac >= lo) & (ac <= hi)
        if mask.sum() < 100:
            continue
        r2_b, rmse_b, r_b = metrics(predicted[mask], expected[mask])
        print(f"  {label:<25} {mask.sum():>10,} {r2_b:>8.4f} {rmse_b:>8.4f} "
              f"{r_b:>8.4f}")

    # 8) Optional dump
    if args.out_tsv:
        out = pd.DataFrame({
            "chrom": meta["chrom"],
            "pos": meta["pos"],
            "ref_len": meta["ref_len"],
            "alt_len": meta["alt_len"],
            "ac_panel": ac,
            "expected_alt_freq": expected,
            "predicted_alt_freq": predicted,
        })
        out.to_csv(args.out_tsv, sep="\t", index=False)
        print(f"\n  wrote per-record TSV → {args.out_tsv}")


if __name__ == "__main__":
    main()
