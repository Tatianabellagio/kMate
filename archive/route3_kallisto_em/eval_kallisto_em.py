"""Quick R² comparison helper for kallisto-EM vs baselines on the recomb sims.

Usage:
    python eval_kallisto_em.py --sim cov50_n50_g3_s42_hotspots_p231_chr1 \\
        --kal-tsv /tmp/kallisto_em_full.tsv --label kallisto_full

Compares against ground truth and the existing window_10kb / window_200kb
baselines for the same sim. Prints R² and MAE per (var_type, ld_tier) tile.
"""
from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np
import pandas as pd

RECOMB_ROOT = Path("/carnegie/nobackup/scratch/tbellagio/hapfire_sv/sims/visor_freqk/pool_sweep_82_recomb")


def load_truth(sim: str) -> pd.DataFrame:
    return pd.read_csv(RECOMB_ROOT / sim / "recomb_truth.tsv.gz", sep="\t")


def load_method(sim: str, name: str) -> pd.DataFrame:
    return pd.read_csv(RECOMB_ROOT / sim / f"{name}.tsv", sep="\t")


def stats(df: pd.DataFrame, col: str) -> tuple[int, float, float]:
    t = df.truth_af.values
    p = df[col].values
    f = np.isfinite(t) & np.isfinite(p)
    if f.sum() < 5:
        return int(f.sum()), float("nan"), float("nan")
    r2 = float(np.corrcoef(t[f], p[f])[0, 1] ** 2)
    mae = float(np.mean(np.abs(p[f] - t[f])))
    return int(f.sum()), r2, mae


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sim", required=True,
                    help="e.g. cov50_n50_g3_s42_hotspots_p231_chr1")
    ap.add_argument("--kal-tsv", required=True)
    ap.add_argument("--label", default="kallisto",
                    help="column label for the kal-tsv (default kallisto)")
    ap.add_argument("--baselines", nargs="*",
                    default=["cactus_em_recomb_window_10kb",
                            "cactus_em_recomb_window_100kb",
                            "cactus_em_recomb_window"])
    args = ap.parse_args()

    truth = load_truth(args.sim)
    df = truth.copy()

    kal = pd.read_csv(args.kal_tsv, sep="\t")
    df[args.label] = kal.alt_freq.values

    methods = [args.label]
    for b in args.baselines:
        try:
            df[b] = load_method(args.sim, b).alt_freq.values
            methods.append(b)
        except Exception as e:
            print(f"  WARN: couldn't load {b}: {e}", flush=True)

    poly = df[(df.truth_af > 0) & (df.truth_af < 1)
              & (df.chrom == "Chr1")].copy()
    poly["indel"] = (poly.alt_len - poly.ref_len).abs()
    poly["vt"] = np.where(
        (poly.ref_len == 1) & (poly.alt_len == 1), "SNP",
        np.where(poly.indel >= 50, "big_SV", "OTHER"))

    print(f"=== {args.sim} (Chr1 polymorphic) ===")
    for vt, name_vt in [("SNP", "SNP"), ("big_SV", "big_SV (≥50bp)")]:
        sub = poly[poly.vt == vt]
        print(f"\n-- {name_vt} (n={len(sub):,}) --")
        print(f"  {'method':<35s} {'n':>10s} {'R²':>8s} {'MAE':>8s}")
        for m in methods:
            n, r2, mae = stats(sub, m)
            print(f"  {m:<35s} {n:>10,} {r2:>8.3f} {mae:>8.3f}")


if __name__ == "__main__":
    main()
