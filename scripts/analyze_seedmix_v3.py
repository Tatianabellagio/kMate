#!/usr/bin/env python
"""SEEDMIX_S1 v3 sanity check.

Compares cactus_em v3 SNP alt_freq against hapFIRE per-SNP frequencies.
Headline: slope (should be ~1.0 without calibration), R², Pearson r, MAE, RMSE.
Also prints v2 numbers side-by-side if its TSV is present.

Caveat: hapFIRE is not absolute truth (it's another method on the same reads),
but it is the established cactus_em ↔ hapFIRE benchmark — v2 hit R²=0.996
slope=1.013 against this same source. v3's job here is to keep the R² and
get rid of v2's 1.43× scale bias.
"""
import argparse
import sys
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]


def load_cactus_em(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    df = df[(df["ref_len"] == 1) & (df["alt_len"] == 1)].copy()  # SNPs only
    df["chrom_n"] = df["chrom"].str.replace("Chr", "", regex=False)
    return df[["chrom_n", "pos", "alt_freq"]].rename(columns={"alt_freq": "pred"})


def load_hapfire(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", header=None, names=["chrom_n", "pos", "hf_af"], dtype={"chrom_n": str})
    return df


def stats(pred: np.ndarray, truth: np.ndarray) -> dict:
    pred = np.asarray(pred, dtype=np.float64)
    truth = np.asarray(truth, dtype=np.float64)
    r = np.corrcoef(pred, truth)[0, 1]
    slope, intercept = np.polyfit(truth, pred, 1)
    resid = pred - truth
    return dict(
        n=len(pred),
        r=r,
        r2=r * r,
        slope=slope,
        intercept=intercept,
        mae=np.mean(np.abs(resid)),
        rmse=np.sqrt(np.mean(resid * resid)),
        mean_pred=np.mean(pred),
        mean_truth=np.mean(truth),
    )


def report(label: str, s: dict) -> None:
    print(
        f"{label:<8s}  n={s['n']:>9,}  r={s['r']:.4f}  R²={s['r2']:.4f}  "
        f"slope={s['slope']:.4f}  intercept={s['intercept']:+.5f}  "
        f"MAE={s['mae']:.4f}  RMSE={s['rmse']:.4f}"
    )


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--v3", type=Path, default=ROOT / "results/seedmix_231_v3_perchrom/SEEDMIX_S1.tsv")
    ap.add_argument("--v2", type=Path, default=ROOT / "results/seedmix_231_v2_perchrom/SEEDMIX_S1.tsv")
    ap.add_argument("--hapfire", type=Path, default=Path("/global/scratch/projects/fc_moilab/projects/grenenet-phase1/frequency/hapFIRE_frequencies/seed_mix/s1_snp_frequency.txt"))
    ap.add_argument("--restrict-chrom", default="1", help="evaluate on this chrom only ('1' for Chr1, '' for all)")
    args = ap.parse_args()

    if not args.v3.exists():
        sys.exit(f"v3 TSV not found: {args.v3}")
    if not args.hapfire.exists():
        sys.exit(f"hapFIRE truth not found: {args.hapfire}")

    print(f"loading hapFIRE truth: {args.hapfire}")
    hf = load_hapfire(args.hapfire)
    if args.restrict_chrom:
        hf = hf[hf["chrom_n"] == args.restrict_chrom].copy()
    print(f"  hapFIRE SNPs: {len(hf):,}")

    print(f"loading cactus_em v3: {args.v3}")
    v3 = load_cactus_em(args.v3)
    if args.restrict_chrom:
        v3 = v3[v3["chrom_n"] == args.restrict_chrom].copy()
    j3 = v3.merge(hf, on=["chrom_n", "pos"], how="inner")
    print(f"  v3 ∩ hapFIRE: {len(j3):,}  (v3 SNPs = {len(v3):,})")
    s3 = stats(j3["pred"].to_numpy(), j3["hf_af"].to_numpy())

    s2 = None
    if args.v2.exists():
        print(f"loading cactus_em v2: {args.v2}")
        v2 = load_cactus_em(args.v2)
        if args.restrict_chrom:
            v2 = v2[v2["chrom_n"] == args.restrict_chrom].copy()
        j2 = v2.merge(hf, on=["chrom_n", "pos"], how="inner")
        print(f"  v2 ∩ hapFIRE: {len(j2):,}")
        s2 = stats(j2["pred"].to_numpy(), j2["hf_af"].to_numpy())

    print("")
    print("=== SEEDMIX_S1 vs hapFIRE (SNPs) ===")
    report("v3", s3)
    if s2 is not None:
        report("v2", s2)
    print("")
    print("Targets:")
    print("  slope ≈ 1.000     (v2 was ~1.013 raw, ~1.43× w/ different cn_var build)")
    print("  R²    ≥ 0.99      (v2 hit 0.996)")
    print("  R² delta v3-v2:   {:+.4f}".format(s3["r2"] - s2["r2"]) if s2 else "")
    print("  slope delta v3-v2:{:+.4f}".format(s3["slope"] - s2["slope"]) if s2 else "")


if __name__ == "__main__":
    main()
