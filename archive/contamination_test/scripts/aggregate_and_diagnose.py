#!/usr/bin/env python
"""
Aggregate per-chrom hapFIRE outputs and compute the contamination diagnostic.

For each SEEDMIX_S{1..8} sample we have 5 per-chrom ecotype_frequency.txt files,
each summing to 1 over the 1141-panel. We average across chroms (simple mean,
matching hapFIRE.py's own across-chrom averaging) to get one 1141-vector per
sample, then split it into:

  in231_mass     = sum over the 231 GrENE-Net founders     (expected ~ 1.0)
  non231_mass    = sum over the remaining 910 ecotypes     (contamination signal)
  per-ecotype freq  = the full 1141-vector for inspection

Outputs:
  results/per_sample_ecotype_freq_1141.tsv   sample x ecotype freq matrix
  results/contamination_summary.tsv          per-sample headline numbers
  results/non231_top.tsv                     top non-231 ecotypes ranked by mean mass
"""

import argparse
from pathlib import Path
import sys

import pandas as pd


_ROOT = Path(__file__).resolve().parents[2]


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--results-dir", default=str(_ROOT / "contamination_test/results"))
    p.add_argument("--grene231-list", default=str(_ROOT / "data/vcf_samples_231.txt"))
    p.add_argument("--out-prefix", default=str(_ROOT / "contamination_test/results"))
    p.add_argument("--samples", nargs="+", default=[f"SEEDMIX_S{i}" for i in range(1, 9)])
    p.add_argument("--chroms", nargs="+", default=["1", "2", "3", "4", "5"])
    return p.parse_args()


def read_ecotype_freq(path):
    df = pd.read_csv(path, sep="\t", header=None, names=["ecotype", "freq"], dtype={"ecotype": str})
    return df.set_index("ecotype")["freq"]


def main():
    args = parse_args()
    results_dir = Path(args.results_dir)

    grene231 = set(Path(args.grene231_list).read_text().split())
    print(f"loaded {len(grene231)} GrENE-Net 231 IDs", file=sys.stderr)

    per_sample_avg = {}
    missing = []
    for sample in args.samples:
        chrom_vecs = []
        for ch in args.chroms:
            d = results_dir / f"{sample}_chr{ch}"
            # default hapFIRE output is "*_ecotype_frequency_selected.txt"
            # (diversity-weighted top-20% blocks); fall back to non-_selected.
            for name in (
                f"{sample}_chr{ch}_ecotype_frequency_selected.txt",
                f"{sample}_chr{ch}_ecotype_frequency.txt",
            ):
                f = d / name
                if f.exists():
                    break
            else:
                missing.append(str(d / f"{sample}_chr{ch}_ecotype_frequency*.txt"))
                continue
            v = read_ecotype_freq(f)
            chrom_vecs.append(v)
        if not chrom_vecs:
            print(f"WARN: no chrom outputs for {sample}", file=sys.stderr)
            continue
        # simple mean across chroms (mirrors hapFIRE.py's internal averaging)
        df = pd.concat(chrom_vecs, axis=1).fillna(0.0)
        per_sample_avg[sample] = df.mean(axis=1)
        print(f"{sample}: averaged {len(chrom_vecs)} chroms, {len(per_sample_avg[sample])} ecotypes", file=sys.stderr)

    if missing:
        print(f"\nMISSING ({len(missing)}):", file=sys.stderr)
        for m in missing:
            print("  " + m, file=sys.stderr)

    if not per_sample_avg:
        sys.exit("no samples completed")

    M = pd.DataFrame(per_sample_avg)  # rows=ecotype, cols=sample
    M.index.name = "ecotype"

    # write the full matrix
    out_full = Path(args.out_prefix) / "per_sample_ecotype_freq_1141.tsv"
    M.to_csv(out_full, sep="\t")
    print(f"wrote {out_full}", file=sys.stderr)

    # split mass and write contamination summary
    in231_idx = M.index.isin(grene231)
    summary = pd.DataFrame({
        "sum_total":      M.sum(axis=0),
        "in231_mass":     M[in231_idx].sum(axis=0),
        "non231_mass":    M[~in231_idx].sum(axis=0),
        "in231_n_called": (M[in231_idx] > 1e-4).sum(axis=0),
        "non231_n_called": (M[~in231_idx] > 1e-4).sum(axis=0),
    })
    out_sum = Path(args.out_prefix) / "contamination_summary.tsv"
    summary.to_csv(out_sum, sep="\t")
    print(f"wrote {out_sum}", file=sys.stderr)
    print("\nContamination summary:")
    print(summary.to_string(float_format="%.4f"))

    # rank non-231 ecotypes by mean mass across samples
    non231 = M[~in231_idx].copy()
    non231["mean_mass"] = non231.mean(axis=1)
    non231["max_mass"]  = non231.max(axis=1)
    non231 = non231.sort_values("mean_mass", ascending=False)
    out_top = Path(args.out_prefix) / "non231_top.tsv"
    non231.to_csv(out_top, sep="\t")
    print(f"wrote {out_top}", file=sys.stderr)
    print(f"\nTop 20 non-231 ecotypes by mean mass:")
    print(non231.head(20)[["mean_mass", "max_mass"] + list(M.columns)].to_string(float_format="%.4f"))


if __name__ == "__main__":
    main()
