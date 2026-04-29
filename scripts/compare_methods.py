#!/usr/bin/env python3
"""
Sweep over the rep9 simulation frequencies, project SV freqs from hapFIRE
ecotype outputs, and compare to truth and freqk.

Usage:
    python compare_methods.py
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path("/home/tbellagio/scratch/hapfire_sv")
VCF_SAMPLES = ROOT / "data" / "vcf_samples_231.txt"
N_SAMPLES = 231
REP = "rep9"
COV = 50
SIZE = "1kb"
K = 31
FREQS = [10, 30, 50, 70, 90]


def load_hapfire_ecotype_freq(path: Path) -> dict:
    df = pd.read_csv(path, sep="\t", header=None, names=["ecotype", "freq"], dtype={"ecotype": str})
    return dict(zip(df["ecotype"], df["freq"]))


def load_freqk_af(path: Path) -> tuple[float, float]:
    line = path.read_text().strip().split("\n")[-1].strip()
    parts = line.split("|")
    return float(parts[0]), float(parts[1])


def main():
    vcf_samples = [s.strip() for s in VCF_SAMPLES.read_text().splitlines() if s.strip()]
    assert len(vcf_samples) == N_SAMPLES

    rows = []
    for f in FREQS:
        n_sv = round(N_SAMPLES * (f / 100.0))
        sv_carriers = set(vcf_samples[:n_sv])
        truth = n_sv / N_SAMPLES

        # hapFIRE result
        hf_path = ROOT / "results" / f"{REP}_cov{COV}_{SIZE}_f{f}" / f"{REP}_cov{COV}_{SIZE}_f{f}_ecotype_frequency_selected.txt"
        hf_proj = np.nan
        hf_panel_sum = np.nan
        if hf_path.exists():
            ef = load_hapfire_ecotype_freq(hf_path)
            hf_panel_sum = sum(ef.values())
            hf_proj = sum(ef.get(c, 0.0) for c in sv_carriers)

        # freqk result
        freqk_path = (
            Path("/home/tbellagio/scratch/visor_freqk/results/del")
            / REP / "var" / f"cov{COV}_err001" / SIZE / f"n{N_SAMPLES}" / f"f{f}" / f"k{K}"
            / f"var_del_{SIZE}_n{N_SAMPLES}_f{f}_err001.allele_frequencies.k{K}.tsv"
        )
        af_freqk = np.nan
        if freqk_path.exists():
            _, af_freqk = load_freqk_af(freqk_path)

        rows.append({
            "rep": REP,
            "cov": COV,
            "size": SIZE,
            "k": K,
            "sv_freq_pct": f,
            "n_sv": n_sv,
            "truth_af_alt": truth,
            "hapfire_panel_sum": hf_panel_sum,
            "af_alt_hapfire_proj": hf_proj,
            "af_alt_freqk": af_freqk,
            "err_hapfire": hf_proj - truth if not np.isnan(hf_proj) else np.nan,
            "err_freqk": af_freqk - truth if not np.isnan(af_freqk) else np.nan,
            "abs_err_hapfire": abs(hf_proj - truth) if not np.isnan(hf_proj) else np.nan,
            "abs_err_freqk": abs(af_freqk - truth) if not np.isnan(af_freqk) else np.nan,
            "hapfire_done": hf_path.exists(),
        })

    df = pd.DataFrame(rows)
    print(df.to_string(index=False))
    out = ROOT / "results" / f"compare_{REP}_cov{COV}_{SIZE}.tsv"
    df.to_csv(out, sep="\t", index=False)
    print(f"\nwrote: {out}")

    # Summary
    finished = df[df["hapfire_done"]]
    if len(finished):
        print(f"\nFinished {len(finished)}/{len(df)} hapFIRE runs")
        print(f"  hapFIRE mean abs err: {finished['abs_err_hapfire'].mean():.4f}")
        print(f"  freqk   mean abs err: {finished['abs_err_freqk'].mean():.4f}")


if __name__ == "__main__":
    main()
