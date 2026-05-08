#!/usr/bin/env python3
"""
Posterior-combine method (option A from the LD-empowered k-mer discussion).

Given:
  - hapFIRE-projection point estimate `f_proj` with sd `sd_proj`
  - freqk point estimate `f_freqk` with sd `sd_freqk`

Combine as inverse-variance-weighted average:
    f_combined = (f_proj/var_proj + f_freqk/var_freqk) / (1/var_proj + 1/var_freqk)
    var_combined = 1 / (1/var_proj + 1/var_freqk)

Variance models (rough first-pass):
  sd_proj   ≈ sqrt(N_carriers) · σ_eco          # σ_eco ~ per-ecotype hapFIRE sd
  sd_freqk  ≈ sqrt(p(1-p) / n_kmer_obs) + α·repeat_score
                                                 # binomial SE + repeat-region bias
For SVs not in the panel, set var_proj = ∞ → combined ≡ freqk.

This is a demo on the simulation data — both for reps where hapFIRE & freqk
both worked, and for the ones where freqk failed.
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path("/home/tbellagio/scratch/hapfire_sv")
SIGMA_ECO = 1e-4   # per-ecotype hapFIRE sd, calibrated from rep9 cov50 f30 results
ALPHA_REP = 0.05   # repeat-bias coefficient for freqk SD


def main():
    df = pd.read_csv(ROOT / "results" / "simulation_summary.tsv", sep="\t")
    df = df[df["hapfire_done"] & df["freqk"].notna()].copy()
    df["abs_err_hapfire"] = (df["hapfire_proj"] - df["truth"]).abs()
    df["abs_err_freqk"]   = (df["freqk"]        - df["truth"]).abs()

    # rough variance models
    df["sd_hapfire"] = np.sqrt(df["n_sv"]) * SIGMA_ECO
    n_kmer = 100  # placeholder; in practice would be n_unique_kmers from freqk index
    df["sd_freqk"] = np.sqrt(df["freqk"] * (1 - df["freqk"]) / n_kmer) + 0.0  # could add ALPHA*repeat_score
    # at low coverage the SE is bigger because n_kmer_observed scales w/ coverage
    df["sd_freqk"] = df["sd_freqk"] * np.sqrt(50.0 / df["cov"])

    # combine
    var_h = df["sd_hapfire"]**2
    var_f = df["sd_freqk"]**2
    df["w_hapfire"] = (1/var_h) / (1/var_h + 1/var_f)
    df["w_freqk"]   = (1/var_f) / (1/var_h + 1/var_f)
    df["combined"]  = df["w_hapfire"]*df["hapfire_proj"] + df["w_freqk"]*df["freqk"]
    df["abs_err_combined"] = (df["combined"] - df["truth"]).abs()

    print(df[["rep","cov","sv_freq_pct","truth","hapfire_proj","freqk","combined",
              "abs_err_hapfire","abs_err_freqk","abs_err_combined",
              "w_hapfire"]].round(4).to_string(index=False))
    print()
    agg = df.groupby("cov").agg(
        n=("rep","size"),
        mae_hapfire=("abs_err_hapfire","mean"),
        mae_freqk=("abs_err_freqk","mean"),
        mae_combined=("abs_err_combined","mean"),
    ).round(4)
    print("MAE by coverage:")
    print(agg.to_string())

    df.to_csv(ROOT / "results" / "simulation_combined.tsv", sep="\t", index=False)
    print(f"\nwrote: {ROOT / 'results' / 'simulation_combined.tsv'}")


if __name__ == "__main__":
    main()
