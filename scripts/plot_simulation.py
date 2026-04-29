#!/usr/bin/env python3
"""
Plot the simulation truth vs hapFIRE projection vs freqk for various reps and
coverages at 1kb DEL.

Truth: pool of 231 ecotypes, the FIRST N_SV = round(231*f) carry the 1-kb DEL.
       So truth_af_alt = N_SV / 231 (for f in {0.10, 0.30, 0.50, 0.70, 0.90}).

For hapFIRE: projection = sum over first N_SV ecotypes (in VCF order) of f_eco_hapfire.

The key claim being tested: at LOW coverage (cov10), freqk struggles because
breakpoint k-mers aren't sampled enough, but hapFIRE-projection uses
genome-wide SNPs which are abundant even at 10x — so it should be ~unaffected.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pathlib import Path

ROOT = Path("/home/tbellagio/scratch/hapfire_sv")
N_SAMPLES = 231
SIZE = "1kb"
K = 31
FREQS = [10, 30, 50, 70, 90]
VCF_SAMPLES = ROOT / "data" / "vcf_samples_231.txt"

vcf_samples = [s.strip() for s in VCF_SAMPLES.read_text().splitlines() if s.strip()]
assert len(vcf_samples) == N_SAMPLES


def load_hapfire_ecotype(path: Path):
    df = pd.read_csv(path, sep="\t", header=None, names=["ecotype","freq"], dtype={"ecotype":str})
    return dict(zip(df["ecotype"], df["freq"]))


def load_freqk_af(path: Path):
    parts = path.read_text().strip().split("\n")[-1].strip().split("|")
    return float(parts[0]), float(parts[1])


def gather_one(rep: str, cov: int):
    rows = []
    for f in FREQS:
        n_sv = round(N_SAMPLES * (f / 100.0))
        sv_carriers = vcf_samples[:n_sv]
        truth = n_sv / N_SAMPLES

        hf_path = ROOT / "results" / f"{rep}_cov{cov}_{SIZE}_f{f}" / f"{rep}_cov{cov}_{SIZE}_f{f}_ecotype_frequency_selected.txt"
        hf_proj = np.nan
        hf_panel_sum = np.nan
        if hf_path.exists() and hf_path.stat().st_size > 100:
            ef = load_hapfire_ecotype(hf_path)
            hf_panel_sum = sum(ef.values())
            hf_proj = sum(ef.get(c, 0.0) for c in sv_carriers)

        freqk_path = (
            Path("/home/tbellagio/scratch/visor_freqk/results/del")
            / rep / "var" / f"cov{cov}_err001" / SIZE / f"n{N_SAMPLES}" / f"f{f}" / f"k{K}"
            / f"var_del_{SIZE}_n{N_SAMPLES}_f{f}_err001.allele_frequencies.k{K}.tsv"
        )
        af_freqk = np.nan
        if freqk_path.exists():
            _, af_freqk = load_freqk_af(freqk_path)

        rows.append({
            "rep": rep, "cov": cov, "sv_freq_pct": f, "n_sv": n_sv,
            "truth": truth,
            "hapfire_proj": hf_proj,
            "freqk": af_freqk,
            "hapfire_panel_sum": hf_panel_sum,
            "hapfire_done": not np.isnan(hf_proj),
        })
    return pd.DataFrame(rows)


def main():
    reps = ["rep9", "rep1", "rep10", "rep11"]
    covs = [10, 20, 50]
    all_df = []
    for r in reps:
        for c in covs:
            all_df.append(gather_one(r, c))
    df = pd.concat(all_df, ignore_index=True)
    out_tsv = ROOT / "results" / "simulation_summary.tsv"
    df.to_csv(out_tsv, sep="\t", index=False)
    print(df[df["hapfire_done"] | df["freqk"].notna()].to_string(index=False))
    print(f"\nwrote: {out_tsv}")

    finished = df[df["hapfire_done"]]
    if len(finished) == 0:
        print("No hapFIRE simulation results yet — skipping plot.")
        return

    finished = finished.copy()
    finished["abs_err_hapfire"] = (finished["hapfire_proj"] - finished["truth"]).abs()
    finished["abs_err_freqk"]   = (finished["freqk"]        - finished["truth"]).abs()

    # Figure 1: truth vs estimate scatter, colored by coverage, by-rep marker
    fig, axes = plt.subplots(1, len(covs), figsize=(5*len(covs), 5), sharex=True, sharey=True)
    if len(covs) == 1:
        axes = [axes]
    for ax, cov in zip(axes, covs):
        sub = finished[finished["cov"] == cov]
        if len(sub) == 0:
            ax.set_title(f"cov{cov}: no hapFIRE results yet")
            continue
        for rep in sub["rep"].unique():
            s = sub[sub["rep"] == rep]
            ax.plot(s["truth"], s["hapfire_proj"], "o-", ms=8, alpha=0.85, label=f"hapFIRE-proj {rep}")
            ax.plot(s["truth"], s["freqk"], "s--", ms=6, alpha=0.6, label=f"freqk {rep}")
        ax.plot([0,1],[0,1], "k--", lw=0.5)
        ax.set_xlim(-0.02, 1.02); ax.set_ylim(-0.02, 1.02)
        ax.set_xlabel("Truth SV ALT freq")
        ax.set_ylabel("Estimated freq")
        ax.set_title(f"cov{cov}× — rep9 1kb DEL on uniform 231-ecotype pool")
        ax.legend(loc="upper left", fontsize=7)
        ax.grid(True, alpha=0.3)
    fig.tight_layout()
    out = ROOT / "results" / "simulation_truth_vs_estimate_by_cov.png"
    fig.savefig(out, dpi=120)
    print(f"wrote: {out}")

    # Figure 2: mean abs error by coverage and freq (this is the key plot)
    grouped = finished.groupby(["cov", "sv_freq_pct"]).agg(
        n=("rep", "size"),
        mae_hapfire=("abs_err_hapfire", "mean"),
        mae_freqk=("abs_err_freqk", "mean"),
        std_hapfire=("abs_err_hapfire", "std"),
        std_freqk=("abs_err_freqk", "std"),
    ).reset_index()
    print("\nGrouped MAE by (cov, sv_freq_pct):")
    print(grouped.round(4).to_string(index=False))

    fig2, axes = plt.subplots(1, 2, figsize=(13, 5))
    # MAE by coverage (averaged over freqs)
    ax = axes[0]
    by_cov = grouped.groupby("cov").agg(mae_hapfire=("mae_hapfire", "mean"),
                                         mae_freqk=("mae_freqk", "mean")).reset_index()
    x = np.arange(len(by_cov)); w = 0.35
    ax.bar(x-w/2, by_cov["mae_hapfire"], w, label="hapFIRE-proj", color="#1f77b4")
    ax.bar(x+w/2, by_cov["mae_freqk"],   w, label="freqk",        color="#ff7f0e")
    ax.set_xticks(x); ax.set_xticklabels([f"cov{c}×" for c in by_cov["cov"]])
    ax.set_ylabel("Mean |estimate - truth|"); ax.set_title("Accuracy vs coverage (mean over 5 freqs × replicates)")
    ax.grid(True, axis="y", alpha=0.3); ax.legend()

    # MAE by freq, lines per coverage (showing low-cov is where the gap matters)
    ax = axes[1]
    for cov in sorted(grouped["cov"].unique()):
        sub = grouped[grouped["cov"] == cov]
        ax.plot(sub["sv_freq_pct"], sub["mae_hapfire"], "o-", label=f"hapFIRE-proj cov{cov}×")
        ax.plot(sub["sv_freq_pct"], sub["mae_freqk"], "s--", label=f"freqk cov{cov}×", alpha=0.7)
    ax.set_xlabel("SV ALT freq (%)"); ax.set_ylabel("Mean |estimate - truth|")
    ax.set_title("Accuracy vs SV freq, by coverage")
    ax.legend(fontsize=7); ax.grid(True, alpha=0.3)
    fig2.tight_layout()
    out2 = ROOT / "results" / "simulation_mae_by_cov.png"
    fig2.savefig(out2, dpi=120)
    print(f"wrote: {out2}")


if __name__ == "__main__":
    main()
