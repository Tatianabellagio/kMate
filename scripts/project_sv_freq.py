#!/usr/bin/env python3
"""
Project hapFIRE ecotype frequencies → SV allele frequency.

The math (exact mirror of how hapFIRE computes per-SNP freqs internally):
    f_SV  =  G_SV  @  f_ecotype
where
    G_SV       (n_ecotypes,) 0/1 indicator: which founders carry the SV ALT
    f_ecotype  (n_ecotypes,) hapFIRE's per-pool ecotype-frequency estimate
                (one row per ecotype, frequencies sum to ~1)

For visor_freqk simulation truth:
    The simulation pools the FIRST N_SAMPLES (=231) ecotypes from the GrENET VCF
    header at uniform 1/N each. The FIRST N_SV = round(N_SAMPLES * SV_FREQ) of
    those carry a single 1-kb deletion. Everyone else is reference.
    Truth:  f_ecotype_i = 1/N for all i in first N_SAMPLES, 0 otherwise.
            f_SV_alt    = N_SV / N_SAMPLES

Outputs a TSV with one row per (truth_freq, projection): truth, freqk_estimate,
hapFIRE_projection.
"""

import argparse
import sys
from pathlib import Path
import numpy as np
import pandas as pd


def load_hapfire_ecotype_freq(path: Path) -> pd.DataFrame:
    """hapFIRE _ecotype_frequency_selected.txt: 2-col TSV (ecotype_id, freq)."""
    df = pd.read_csv(path, sep="\t", header=None, names=["ecotype", "freq"], dtype={"ecotype": str})
    return df


def load_freqk_af(path: Path) -> tuple[float, float]:
    """freqk allele_frequencies.{k}.tsv: single line af_ref|af_alt."""
    line = Path(path).read_text().strip().split("\n")[-1].strip()
    parts = line.split("|")
    if len(parts) != 2:
        raise ValueError(f"unexpected freqk format in {path}: {line!r}")
    return float(parts[0]), float(parts[1])


def project(
    ecotype_freq_df: pd.DataFrame,
    sv_carrier_ecotypes: list[str],
) -> dict:
    """Compute SV ALT freq = sum of f_ecotype over carriers."""
    fmap = dict(zip(ecotype_freq_df["ecotype"], ecotype_freq_df["freq"]))
    n_panel = len(fmap)
    sum_panel = sum(fmap.values())
    carrier_set = set(sv_carrier_ecotypes)
    missing = carrier_set - set(fmap)
    if missing:
        print(
            f"  warning: {len(missing)}/{len(carrier_set)} SV-carrier ecotypes "
            f"not found in hapFIRE panel (e.g. {sorted(missing)[:3]})",
            file=sys.stderr,
        )
    af_alt_proj = sum(fmap[e] for e in carrier_set if e in fmap)
    return {
        "n_panel": n_panel,
        "sum_panel": sum_panel,
        "n_carriers": len(carrier_set),
        "n_carriers_in_panel": len(carrier_set & set(fmap)),
        "af_alt_proj": af_alt_proj,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ecotype-freq", required=True, type=Path,
                    help="hapFIRE _ecotype_frequency_selected.txt")
    ap.add_argument("--vcf-samples", required=True, type=Path,
                    help="bcftools query -l output (one ecotype per line, VCF header order)")
    ap.add_argument("--n-samples", required=True, type=int,
                    help="number of ecotypes pooled (e.g. 231)")
    ap.add_argument("--sv-freq", required=True, type=float,
                    help="SV ALT freq used in simulation (e.g. 0.5)")
    ap.add_argument("--freqk-af", type=Path, default=None,
                    help="optional: freqk allele_frequencies.tsv to compare against")
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()

    # 1. SV carriers: first N_SV samples in VCF header order
    vcf_samples = [
        s.strip() for s in args.vcf_samples.read_text().splitlines() if s.strip()
    ][: args.n_samples]
    n_sv = round(args.n_samples * args.sv_freq)
    sv_carriers = vcf_samples[:n_sv]
    truth_af_alt = n_sv / args.n_samples

    print(f"truth: N_SAMPLES={args.n_samples}, SV_FREQ={args.sv_freq}, "
          f"N_SV={n_sv}, truth AF_ALT={truth_af_alt:.4f}")

    # 2. hapFIRE ecotype frequencies
    ef = load_hapfire_ecotype_freq(args.ecotype_freq)
    print(f"hapFIRE ecotype panel: {len(ef)} entries; sum = {ef['freq'].sum():.4f}")
    print(f"  top-5 by freq: {ef.nlargest(5, 'freq').to_dict('records')}")

    # 3. Project
    proj = project(ef, sv_carriers)
    print(f"projection: af_alt_proj = {proj['af_alt_proj']:.4f} "
          f"(panel sum={proj['sum_panel']:.4f}; "
          f"{proj['n_carriers_in_panel']}/{proj['n_carriers']} carriers in panel)")

    # 4. freqk comparison
    af_freqk = None
    if args.freqk_af is not None and args.freqk_af.exists():
        _, af_freqk = load_freqk_af(args.freqk_af)
        print(f"freqk:      af_alt_freqk = {af_freqk:.4f}")

    # 5. Write summary
    out = pd.DataFrame([{
        "ecotype_freq_file": str(args.ecotype_freq),
        "n_samples": args.n_samples,
        "sv_freq_input": args.sv_freq,
        "n_sv_carriers": n_sv,
        "truth_af_alt": truth_af_alt,
        "hapfire_panel_sum": proj["sum_panel"],
        "af_alt_hapfire_proj": proj["af_alt_proj"],
        "af_alt_freqk": af_freqk,
        "err_hapfire": proj["af_alt_proj"] - truth_af_alt,
        "err_freqk": (af_freqk - truth_af_alt) if af_freqk is not None else None,
    }])
    args.out.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out, sep="\t", index=False)
    print(f"wrote: {args.out}")


if __name__ == "__main__":
    main()
