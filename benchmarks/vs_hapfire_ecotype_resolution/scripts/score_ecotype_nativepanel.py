#!/usr/bin/env python3
"""Ecotype-RESOLUTION metrics for kMate vs hapFIRE, each on its OWN native panel.

WHY THIS EXISTS
---------------
The original `score_ecotype.py` + `run_hapfire_ecotype.sbatch` pair fed hapFIRE the
SHARED panel `vs_vg_giraffe/work/shared_snps_<panel>_Chr1.vcf.gz`. hapFIRE
errors on missing GTs, so that panel is built with `MISSING '.' -> 0|0` (impute REF),
handing hapFIRE false REF homozygotes at ~7.1% of records. kMate is never exposed to
this: its `var_pa` ships a `var_called` mask that excludes exactly those sites. One
tool may say "unknown", the other is forced to guess REF -- so a shared panel is not a
shared *input*, and the resulting hapFIRE numbers are not attributable to hapFIRE.

That shared-panel arm is archived under
`benchmarks/archive/shared_panel_hapfire_comparison_invalid/`. This script is its
native-panel replacement.

This is the same correction, and the same corrected design, as the 2026-07-15
panel-mismatch fix documented in `../../vs_hapfire/PANEL_MISMATCH_BUG.md`:
each tool runs on the panel its reads were simulated from.

  kMate   : arch3 panel, arch3-simulated reads
            (p231/results/kmate_chrom_poolsize_depth/, existing -- no rerun)
  hapFIRE : greneNet_final_v1.1 panel, greneNet-derived reads
            (vs_hapfire/results/greneNet_fair/, existing -- no rerun)

Matched founder draw + seed => identical pool composition, so founder-frequency
recovery is directly comparable. The truth (`pool_weights.tsv`) is legitimately
SHARED here even though the panels are not, because a founder's pool weight is a
property of the draw, not of any variant-calling lineage.

WHAT THIS ADDS OVER vs_hapfire
------------------------------------
`vs_hapfire/scripts/score_hapfire_vs_kmate.py` already covers h_R2 / h_RMSE /
n_found on this same native-panel basis. It does NOT compute the ecotype-RESOLUTION
metrics below (eff_n, recall, precision, mass_on_true, top_n_jaccard), which are the
reason `vs_hapfire_ecotype_resolution/` exists as a separate benchmark. Those are what this script
adds -- it deliberately does not duplicate the h_R2/speed columns.

Metrics per (tool, N, depth, seed):
  true_n        founders actually in the pool (weight > 0)
  eff_n         effective #founders = 1 / sum(h^2)          (kMate's headline readout)
  n_detected    #founders with h >= 0.5/true_n (half the equal share)
  recall        #true founders detected / true_n
  precision     #detected that are true / n_detected
  mass_on_true  fraction of estimated mass on true founders (1 = perfect)
  top_n_jaccard Jaccard(top-true_n founders by h, true set)   [threshold-free]
  h_rmse        RMSE(estimate vs true weights) over all panel founders

Run with the `basic` env (matplotlib/statsmodels hang in `plotting`):
    /global/home/users/tbellg/miniforge3/envs/basic/bin/python \
        benchmarks/vs_hapfire_ecotype_resolution/scripts/score_ecotype_nativepanel.py
"""
import os
import numpy as np
import pandas as pd

ROOT = "/global/scratch/users/tbellg/kmate"
HF_RES = f"{ROOT}/benchmarks/vs_hapfire/results/greneNet_fair"
HF_SIMS = f"{ROOT}/benchmarks/vs_hapfire/sims_greneNet"
P231 = f"{ROOT}/benchmarks/p231"
OUT = f"{ROOT}/benchmarks/vs_hapfire_ecotype_resolution/results"

# Same grid as the corrected vs_hapfire comparison: every N at depth {1,10},
# plus a depth extension {30,50} at N=50 only. Encoded exactly (rather than as a
# full N x depth cross product) so a genuinely missing run is reported as missing
# instead of being buried among ~50 cells that were never meant to exist.
POOL_SIZES = [2, 5, 20, 50, 150, 231]
SEEDS = [42, 43, 44, 45, 46]
GRID = ([(n, d) for n in POOL_SIZES for d in (1, 10)]
        + [(50, d) for d in (30, 50)])


def resolution_metrics(founders, est, tw, true_set):
    """Ecotype-resolution metrics for one estimate vector.

    founders : array of founder ids (str), aligned with est
    est      : estimated founder frequencies (need not sum to exactly 1)
    tw       : dict founder -> true weight
    true_set : set of founders with weight > 0
    """
    est = np.asarray(est, dtype=float)
    tot = est.sum()
    p = est / tot if tot > 0 else est
    true_n = len(true_set)

    # eff_n is only meaningful on a normalized vector, hence p not est.
    eff_n = 1.0 / np.sum(p ** 2) if np.sum(p ** 2) > 0 else np.nan

    thr = 0.5 / true_n if true_n else np.inf
    det_mask = p >= thr
    detected = set(np.asarray(founders)[det_mask])
    n_detected = int(det_mask.sum())

    recall = len(detected & true_set) / true_n if true_n else np.nan
    precision = len(detected & true_set) / n_detected if n_detected else np.nan

    is_true = np.array([f in true_set for f in founders])
    mass_on_true = float(p[is_true].sum()) if tot > 0 else np.nan

    order = np.argsort(-p)
    top = set(np.asarray(founders)[order[:true_n]])
    union = top | true_set
    top_n_jaccard = len(top & true_set) / len(union) if union else np.nan

    truth_vec = np.array([tw.get(f, 0.0) for f in founders])
    h_rmse = float(np.sqrt(np.mean((p - truth_vec) ** 2)))

    return dict(true_n=true_n, eff_n=eff_n, n_detected=n_detected, recall=recall,
                precision=precision, mass_on_true=mass_on_true,
                top_n_jaccard=top_n_jaccard, h_rmse=h_rmse)


def main():
    rows, missing = [], []
    for n, cov in GRID:
            for seed in SEEDS:
                hf_ef = (f"{HF_RES}/n{n}_cov{cov}_s{seed}/"
                         f"greneNet_hapfire_n{n}_cov{cov}_s{seed}_ecotype_frequency.txt")
                km_npz = (f"{P231}/results/kmate_chrom_poolsize_depth/n{n}_cov{cov}_s{seed}/"
                          f"p231_chrom_psd_n{n}_cov{cov}_s{seed}.h_per_chrom.npz")
                truth_p = f"{HF_SIMS}/cov{cov}_n{n}_g0_s{seed}_greneNet_chr1/pool_weights.tsv"
                if not all(os.path.exists(p) for p in (hf_ef, km_npz, truth_p)):
                    missing.append(f"n{n}_cov{cov}_s{seed}")
                    continue

                twdf = pd.read_csv(truth_p, sep="\t")
                twdf["founder"] = twdf["founder"].astype(str)
                tw = dict(zip(twdf["founder"], twdf["weight"]))
                true_set = set(twdf[twdf.weight > 0]["founder"])

                ef = pd.read_csv(hf_ef, sep="\t", header=None, names=["founder", "freq"])
                ef["founder"] = ef["founder"].astype(str)
                rows.append(dict(tool="hapFIRE", panel="greneNet", N=n, depth=cov, seed=seed,
                                 **resolution_metrics(ef["founder"].values,
                                                      ef["freq"].values, tw, true_set)))

                hd = np.load(km_npz, allow_pickle=True)
                rows.append(dict(tool="kMate", panel="arch3", N=n, depth=cov, seed=seed,
                                 **resolution_metrics(hd["founders"].astype(str),
                                                      hd["Chr1"].astype(float), tw, true_set)))

    if not rows:
        raise SystemExit("no matched (hapFIRE, kMate, truth) triples found")

    df = pd.DataFrame(rows).sort_values(["N", "depth", "seed", "tool"])
    os.makedirs(OUT, exist_ok=True)
    path = f"{OUT}/ecotype_nativepanel_table.tsv"
    df.to_csv(path, sep="\t", index=False)
    print(f"wrote {len(df)} rows ({len(df)//2}/{len(GRID)*len(SEEDS)} designed pools) "
          f"-> {path}")
    # Never let dropped cells pass silently as "covered everything".
    print(f"missing runs: {len(missing)}" + (f" -> {', '.join(missing)}" if missing else "")
          + "\n")

    summ = (df.groupby(["tool", "N", "depth"])[
                ["recall", "precision", "mass_on_true", "top_n_jaccard", "h_rmse"]]
              .mean().round(4).reset_index())
    print(summ.to_string(index=False))


if __name__ == "__main__":
    main()
