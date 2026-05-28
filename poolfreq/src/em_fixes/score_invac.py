"""Score the invac gamma-sweep against g0 truth and SEEDMIX/hapFIRE references.

Reads scratch/h_fixes/invac/{g0_<sim>,seedmix_S<n>}.invac.npz
Outputs a per-gamma table for the 5 g0 sims + SEEDMIX 8-rep average.
"""
from __future__ import annotations
import csv, collections, glob, json, os
import numpy as np

REPO = "/carnegie/nobackup/scratch/tbellagio/hapfire_sv"
INVAC = f"{REPO}/scratch/h_fixes/invac"
SPLIT = json.load(open(f"{REPO}/data/founder_split_cactus_pg.json"))
CACT = set(map(str, SPLIT["cactus"]))
HAPFIRE = "/carnegie/nobackup/scratch/xwu/GrENE_net/hapFIRE_frequencies/seed_mix"

GAMMAS = ["0", "0.25", "0.5", "0.75", "1"]
G0_SIMS = ["g0_n231_rep0_rand", "g0_n200_rep0_rand",
           "g0_n50_rep0_cact", "g0_n50_rep1_bal", "g0_n50_rep2_pg"]
G0_TARGET = {"g0_n231_rep0_rand": 0.346, "g0_n200_rep0_rand": 0.345,
             "g0_n50_rep0_cact": 0.80, "g0_n50_rep1_bal": 0.34,
             "g0_n50_rep2_pg": 0.10}


def load_truth():
    bysim = collections.defaultdict(dict)
    with open(f"{REPO}/scratch/g0_sweep_per_founder.tsv") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            if row["tag"] == "filt2":
                bysim[row["sim"]][row["founder"]] = float(row["truth"])
    return bysim


def cactus_mass(h, founders):
    return float(sum(h[i] for i, f in enumerate(founders) if f in CACT))


def main():
    truth = load_truth()
    F = 231
    print("=" * 110)
    print("G0 SIMS  (target cactus_mass in parens)")
    print("=" * 110)
    hdr = f"{'sim':<20} {'gamma':>6} {'cact_mass':>10} {'MAE':>9} {'RMSE':>9} {'max|d|':>8} {'leak':>8}"
    print(hdr)
    g0_summary = collections.defaultdict(dict)
    for sim in G0_SIMS:
        d = np.load(f"{INVAC}/g0_{sim}.invac.npz", allow_pickle=True)
        founders = np.asarray(d["founders"]).astype(str)
        tvec = np.array([truth[sim][f] for f in founders])
        zero = tvec == 0
        tgt = G0_TARGET[sim]
        for g in GAMMAS:
            h = d[f"h_gamma{g}"]
            cm = cactus_mass(h, founders)
            mae = float(np.mean(np.abs(h - tvec)))
            rmse = float(np.sqrt(np.mean((h - tvec) ** 2)))
            mx = float(np.max(np.abs(h - tvec)))
            leak = float(h[zero].sum())
            g0_summary[sim][g] = dict(cm=cm, mae=mae, rmse=rmse, mx=mx, leak=leak)
            print(f"{sim+' ('+format(tgt,'.3f')+')':<20} {g:>6} {cm:>10.4f} "
                  f"{mae:>9.5f} {rmse:>9.5f} {mx:>8.4f} {leak:>8.4f}")
        print("-" * 110)

    # --- SEEDMIX: 8-rep average ---
    print("\n" + "=" * 110)
    print("SEEDMIX  (8-rep average h; target cactus_mass ~0.35)")
    print("=" * 110)
    # hapFIRE 8-rep avg reference
    hf = {}
    for n in range(1, 9):
        with open(f"{HAPFIRE}/s{n}_ecotype_frequency.txt") as f:
            for line in f:
                fid, val = line.split()
                hf.setdefault(fid, []).append(float(val))
    # load reps
    reps = {}
    founders = None
    for n in range(1, 9):
        d = np.load(f"{INVAC}/seedmix_S{n}.invac.npz", allow_pickle=True)
        if founders is None:
            founders = np.asarray(d["founders"]).astype(str)
        reps[n] = d
    hf_vec = np.array([np.mean(hf[f]) for f in founders])
    unif = np.full(F, 1.0 / F)

    print(f"{'gamma':>6} {'cact_mass':>10} {'RMSE_unif':>11} {'RMSE_hapF':>11} "
          f"{'mean_perrep_std':>16}")
    sm_summary = {}
    for g in GAMMAS:
        H = np.array([reps[n][f"h_gamma{g}"] for n in range(1, 9)])  # (8, F)
        havg = H.mean(0)
        cm = cactus_mass(havg, founders)
        rmse_unif = float(np.sqrt(np.mean((havg - unif) ** 2)))
        rmse_hf = float(np.sqrt(np.mean((havg - hf_vec) ** 2)))
        perrep_std = float(H.std(0).mean())
        sm_summary[g] = dict(cm=cm, rmse_unif=rmse_unif, rmse_hf=rmse_hf,
                             perrep_std=perrep_std)
        print(f"{g:>6} {cm:>10.4f} {rmse_unif:>11.5f} {rmse_hf:>11.5f} "
              f"{perrep_std:>16.6f}")

    print("\n" + "=" * 110)
    print("BASELINES for reference:")
    print("  filt2       : SEEDMIX cact=0.518 RMSE_unif=0.00318 | g0_n231 RMSE=0.00283 cact=0.503")
    print("  subsampMedian: SEEDMIX cact=0.370 RMSE_unif=0.00255 | g0_n231 RMSE=0.00225 cact=0.353")
    print("=" * 110)

    # dump json
    out = {"g0": {s: g0_summary[s] for s in g0_summary}, "seedmix": sm_summary}
    with open(f"{INVAC}/score_summary.json", "w") as f:
        json.dump(out, f, indent=2)
    print(f"\nwrote {INVAC}/score_summary.json")


if __name__ == "__main__":
    main()
