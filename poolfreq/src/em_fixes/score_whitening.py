"""
Score whitening IRLS/GLS h vs baselines.
  g0 : 5 sims, truth from scratch/g0_sweep_per_founder.tsv (tag in {filt2,subsamp})
       metrics per sim: cactus_mass, MAE, RMSE, max|h-truth|, leakage
       (leakage = mass placed on founders whose truth==0)
  SEEDMIX : 8 reps; cactus_mass; RMSE vs uniform(1/231); RMSE vs hapFIRE 8-rep avg;
            per-rep std of cactus_mass.
Compares whitening to filt2 baseline and subsampMedian baseline (the bar to beat).
"""
from __future__ import annotations
import os, csv, json, glob
import numpy as np
import collections

ROOT = "/carnegie/nobackup/scratch/tbellagio/hapfire_sv"
WDIR = os.path.join(ROOT, "scratch/h_fixes/whitening")
CLASSES = json.load(open(os.path.join(ROOT, "data/founder_split_cactus_pg.json")))
CACT = set(CLASSES["cactus"]); PG = set(CLASSES["PG"])
HF_DIR = "/carnegie/nobackup/scratch/xwu/GrENE_net/hapFIRE_frequencies/seed_mix"

G0_SIMS = ["g0_n231_rep0_rand", "g0_n200_rep0_rand",
           "g0_n50_rep0_cact", "g0_n50_rep1_bal", "g0_n50_rep2_pg"]
SEEDMIX = [f"S{i}" for i in range(1, 9)]


def load_truth_and_baselines():
    """Return {sim: {'truth':dict, 'filt2':dict, 'subsamp':dict}} (founder->val)."""
    rows = list(csv.DictReader(open(os.path.join(ROOT, "scratch/g0_sweep_per_founder.tsv")),
                               delimiter="\t"))
    out = collections.defaultdict(lambda: collections.defaultdict(dict))
    for r in rows:
        sim = r["sim"]; f = r["founder"]; tag = r["tag"]
        out[sim]["truth"][f] = float(r["truth"])
        out[sim][tag][f] = float(r["h"])
    return out


def load_whiten(sample):
    p = os.path.join(WDIR, f"{sample}.whiten.npz")
    if not os.path.exists(p):
        return None
    d = np.load(p, allow_pickle=True)
    founders = np.asarray(d["founders"]).astype(str)
    h = np.asarray(d["h"]).astype(float)
    return dict(zip(founders, h)), d


def g0_metrics(hdict, truthdict):
    fs = list(truthdict)
    h = np.array([hdict.get(f, 0.0) for f in fs])
    t = np.array([truthdict[f] for f in fs])
    cm = sum(hdict.get(f, 0.0) for f in fs if f in CACT)
    mae = float(np.mean(np.abs(h - t)))
    rmse = float(np.sqrt(np.mean((h - t) ** 2)))
    mx = float(np.max(np.abs(h - t)))
    leak = float(sum(hdict.get(f, 0.0) for f in fs if truthdict[f] == 0.0))
    return dict(cactus_mass=cm, MAE=mae, RMSE=rmse, maxabs=mx, leakage=leak)


def load_hf_avg():
    accum = collections.defaultdict(list)
    for i in range(1, 9):
        f = os.path.join(HF_DIR, f"s{i}_ecotype_frequency.txt")
        for line in open(f):
            p = line.split()
            accum[p[0]].append(float(p[1]))
    return {k: float(np.mean(v)) for k, v in accum.items()}


def main():
    data = load_truth_and_baselines()
    print("=" * 100)
    print("G0 SIMS  (truth-based)   target cactus_mass: n231=0.346 n200=0.345 n50_cact=0.80 n50_bal=0.34 n50_pg=0.10")
    print("=" * 100)
    targets = {"g0_n231_rep0_rand": 0.346, "g0_n200_rep0_rand": 0.345,
               "g0_n50_rep0_cact": 0.80, "g0_n50_rep1_bal": 0.34, "g0_n50_rep2_pg": 0.10}
    hdr = f"{'sim':22s} {'method':10s} {'cactus':>7s}(tgt) {'MAE':>9s} {'RMSE':>9s} {'max|d|':>8s} {'leak':>7s}"
    print(hdr); print("-" * len(hdr))
    for sim in G0_SIMS:
        tr = data[sim]["truth"]
        tgt = targets[sim]
        for method, hd in [("filt2", data[sim].get("filt2")),
                           ("subsampMed", data[sim].get("subsamp")),
                           ("whiten", (load_whiten(sim) or (None,))[0])]:
            if hd is None:
                print(f"{sim:22s} {method:10s}  (missing)")
                continue
            m = g0_metrics(hd, tr)
            print(f"{sim:22s} {method:10s} {m['cactus_mass']:7.4f}({tgt:.2f}) "
                  f"{m['MAE']:9.6f} {m['RMSE']:9.6f} {m['maxabs']:8.4f} {m['leakage']:7.4f}")
        print()

    print("=" * 100)
    print("SEEDMIX (8 reps)   target cactus_mass ~0.35 ; baselines filt2=0.518 subsampMed=0.370")
    print("=" * 100)
    hf_avg = load_hf_avg()
    fs_panel = sorted(hf_avg)
    unif = 1.0 / len(fs_panel)
    hf_vec = np.array([hf_avg[f] for f in fs_panel])

    cms = []; rmse_unif = []; rmse_hf = []
    hmats = []
    for s in SEEDMIX:
        w = load_whiten(s)
        if w is None:
            print(f"  {s}: missing"); continue
        hd = w[0]
        cm = sum(hd.get(f, 0.0) for f in fs_panel if f in CACT)
        hv = np.array([hd.get(f, 0.0) for f in fs_panel])
        cms.append(cm)
        rmse_unif.append(np.sqrt(np.mean((hv - unif) ** 2)))
        rmse_hf.append(np.sqrt(np.mean((hv - hf_vec) ** 2)))
        hmats.append(hv)
    if cms:
        cms = np.array(cms)
        print(f"  whiten SEEDMIX: cactus_mass mean={cms.mean():.4f}  std(per-rep)={cms.std():.4f}  "
              f"range[{cms.min():.4f},{cms.max():.4f}]")
        print(f"  whiten RMSE vs uniform(1/231)  mean={np.mean(rmse_unif):.5f}")
        print(f"  whiten RMSE vs hapFIRE 8-rep avg mean={np.mean(rmse_hf):.5f}")
        print()
        print("  per-rep cactus_mass: " + " ".join(f"{s}={c:.4f}" for s, c in zip(SEEDMIX, cms)))
    print()
    print("  BASELINE REFERENCE (SEEDMIX): filt2 cactus=0.518 RMSE_unif=0.00318 ; "
          "subsampMed cactus=0.370 RMSE_unif=0.00255")


if __name__ == "__main__":
    main()
