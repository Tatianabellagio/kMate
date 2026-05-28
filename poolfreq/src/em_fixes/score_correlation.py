#!/usr/bin/env python3
"""
Score the correlation (c2 overdispersion) EM outputs vs filt2 + subsampMedian.

g0: per-sim MAE/RMSE/maxErr/leakage + cactus-class est vs true mass, against
    scratch/g0_sweep_per_founder.tsv truth. Sweeps rho.
SEEDMIX: 8-rep avg cactus_mass, RMSE vs uniform (1/231), RMSE vs hapFIRE 8-rep
    avg, per-rep std. Sweeps rho.

Reads correlation outputs from scratch/h_fixes/correlation/*.corr.npz and the
existing filt2/subsamp baselines in scratch/{g0_sweep_h_test,seedmix_h_test}/.
"""
import csv, json, math, os
from collections import defaultdict
from pathlib import Path
import numpy as np

ROOT = str(Path(__file__).resolve().parents[3])
CORR = f"{ROOT}/scratch/h_fixes/correlation"
G0_OUT = f"{ROOT}/scratch/g0_sweep_h_test"
SM_OUT = f"{ROOT}/scratch/seedmix_h_test"
PF = f"{ROOT}/scratch/g0_sweep_per_founder.tsv"
SPLIT = f"{ROOT}/data/founder_split_cactus_pg.json"
HAPFIRE = "/global/scratch/projects/fc_moilab/projects/grenenet-phase1/frequency/hapFIRE_frequencies/seed_mix"

G0_SIMS = ["g0_n231_rep0_rand", "g0_n200_rep0_rand",
           "g0_n50_rep0_cact", "g0_n50_rep1_bal", "g0_n50_rep2_pg"]

split = json.load(open(SPLIT))
CACTUS = set(split["cactus"])


def load_baseline_h(out_dir, tag, sim):
    p = f"{out_dir}/{tag}_{sim}.h_sweep.npz"
    if not os.path.exists(p):
        return None, None
    d = np.load(p, allow_pickle=True)
    return (np.asarray(d["founders"]).astype(str),
            np.asarray(d["h_per_alpha"])[0].astype(float))


def load_corr(tag):
    p = f"{CORR}/{tag}.corr.npz"
    if not os.path.exists(p):
        return None
    d = np.load(p, allow_pickle=True)
    return (np.asarray(d["founders"]).astype(str),
            np.asarray(d["rhos"]).astype(float),
            np.asarray(d["h_per_rho"]).astype(float))


# ---------------- g0 ----------------
truth = defaultdict(dict)
side = defaultdict(dict)
with open(PF) as f:
    for d in csv.DictReader(f, delimiter="\t"):
        truth[d["sim"]][d["founder"]] = float(d["truth"])
        side[d["sim"]][d["founder"]] = d["side"]


def g0_score(founders, h, sim):
    t = np.array([truth[sim].get(f, np.nan) for f in founders])
    ok = ~np.isnan(t)
    h2, t2, fo = h[ok], t[ok], founders[ok]
    ae = np.abs(h2 - t2)
    mae = ae.mean(); rmse = math.sqrt((ae**2).mean()); mx = ae.max()
    leak = float(h2[t2 == 0].sum())
    cmask = np.array([f in CACTUS for f in fo])
    return mae, rmse, mx, leak, float(h2[cmask].sum()), float(t2[cmask].sum())


print("=" * 96)
print("G0 SIMS  (truth cactus targets: n231=0.346 n200=0.345 n50_cact=0.80 "
      "n50_bal=0.34 n50_pg=0.10)")
print("=" * 96)
hdr = f"{'sim':<19}{'method':<14}{'MAE':>9}{'RMSE':>9}{'maxErr':>8}{'leak':>7}{'C_est':>7}{'C_tru':>7}"
for sim in G0_SIMS:
    print("-" * 96)
    print(hdr)
    rows = []
    # baselines
    for tag in ["filt2", "subsamp"]:
        fo, h = load_baseline_h(G0_OUT, tag, sim)
        if h is not None:
            rows.append((tag, g0_score(fo, h, sim)))
    # correlation rhos
    c = load_corr(sim)
    if c is not None:
        fo, rhos, H = c
        for i, rho in enumerate(rhos):
            rows.append((f"corr_rho{rho:g}", g0_score(fo, H[i], sim)))
    base_rmse = rows[0][1][1] if rows else None
    for name, (mae, rmse, mx, leak, ce, ct) in rows:
        flag = f"  RMSE {(rmse/base_rmse-1)*100:+.0f}% vs filt2" if base_rmse else ""
        print(f"{sim:<19}{name:<14}{mae:>9.5f}{rmse:>9.5f}{mx:>8.4f}"
              f"{leak:>7.3f}{ce:>7.3f}{ct:>7.3f}{flag}")
print()

# ---------------- SEEDMIX ----------------
print("=" * 96)
print("SEEDMIX 8-rep  (truth cactus mass ~0.35; lower RMSE_unif/RMSE_hapFIRE = better)")
print("=" * 96)

# hapFIRE 8-rep avg ecotype freq, keyed by founder id
hf = {}
for s in range(1, 9):
    p = f"{HAPFIRE}/s{s}_ecotype_frequency.txt"
    d = {}
    with open(p) as f:
        for line in f:
            a = line.split()
            if len(a) >= 2:
                d[a[0]] = float(a[1])
    hf[s] = d
# build hapFIRE per-founder mean across reps using a reference founder order later


def sm_methods():
    """Yield (method_name, founders, H[reps x F]) for each method."""
    # filt2 + subsamp baselines
    for tag in ["filt2", "subsamp"]:
        Hs, fo_ref = [], None
        ok = True
        for s in range(1, 9):
            fo, h = load_baseline_h(SM_OUT, tag, f"S{s}")
            if h is None:
                ok = False; break
            fo_ref = fo; Hs.append(h)
        if ok:
            yield tag, fo_ref, np.array(Hs)
    # correlation: one file per rep, each with rho axis -> regroup by rho
    corr0 = load_corr("seedmix_S1")
    if corr0 is not None:
        _, rhos, _ = corr0
        for i, rho in enumerate(rhos):
            Hs, fo_ref = [], None
            ok = True
            for s in range(1, 9):
                c = load_corr(f"seedmix_S{s}")
                if c is None:
                    ok = False; break
                fo, _, H = c
                fo_ref = fo; Hs.append(H[i])
            if ok:
                yield f"corr_rho{rho:g}", fo_ref, np.array(Hs)


print(f"{'method':<14}{'C_mass':>9}{'C_std':>8}{'RMSE_unif':>11}{'RMSE_hapF':>11}"
      f"{'reproStd':>10}")
print("-" * 96)
base_unif = None
for name, fo, H in sm_methods():
    cmask = np.array([f in CACTUS for f in fo])
    cmass = H[:, cmask].sum(axis=1)         # per-rep cactus mass
    h_avg = H.mean(axis=0)
    unif = 1.0 / len(fo)
    rmse_unif = math.sqrt(((h_avg - unif) ** 2).mean())
    # hapFIRE per-founder 8-rep avg aligned to fo
    hf_avg = np.array([np.mean([hf[s].get(f, np.nan) for s in range(1, 9)])
                       for f in fo])
    m = ~np.isnan(hf_avg)
    rmse_hf = math.sqrt(((h_avg[m] - hf_avg[m]) ** 2).mean())
    # per-rep reproducibility: mean over founders of std across reps
    repro = H.std(axis=0).mean()
    if base_unif is None:
        base_unif = rmse_unif
    print(f"{name:<14}{cmass.mean():>9.4f}{cmass.std():>8.4f}"
          f"{rmse_unif:>11.5f}{rmse_hf:>11.5f}{repro:>10.6f}")
print()
print("Baselines (from task spec): filt2 SEEDMIX cactus 0.518 RMSE_unif 0.00318;")
print("                            subsampMedian cactus 0.370 RMSE_unif 0.00255.")
