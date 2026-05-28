"""
Score the balanced sub-design EM (and baselines) for the bias-free-anchor test.

Reads:
  scratch/h_fixes/balanced/balanced_<SIM>.npz   (h from balanced SNP sub-design)
  scratch/g0_sweep_h_test/{filt2,f2subsamp}_<SIM>.h_sweep.npz   (baselines, alpha=0)
  scratch/seedmix_h_test/{filt2,f2subsamp}_S{i}.h_sweep.npz
  scratch/g0_sweep_per_founder.tsv   (g0 truth)
  hapFIRE seed_mix s{i}_ecotype_frequency.txt

Writes:
  scratch/h_fixes/balanced/balanced_scores.tsv
and prints the score tables.
"""
from __future__ import annotations
import json, os
import numpy as np
import pandas as pd

ROOT = "/carnegie/nobackup/scratch/tbellagio/hapfire_sv"
BAL = os.path.join(ROOT, "scratch/h_fixes/balanced")
G0D = os.path.join(ROOT, "scratch/g0_sweep_h_test")
SMD = os.path.join(ROOT, "scratch/seedmix_h_test")
HAPFIRE = "/carnegie/nobackup/scratch/xwu/GrENE_net/hapFIRE_frequencies/seed_mix"

G0_SIMS = ["g0_n231_rep0_rand", "g0_n200_rep0_rand",
           "g0_n50_rep0_cact", "g0_n50_rep1_bal", "g0_n50_rep2_pg"]
G0_TARGET = {"g0_n231_rep0_rand": 0.346, "g0_n200_rep0_rand": 0.345,
             "g0_n50_rep0_cact": 0.800, "g0_n50_rep1_bal": 0.340,
             "g0_n50_rep2_pg": 0.100}


def load_classes():
    cls = json.load(open(os.path.join(ROOT, "data/founder_split_cactus_pg.json")))
    return set(cls["cactus"]), set(cls["PG"])


def h_baseline(prefix_dir, prefix, sim):
    """Load alpha=0 (standard unweighted EM, full or subsamp design) h."""
    p = os.path.join(prefix_dir, f"{prefix}_{sim}.h_sweep.npz")
    if not os.path.exists(p):
        return None, None
    d = np.load(p, allow_pickle=True)
    fo = np.asarray(d["founders"]).astype(str)
    h = np.asarray(d["h_per_alpha"][0], dtype=np.float64)
    return fo, h


def h_balanced(sim):
    p = os.path.join(BAL, f"balanced_{sim}.npz")
    if not os.path.exists(p):
        return None, None
    d = np.load(p, allow_pickle=True)
    return np.asarray(d["founders"]).astype(str), np.asarray(d["h"], dtype=np.float64)


def metrics_vs_truth(fo, h, truth_map, cact_set):
    """truth_map: dict founder->truth h. Returns dict of metrics."""
    truth = np.array([truth_map.get(str(f), 0.0) for f in fo])
    isc = np.array([str(f) in cact_set for f in fo])
    diff = h - truth
    cm = float(h[isc].sum())
    # leakage: estimated mass placed on founders with truth==0
    leak = float(h[truth == 0].sum())
    return {
        "cactus_mass": cm,
        "MAE": float(np.abs(diff).mean()),
        "RMSE": float(np.sqrt((diff**2).mean())),
        "max_abs": float(np.abs(diff).max()),
        "leakage": leak,
    }


def main():
    cact_set, pg_set = load_classes()
    rows = []

    # ---------- g0 truth ----------
    tdf = pd.read_csv(os.path.join(ROOT, "scratch/g0_sweep_per_founder.tsv"), sep="\t")
    for sim in G0_SIMS:
        sub = tdf[tdf["sim"] == sim].drop_duplicates("founder")
        truth_map = dict(zip(sub["founder"].astype(str), sub["truth"].astype(float)))
        tgt = G0_TARGET[sim]
        for method, loader in [
            ("balanced_SNP", lambda s=sim: h_balanced(s)),
            ("filt2_full",   lambda s=sim: h_baseline(G0D, "filt2", s)),
            ("subsampMedian", lambda s=sim: h_baseline(G0D, "f2subsamp", s)),
        ]:
            fo, h = loader()
            if h is None:
                continue
            m = metrics_vs_truth(fo, h, truth_map, cact_set)
            rows.append(dict(group="g0", sim=sim, method=method,
                             target_cactus=tgt, **m))

    # ---------- SEEDMIX ----------
    # hapFIRE 8-rep avg (231 founders)
    hf = {}
    for i in range(1, 9):
        d = pd.read_csv(os.path.join(HAPFIRE, f"s{i}_ecotype_frequency.txt"),
                        sep="\t", header=None, names=["founder", "freq"])
        hf[i] = dict(zip(d["founder"].astype(str), d["freq"].astype(float)))
    # union founder order from one rep
    hf_founders = list(pd.read_csv(os.path.join(HAPFIRE, "s1_ecotype_frequency.txt"),
                                   sep="\t", header=None, names=["f", "v"])["f"].astype(str))
    HF = np.array([[hf[i].get(f, 0.0) for f in hf_founders] for i in range(1, 9)])
    hf_avg = HF.mean(0)
    isc_hf = np.array([f in cact_set for f in hf_founders])
    uniform = np.full(len(hf_founders), 1.0 / len(hf_founders))

    sm_methods = {"balanced_SNP": [], "filt2_full": [], "subsampMedian": [],
                  "hapFIRE": []}
    sm_founders_ref = None
    for i in range(1, 9):
        sim = f"S{i}"
        for method, loader in [
            ("balanced_SNP", lambda s=sim: h_balanced(s)),
            ("filt2_full",   lambda s=sim: h_baseline(SMD, "filt2", s)),
            ("subsampMedian", lambda s=sim: h_baseline(SMD, "f2subsamp", s)),
        ]:
            fo, h = loader()
            if h is None:
                continue
            sm_founders_ref = fo
            # align to founder order
            order = {f: j for j, f in enumerate(fo)}
            sm_methods[method].append(h)
        # hapFIRE for this rep aligned to our founder order
    # build per-method 8-rep stacks aligned to our founders
    if sm_founders_ref is not None:
        our_fo = sm_founders_ref
        isc_our = np.array([f in cact_set for f in our_fo])
        # hapFIRE aligned to our founder order
        hf_avg_our = np.array([hf_avg[hf_founders.index(f)] if f in hf_founders else 0.0
                               for f in our_fo])
        unif_our = np.full(len(our_fo), 1.0 / len(our_fo))
        for method in ["balanced_SNP", "filt2_full", "subsampMedian"]:
            stk = sm_methods[method]
            if len(stk) == 0:
                continue
            M = np.array(stk)  # 8 x 231
            avg = M.mean(0)
            cm_per = M[:, isc_our].sum(1)            # per-rep cactus mass
            rmse_unif = float(np.sqrt(((avg - unif_our)**2).mean()))
            rmse_hf = float(np.sqrt(((avg - hf_avg_our)**2).mean()))
            # per-rep std of cactus mass + per-founder std
            rows.append(dict(group="SEEDMIX", sim="S1-8_avg", method=method,
                             target_cactus=0.35,
                             cactus_mass=float(cm_per.mean()),
                             cactus_mass_std=float(cm_per.std()),
                             RMSE_vs_uniform=rmse_unif,
                             RMSE_vs_hapFIRE=rmse_hf,
                             perfounder_std=float(M.std(0).mean()),
                             n_reps=len(stk)))
        # hapFIRE self row
        cm_hf = HF[:, isc_hf].sum(1)
        rows.append(dict(group="SEEDMIX", sim="S1-8_avg", method="hapFIRE",
                         target_cactus=0.35,
                         cactus_mass=float(cm_hf.mean()),
                         cactus_mass_std=float(cm_hf.std()),
                         RMSE_vs_uniform=float(np.sqrt(((hf_avg - 1.0/len(hf_founders))**2).mean())),
                         RMSE_vs_hapFIRE=0.0,
                         perfounder_std=float(HF.std(0).mean()),
                         n_reps=8))

    out = pd.DataFrame(rows)
    op = os.path.join(BAL, "balanced_scores.tsv")
    out.to_csv(op, sep="\t", index=False)

    pd.set_option("display.width", 200)
    pd.set_option("display.max_columns", 30)
    pd.set_option("display.float_format", lambda x: f"{x:.4f}")
    print("\n================ G0 (vs per-founder truth) ================")
    g = out[out.group == "g0"][["sim", "method", "target_cactus", "cactus_mass",
                                "MAE", "RMSE", "max_abs", "leakage"]]
    print(g.to_string(index=False))
    print("\n================ SEEDMIX (8-rep avg) ================")
    s = out[out.group == "SEEDMIX"][["method", "target_cactus", "cactus_mass",
                                     "cactus_mass_std", "RMSE_vs_uniform",
                                     "RMSE_vs_hapFIRE", "perfounder_std"]]
    print(s.to_string(index=False))
    print(f"\nWrote {op}")


if __name__ == "__main__":
    main()
