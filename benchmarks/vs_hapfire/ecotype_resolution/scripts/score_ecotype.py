#!/usr/bin/env python3
"""Score kMate --h-only against the true founder mixture for an ecotype-count pool.

Truth = pool_weights.tsv (founder, count, weight); true founders = weight>0.
Estimate = <out>.h_per_chrom.npz (founders array + per-chrom h, here Chr1).

Metrics (one row appended per pool):
  true_n        number of founders actually in the pool
  eff_n         estimated effective #founders = 1 / sum(h^2)   (kMate's readout)
  n_detected    #founders with h >= thr  (thr = 0.5/true_n, i.e. half the true share)
  recall        #true founders detected / true_n
  precision     #detected that are true / n_detected
  mass_on_true  fraction of estimated h mass on the true founders (1 = perfect)
  top_n_jaccard Jaccard(top-true_n founders by h, true set)
  h_rmse        RMSE(h_est vs true weights) over all panel founders
"""
import argparse, sys
from pathlib import Path
import numpy as np, pandas as pd


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--h", help="kMate <out>.h_per_chrom.npz")
    ap.add_argument("--hapfire-ecotype", help="hapFIRE _ecotype_frequency.txt (founder<TAB>freq)")
    ap.add_argument("--tool", default="kMate")
    ap.add_argument("--weights", required=True)
    ap.add_argument("--pool", required=True)
    ap.add_argument("--panel", required=True)
    ap.add_argument("--chrom", default="Chr1")
    ap.add_argument("--table", required=True)
    a = ap.parse_args()

    if a.hapfire_ecotype:
        e = pd.read_csv(a.hapfire_ecotype, sep="\t", header=None, names=["founder", "h"],
                        dtype={"founder": str})
        founders = e["founder"].values
        h = e["h"].astype(float).values
    else:
        z = np.load(a.h, allow_pickle=True)
        founders = np.array([str(x) for x in z["founders"]])
        h = np.asarray(z[a.chrom], dtype=float) if a.chrom in z.files else np.asarray(z[z.files[-1]], dtype=float)
    h = h / h.sum() if h.sum() > 0 else h

    w = pd.read_csv(a.weights, sep="\t", dtype={"founder": str})
    truth = dict(zip(w["founder"], w["weight"].astype(float)))
    w_true = np.array([truth.get(f, 0.0) for f in founders])
    w_true = w_true / w_true.sum() if w_true.sum() > 0 else w_true
    is_true = w_true > 0
    true_n = int(is_true.sum())

    eff_n = float(1.0 / np.sum(h ** 2)) if np.sum(h ** 2) > 0 else np.nan
    thr = 0.5 / max(true_n, 1)
    det = h >= thr
    n_detected = int(det.sum())
    recall = float((det & is_true).sum() / true_n) if true_n else np.nan
    precision = float((det & is_true).sum() / n_detected) if n_detected else np.nan
    mass_on_true = float(h[is_true].sum())
    topn = set(np.argsort(h)[::-1][:true_n].tolist())
    trueset = set(np.where(is_true)[0].tolist())
    jac = len(topn & trueset) / len(topn | trueset) if (topn | trueset) else np.nan
    h_rmse = float(np.sqrt(np.mean((h - w_true) ** 2)))

    row = dict(tool=a.tool, pool=a.pool, panel=a.panel, true_n=true_n,
               seed=int(a.pool.split("_s")[1].split("_")[0]) if "_s" in a.pool else -1,
               eff_n=round(eff_n, 3), n_detected=n_detected, recall=round(recall, 3),
               precision=round(precision, 3), mass_on_true=round(mass_on_true, 4),
               top_n_jaccard=round(jac, 3), h_rmse=round(h_rmse, 5))
    print(f"[ecotype] {a.pool}: true_n={true_n} eff_n={eff_n:.2f} detected={n_detected} "
          f"recall={recall:.2f} prec={precision:.2f} mass_true={mass_on_true:.3f} jac={jac:.2f}",
          file=sys.stderr)
    Path(a.table).parent.mkdir(parents=True, exist_ok=True)
    cols = ["tool", "pool", "panel", "true_n", "seed", "eff_n", "n_detected", "recall",
            "precision", "mass_on_true", "top_n_jaccard", "h_rmse"]
    pd.DataFrame([row])[cols].to_csv(a.table, sep="\t", index=False)


if __name__ == "__main__":
    main()
