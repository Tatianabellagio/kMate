#!/usr/bin/env python3
"""Quick single-pool vg-SV scorer — replicates build_4tool_table.score_sv for vg only.
Usage: _score_one_vgsv.py <pool> [panel=p80] [svlen=50]"""
import sys
from pathlib import Path
import numpy as np, pandas as pd, scipy.sparse as sp

sys.path.insert(0, str(Path(__file__).resolve().parent))
import score_snp_fair as ssf  # metrics()

pool = sys.argv[1]; panel = sys.argv[2] if len(sys.argv) > 2 else "p80"
svlen = int(sys.argv[3]) if len(sys.argv) > 3 else 50
root = Path("/global/scratch/users/tbellg/kmate")
truth = root / f"benchmarks/{panel}/sims/{pool}/recomb_truth.tsv.gz"
meta  = root / f"benchmarks/{panel}/data/var_pa_{panel}.meta.npz"
called= root / f"benchmarks/{panel}/data/var_pa_{panel}.var_called.npz"
vgsv  = root / f"benchmarks/accuracy_vs_competitors/work/vg_sv_{pool}.tsv"

tr = pd.read_csv(truth, sep="\t")
m = np.load(meta, allow_pickle=True)
assert len(tr) == len(m["pos"]), f"{len(tr)} != {len(m['pos'])}"
svmask = np.abs(m["alt_len"].astype(int) - m["ref_len"].astype(int)) >= svlen
vc = sp.load_npz(called).tocsr()
F = vc.shape[0]
n_called = np.asarray(vc[:, svmask].sum(axis=0)).ravel().astype(int)
trV = tr[svmask].reset_index(drop=True).copy()
trV["svidx"] = np.arange(len(trV)); trV["n_called"] = n_called

est = pd.read_csv(vgsv, sep="\t")[["svidx", "est"]]
print(f"pool={pool}  panel SVs={len(trV):,}  F={F}  vg est rows={len(est):,}")
for basis, sub in [("allrec", trV), ("fullcalled", trV[trV["n_called"] == F])]:
    j = sub.merge(est, on="svidx", how="inner")
    mt = ssf.metrics(j["truth_af"].values, j["est"].values)
    print(f"  vg-SV [{basis:10s}] n={mt['n']:>6}  MAE={mt['MAE']:.4f}  RMSE={mt['RMSE']:.4f}  "
          f"R2={mt['R2']:.4f}  pearson_r={mt['pearson_r']:.4f}")
