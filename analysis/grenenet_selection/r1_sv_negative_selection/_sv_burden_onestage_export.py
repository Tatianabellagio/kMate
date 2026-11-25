#!/usr/bin/env python
"""Export the founder x garden data for the ONE-STAGE burden model (`_sv_burden_onestage.R`).

Long table, one row per founder x garden (231 x 31 = 7,161): founder log-fitness y (as `_sv_burden_lmm.py`),
garden id and bio1, and the founder's standardised burden rates for SNPs, insertions and deletions in
coding, reference TE, promoter and intergenic sequence (as `_sv_burden_full_model.py`, CONTEXT=CDS,promoter),
plus the assembly route. TE added 2026-09-22, prompted by a TE enrichment among the GEA hits
(`gea_hit_landscape.ipynb`) -- was tested for SFS/Tajima's D shift (Fig. 1 / S1) but not in this burden
model. Also the SNP kinship matrix (founder_kinship.npy) with founder ids.
Env: kmate, SV_UNIT=event. Writes results/sv_adaptive/onestage/{long.csv,kinship.csv}.
"""
import os, sys
import numpy as np, pandas as pd
os.environ.setdefault("SV_UNIT", "event")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _sv_burden_lmm as BL

# BURDEN_ORIGIN selects which origin class the burden rates are built from (see _sv_burden_build.origin_masks):
#   shared (default, what the reported models use) = carried by >=1 cold-origin AND >=1 warm-origin founder
#   all = no origin filter, i.e. every variant passing call-rate + MAC. Writes to onestage_<origin>/.
ORIGIN = os.environ.get("BURDEN_ORIGIN", "shared")
assert ORIGIN in ("shared", "all", "cold", "warm")
OUT = f"{BL.OUT}/onestage" + ("" if ORIGIN == "shared" else f"_{ORIGIN}")
TERMS = [f"{k}_{x}" for k in ("snp", "ins", "del") for x in ("CDS", "TE", "promoter", "intergenic")]


def main():
    os.makedirs(OUT, exist_ok=True)
    D = pd.read_csv(f"{BL.OUT}/founder_burden{BL.TAG}.csv", dtype={"founder": str})
    Zf = np.load(BL.EF, allow_pickle=True); assert np.array_equal(Zf["founders"].astype(str), D.founder.to_numpy())
    Y = np.log((Zf["H0"][None, :] + Zf["DH_rel"] + BL.EPS) / (Zf["H0"][None, :] + BL.EPS))       # gardens x founders
    F = pd.DataFrame({"founder": D.founder, "cactus": D.is_cactus.astype(int)})
    for t in TERMS:
        k, x = t.split("_"); r = (D[f"{k}_{ORIGIN}_{x}_all"] / D[f"called_{k}_{ORIGIN}_{x}_all"]).to_numpy(float)
        F[t] = (r - r.mean()) / r.std()
    rows = []
    for s in range(Y.shape[0]):
        rows.append(F.assign(garden=int(Zf["sites"][s]), bio1=float(Zf["bio1"][s]), y=Y[s]))
    L = pd.concat(rows, ignore_index=True)
    L.to_csv(f"{OUT}/long.csv", index=False)
    K = np.load(f"{BL.OUT}/founder_kinship.npy")
    pd.DataFrame(K, index=D.founder, columns=D.founder).to_csv(f"{OUT}/kinship.csv")
    print(f"origin={ORIGIN} -> {OUT}")
    print(f"long table {L.shape}, {L.founder.nunique()} founders x {L.garden.nunique()} gardens; kinship {K.shape}")


if __name__ == "__main__":
    main()
