#!/usr/bin/env python
"""Write the content x location burdens in the column convention `_sv_burden_pooled_models.R` reads.

`_sv_burden_cargo_build.py` counts, per founder, carried / called SV events by content (TE / nonTE) x
location (CDS / UTR / intron / promoter / intergenic, reference-TE layer removed), per direction. The
pooled model reads `{ins,del}_<origin>_<context>_all` and `called_{ins,del}_<origin>_<context>_all` and sums
insertions and deletions itself, so this renames the classes to contexts the model's POOL_SET=cargo10 names:

    TEcds TEutr TEintron TEpromoter TEintergenic  nTEcds nTEutr nTEintron nTEpromoter nTEintergenic

Founder-level covariates the model also reads from the burden file (bio1, callrate_sv, is_cactus) are copied
from the reported `founder_burden_event.csv`, and a check confirms both tables are built from the SAME event
set: per founder, total carried and total called over all classes equal the reported table's totals over its
six contexts, for both origins.

Env: basic or kmate. Writes results/sv_adaptive/founder_burden_event_cargo.csv.
"""
import os
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
SA = f"{HERE}/results/sv_adaptive"
LOC = {"CDS": "cds", "UTR": "utr", "intron": "intron", "promoter": "promoter", "intergenic": "intergenic"}
CONT = {"TE": "TE", "nonTE": "nTE"}
CTX6 = ["CDS", "UTR", "intron", "TE", "promoter", "intergenic"]


def main():
    D = pd.read_csv(f"{SA}/founder_burden_cargo.csv", dtype={"founder": str})
    R = pd.read_csv(f"{SA}/founder_burden_event.csv", dtype={"founder": str})
    assert D.founder.tolist() == R.founder.tolist()
    O = R[["founder", "bio1", "callrate_sv", "is_cactus"]].copy()
    for o in ("all", "shared"):
        for k in ("ins", "del"):
            tot_c = tot_q = 0
            for c, cn in CONT.items():
                for l, ln in LOC.items():
                    ctx = f"{cn}{ln}"
                    O[f"{k}_{o}_{ctx}_all"] = D[f"{k}_{c}_{l}_{o}"]
                    O[f"called_{k}_{o}_{ctx}_all"] = D[f"called_{k}_{c}_{l}_{o}"]
                    tot_c = tot_c + D[f"{k}_{c}_{l}_{o}"]; tot_q = tot_q + D[f"called_{k}_{c}_{l}_{o}"]
            ref_c = sum(R[f"{k}_{o}_{x}_all"] for x in CTX6)
            ref_q = sum(R[f"called_{k}_{o}_{x}_all"] for x in CTX6)
            assert np.allclose(tot_c, ref_c) and np.allclose(tot_q, ref_q), \
                f"{k}/{o}: cargo classes do not add up to the reported six contexts"
            print(f"  {k}/{o}: totals match the reported burden table for all {len(D)} founders")
    O.to_csv(f"{SA}/founder_burden_event_cargo.csv", index=False)
    print(f"wrote {SA}/founder_burden_event_cargo.csv ({O.shape[1]} columns)")


if __name__ == "__main__":
    main()
