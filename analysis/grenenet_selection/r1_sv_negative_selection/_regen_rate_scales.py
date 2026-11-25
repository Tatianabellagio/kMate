#!/usr/bin/env python
"""Regenerate the per-variant scale files (`pooled_rate_scale*.csv`) for the pooled burden arms.

Why this exists. Until 2026-09-29 `_sv_burden_pooled_models.R` wrote that file under a name tagged only
with SUBSET, so a stratified run (which redefines the contexts, sv_CDS -> sv_CDS_rare etc.) silently
clobbered the previous arm's copy. The script now tags the filename with the full scale-relevant
specification; this rebuilds the arms whose files were lost before that fix.

The quantities are deterministic functions of the burden tables -- sd of the carried/called rate across the
founders in scope, and the medians of the called and carried counts -- so nothing needs refitting. The
script VALIDATES the method by reproducing the two files that survived, and refuses to write if they do not
match to machine precision. `plotting` env.
"""
import os, sys
import numpy as np, pandas as pd

G = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results", "sv_adaptive")
OUT = f"{G}/final_models_pooled_all"
CTX = ["CDS", "TE", "promoter", "intergenic"]

B = pd.read_csv(f"{G}/founder_burden_event.csv", dtype={"founder": str})
CAC = (B.is_cactus.astype(str) == "True").to_numpy()
SUBSETS = {"": np.ones(len(B), bool), "cactus": CAC, "shortread": ~CAC}


def scale(rows, keep):
    out = []
    for name, car, cal in rows:
        r = (car[keep] / cal[keep]).to_numpy()
        out.append(dict(context=name, sd_rate=r.std(ddof=1),
                        med_called=float(np.median(cal[keep])), med_carried=float(np.median(car[keep]))))
    return pd.DataFrame(out)


def rows_for(stratify, sctx, bins):
    """Context rows exactly as the R script builds them: stratified contexts come from the strata table,
    every other context stays pooled over all frequencies/sizes from the main burden table."""
    rows = []
    FR = None
    if stratify:
        FR = pd.read_csv(f"{G}/founder_burden_event_{stratify}.csv", dtype={"founder": str})
        assert (FR.founder.to_numpy() == B.founder.to_numpy()).all(), "founder order differs"
    for c in CTX:
        if stratify and c in sctx:
            for b in bins:
                rows.append((f"{c}_{b}",
                             FR[f"ins_{c}_{b}"] + FR[f"del_{c}_{b}"],
                             FR[f"called_ins_{c}_{b}"] + FR[f"called_del_{c}_{b}"]))
        else:
            rows.append((c, B[f"ins_all_{c}_all"] + B[f"del_all_{c}_all"],
                         B[f"called_ins_all_{c}_all"] + B[f"called_del_all_{c}_all"]))
    return rows


# (filename suffix, subset, stratify, stratified contexts, bins)
ARMS = [
    # the plain filename is what an unstratified, whole-panel run writes -- restored because the model
    # notebook reads it, and because EXTRA_FIXED does not change the burdens so it serves REDUCED and FULL alike
    ("",                                                       "",          None,   [],      []),
    ("_x_obio1cobio1cbio1c",                                   "",          None,   [],      []),
    ("_x_obio1cobio1cbio1c_cactus",                            "cactus",    None,   [],      []),
    ("_x_obio1cobio1cbio1c_shortread",                         "shortread", None,   [],      []),
    ("_x_obio1cobio1cbio1ccallratesv_freqCDSrarecommon",       "",          "freq", ["CDS"], ["rare", "common"]),
    ("_x_obio1cobio1cbio1ccallratesv_cactus_freqCDSprivaterarecommon",
                                                               "cactus",    "freq", ["CDS"], ["private", "rare", "common"]),
    ("_x_obio1cobio1cbio1c_freqCDSsmallmediumlarge",           "",          "size", ["CDS"], ["small", "medium", "large"]),
]
# files that survived the clobbering, used to prove the method before anything is written
VALIDATE = [("pooled_rate_scale_cactus.csv",    "cactus", "freq", ["CDS"], ["private", "rare", "common"]),
            ("pooled_rate_scale_shortread.csv", "shortread", None, [], [])]

ok = True
for fn, sub, st, sctx, bins in VALIDATE:
    p = f"{OUT}/{fn}"
    if not os.path.exists(p):
        print(f"  [skip] {fn} absent"); continue
    got = scale(rows_for(st, sctx, bins), SUBSETS[sub]).set_index("context")
    exp = pd.read_csv(p).set_index("context")
    m = np.allclose(got.loc[exp.index, ["sd_rate", "med_called", "med_carried"]].to_numpy(),
                    exp[["sd_rate", "med_called", "med_carried"]].to_numpy(), rtol=0, atol=1e-12)
    print(f"  validate {fn:<44} reproduced exactly: {m}")
    ok &= m
if not ok:
    sys.exit("VALIDATION FAILED -- not writing anything")

print("\nwriting:")
for suf, sub, st, sctx, bins in ARMS:
    df = scale(rows_for(st, sctx, bins), SUBSETS[sub])
    f = f"{OUT}/pooled_rate_scale{suf}.csv"
    df.to_csv(f, index=False)
    print(f"  {os.path.basename(f):<66} {', '.join(df.context)}")
