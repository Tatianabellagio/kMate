#!/usr/bin/env python
"""Aggregate rc_site_rank_count summaries into one table: per test x class x axis
lambda_GC and #Bonferroni/#FDR clq0.9 blocks."""
import glob, json, os
import pandas as pd

OUTD = "/global/scratch/users/tbellg/kmate/results/grenenet_gea/gea_newpanel/rank_count"
rows = []
for f in sorted(glob.glob(f"{OUTD}/summary_*.json")):
    d = json.load(open(f))
    for test in ("kendall", "binomial"):
        t = d[test]
        rows.append(dict(test=test, cls=d["cls"], axis=d["axis"],
                         n_tested=t["n_tested"], lambda_gc=t["lambda_gc"],
                         min_p=t["min_p"], n_var_bonf=t["n_var_bonf"],
                         n_blk_bonf=t["n_blocks_bonf"], n_var_fdr=t["n_var_fdr"],
                         n_blk_fdr=t["n_blocks_fdr"]))
df = pd.DataFrame(rows).sort_values(["test", "cls", "axis"]).reset_index(drop=True)
df.to_csv(f"{OUTD}/rc_summary_all.csv", index=False)
pd.set_option("display.width", 200, "display.max_rows", 200)
print(df.to_string(index=False))
print(f"\n-> {OUTD}/rc_summary_all.csv")
