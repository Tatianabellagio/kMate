#!/usr/bin/env python
"""Which af_store record did the GEA actually flag, where same-length alleles collide?

`screen_3criteria.resolve_records` expands every (chrom, pos, ref_len, alt_len) key that
matches more than one af_store row (1,111 of 2,529 screen records) and the shortlists then
keep the colliding record that best passes the movement / climate criteria -- i.e. the
allele is chosen on the outcome being judged. The GEA input (`wza_in_clq09_tile`) does
distinguish the colliders: its `MAF` column is exactly the gen9 pool-mean minor allele
frequency (checked on 400 unambiguous records, max |diff| 3e-8). So each screen record is
flagged by whether its own gen9 MAF matches a Bonferroni-significant GEA record at its key.

Output -> results/screen_record_gea_match.csv
  store_row, target_gene, key_ambiguous, maf_g9, gea_sig_record (True / False;
  NaN for records whose key is not in the GEA pool, i.e. GWAS-only)
env kmate, compute node, ~3 min.
"""
import os, sys
import numpy as np
import pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(HERE, "..", "..")))
import lib                                                      # noqa: E402

K = ["chrom", "pos", "ref_len", "alt_len"]
S = pd.read_csv(f"{HERE}/results/screen_3criteria.csv",
                usecols=K + ["store_row", "target_gene", "key_ambiguous"])
P = pd.read_csv(f"{HERE}/results/gea_pool_records.csv")[K + ["MAF"]].drop_duplicates()

rows = np.unique(S.store_row.to_numpy())
M = np.load(f"{lib.GEA}/common/results/pool_matrices/pool_gen9_nonsnp_af.npy", mmap_mode="r")
A = np.stack([np.asarray(M[i])[rows] for i in range(M.shape[0])])
mu = np.nanmean(A, 0)
maf = pd.Series(np.minimum(mu, 1 - mu), index=rows)
S["maf_g9"] = S.store_row.map(maf)

inpool = S.set_index(K).index.isin(P.set_index(K).index)
m = S.merge(P, on=K, how="left")
m["hit"] = (m.MAF - m.maf_g9).abs() < 1e-6
flag = m.groupby("store_row").hit.any()
S["gea_sig_record"] = S.store_row.map(flag).where(inpool)
S.to_csv(f"{HERE}/results/screen_record_gea_match.csv", index=False)
a = S[S.key_ambiguous]
print(f"screen records {len(S)}; ambiguous {len(a)}: GEA-significant "
      f"{int((a.gea_sig_record == True).sum())}, not {int((a.gea_sig_record == False).sum())}, "
      f"not in pool {int(a.gea_sig_record.isna().sum())}; unambiguous in pool but unmatched "
      f"{int(((~S.key_ambiguous) & (S.gea_sig_record == False)).sum())}")
