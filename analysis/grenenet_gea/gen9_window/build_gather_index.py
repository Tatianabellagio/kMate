#!/usr/bin/env python
"""Precompute, ONCE, the panel-row gather index that pulls the gen9 class-matrix
records' AF out of a full per-sample genome-wide TSV (8.49M rows, panel order).

The global gen9 class matrices were filtered from the af_store; each record stores
`col` = its index within the snp / nonsnp source array. The af_store snp/nonsnp
arrays are the panel split by snp_mask, so:
    panel_row(record) = where(snp_mask)[0][col]    (snp class)
                      = where(~snp_mask)[0][col]    (sv / smallindel class)
A genome-wide window TSV is in panel order, so af_full[panel_row] = that record's AF.

Outputs (results/grenenet_gea/gen9_window/):
  gather_rows.npy   int64 [N]      panel row for each gen9 record (snp|sv|smallindel order)
  records.csv       chrom,pos,cls  row-aligned to gather_rows (the evolved AF columns)
"""
import os, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import lib

CM = f"{lib.GEA}/phase1_replication/class_matrices"
STORE = lib.AF_STORE
OUT = f"{lib.GEA}/gen9_window"
CLASSES = ["snp", "sv", "smallindel"]
SRC = {"snp": "snp", "sv": "nonsnp", "smallindel": "nonsnp"}

snp_mask = np.load(f"{STORE}/snp_mask.npy")
prows = {"snp": np.where(snp_mask)[0], "nonsnp": np.where(~snp_mask)[0]}

gather, chrom, pos, cls = [], [], [], []
for c in CLASSES:
    recs = pd.read_csv(f"{CM}/{c}_gen9.records.csv", usecols=["chrom", "pos", "col"])
    pr = prows[SRC[c]][recs.col.to_numpy()]
    gather.append(pr); chrom.append(recs.chrom.to_numpy())
    pos.append(recs.pos.to_numpy()); cls += [c] * len(recs)
    print(f"  {c}: {len(recs):,} gen9 records")

gather = np.concatenate(gather)
df = pd.DataFrame(dict(chrom=np.concatenate(chrom), pos=np.concatenate(pos), cls=cls))
os.makedirs(OUT, exist_ok=True)
np.save(f"{OUT}/gather_rows.npy", gather)
df.to_csv(f"{OUT}/records.csv", index=False)
print(f"total gen9 records: {len(gather):,}  -> {OUT}/gather_rows.npy, records.csv")

# --- validate: global tsv gathered == global af_store decoded, for one sample ---
s = str(np.load(f"{STORE}/samples.npy", allow_pickle=True)[0])
af_full = pd.read_csv(f"{lib.PROJ}/results/grenenet_kmate_arch3/{s}.tsv", sep="\t",
                      usecols=["alt_freq"]).alt_freq.to_numpy(np.float64)
# snp-class check
recs_snp = pd.read_csv(f"{CM}/snp_gen9.records.csv", usecols=["col"]).col.to_numpy()
store_snp = lib.decode_af(np.load(f"{STORE}/af_snp/{s}.npy"))[recs_snp]
gathered_snp = af_full[prows["snp"][recs_snp]]
ok = np.allclose(store_snp, gathered_snp, equal_nan=True, atol=2e-4)
print(f"VALIDATION (sample {s}, snp class): gathered==af_store -> {ok} "
      f"(max|diff|={np.nanmax(np.abs(store_snp-gathered_snp)):.2e})")
