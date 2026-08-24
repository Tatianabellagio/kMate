#!/usr/bin/env python
"""Extract one gen9 member sample's WINDOW-mode AF at the gen9 class records.
Reads its genome-wide window TSV, gathers at the precomputed panel rows, saves a
float32 vector (NaN preserved) row-aligned to gen9_window/records.csv."""
import os, sys
import numpy as np, pandas as pd

ROOT = "/global/scratch/users/tbellg/kmate"
WIN = f"{ROOT}/results/grenenet_kmate_window"
GW = f"{ROOT}/results/grenenet_gea/gen9_window"
sample = sys.argv[1]
out = f"{GW}/persample/{sample}.npy"
if os.path.exists(out):
    print(f"{sample}: exists, skip"); sys.exit(0)
gather = np.load(f"{GW}/gather_rows.npy")
af = pd.read_csv(f"{WIN}/{sample}.tsv", sep="\t", usecols=["alt_freq"]).alt_freq.to_numpy(np.float32)
if af.shape[0] != 8_489_646:
    sys.exit(f"{sample}: {af.shape[0]} rows != 8.49M panel")
os.makedirs(f"{GW}/persample", exist_ok=True)
tmp = out + ".tmp.npy"
np.save(tmp, af[gather])
os.replace(tmp, out)
print(f"{sample}: wrote {gather.size:,} window-AF values")
