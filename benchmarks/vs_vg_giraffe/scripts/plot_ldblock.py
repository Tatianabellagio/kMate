#!/usr/bin/env python3
"""Truth-vs-estimate scatters for an LD-block window kMate run:
panel 1 = ALL records; panel 2 = records in LOCAL-FIT blocks only (status==0)."""
import argparse, sys
from pathlib import Path
import numpy as np, pandas as pd
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))  # benchmarks/
from _accuracy_panel import grid_panel

ap = argparse.ArgumentParser()
ap.add_argument("--est", required=True)
ap.add_argument("--truth", required=True)
ap.add_argument("--pack", required=True)
ap.add_argument("--chrom", default="Chr1")
ap.add_argument("--out", required=True)
ap.add_argument("--snps-only", action="store_true")
a = ap.parse_args()

est = pd.read_csv(a.est, sep="\t")
tr = pd.read_csv(a.truth, sep="\t")
assert len(est) == len(tr), (len(est), len(tr))
t = tr["truth_af"].values.astype(np.float32)
e = est["alt_freq"].values.astype(np.float32)
pos = est["pos"].values

pk = np.load(a.pack, allow_pickle=True)
status = pk[f"{a.chrom}_status"]
bstart = pk[f"{a.chrom}_block_start"].astype(np.int64)
bend = pk[f"{a.chrom}_block_end"].astype(np.int64)
order = np.argsort(bstart); bstart, bend, status = bstart[order], bend[order], status[order]

# assign each record to a block via interval search (disjoint sorted blocks)
j = np.searchsorted(bstart, pos, side="right") - 1
in_block = (j >= 0) & (j < len(bstart))
j_clip = np.clip(j, 0, len(bstart) - 1)
in_block &= pos <= bend[j_clip]
rec_status = np.full(len(pos), -1, dtype=np.int64)
rec_status[in_block] = status[j_clip[in_block]]

n_loc = int((status == 0).sum()); n_fb = int((status == 1).sum())
print(f"blocks: {n_loc} local-fit (status0), {n_fb} fallback (status1), {len(status)} total",
      file=sys.stderr)
print(f"records: {(rec_status==0).sum():,} in local-fit blocks, "
      f"{(rec_status==1).sum():,} in fallback, {(rec_status<0).sum():,} outside", file=sys.stderr)

if a.snps_only:
    snp = (est.ref_len == 1) & (est.alt_len == 1)
    t, e, rec_status = t[snp.values], e[snp.values], rec_status[snp.values]

mask_local = rec_status == 0
cells = [
    (f"LD-block window — ALL records", t, e),
    (f"LD-block window — LOCAL-FIT blocks only", t[mask_local], e[mask_local]),
]
grid_panel(cells, "", a.out, ncols=2, cell=4.2)
