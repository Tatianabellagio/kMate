#!/usr/bin/env python
"""The GEA's pools: each plot at its last available generation, as pseudo-generation 9 (run before _submit_gea.sh).

The GrENE-net phase-1 GEA unit is one pool per plot, taken at the plot's last available generation (1, 2 or 3; plots
go extinct). The list is phase-1's key_files/final_gen.csv (sample_name = site_gen_plot pool id, gen = its
generation); site 33 is dropped (lib.EXCLUDED_SITES). This is a pure merge: each pool's row is copied from kMate's
flower-weighted pool matrix of its own generation, nothing is re-pooled.

Writes in lib.POOL_MATRICES, for --kind snp and --kind nonsnp:
  pool_gen9_{kind}_af.npy    float32 [pools x records]
  pool_gen9_{kind}.meta.csv  the same schema as pool_gen{1,2,3}_{kind}.meta.csv (pool, site, plot, ..., bio1-19)
Environment: kmate.  Usage: assemble_lastgen.py --kind snp; assemble_lastgen.py --kind nonsnp
"""
from __future__ import annotations
import argparse, os, sys
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import lib

FINAL_GEN = "/global/scratch/users/tbellg/gea_grene-net/key_files/final_gen.csv"
PM = lib.POOL_MATRICES          # the panel version's (KMATE_PANEL; arch3: common/results/pool_matrices)
OUT_GEN = 9  # pseudo-generation label for the merged last-gen set


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--kind", default="snp", choices=["snp", "nonsnp"])
    KIND = ap.parse_args().kind
    fg = pd.read_csv(FINAL_GEN)
    fg["pool"] = fg["sample_name"].astype(str)
    fg["gen"] = fg["gen"].astype(int)
    fg = fg[~fg["site"].astype(int).isin(lib.EXCLUDED_SITES)]     # ⛔ lib.EXCLUDED_SITES (site 33)
    # canonical pool order (same key build_pool_matrix uses)
    fg = fg.sort_values("pool", key=lambda s: s.map(
        lambda p: tuple(int(x) for x in p.split("_")))).reset_index(drop=True)
    print(f"final_gen pools: {len(fg)} | by gen:\n{fg.gen.value_counts().sort_index().to_string()}")

    # load each generation's pool meta + memmap, build pool->row maps
    gmeta, gmat, grow = {}, {}, {}
    for g in (1, 2, 3):
        m = pd.read_csv(f"{PM}/pool_gen{g}_{KIND}.meta.csv")
        m["pool"] = m["pool"].astype(str)
        gmeta[g] = m
        gmat[g] = np.load(f"{PM}/pool_gen{g}_{KIND}_af.npy", mmap_mode="r")
        grow[g] = {p: i for i, p in enumerate(m["pool"])}
        print(f"  gen{g}: {gmat[g].shape[0]} pools x {gmat[g].shape[1]:,} records")

    n_rec = gmat[1].shape[1]
    assert all(gmat[g].shape[1] == n_rec for g in (1, 2, 3)), "record axis mismatch across gens"

    # coverage check
    missing = [r.pool for _, r in fg.iterrows() if r.pool not in grow[r.gen]]
    print(f"\nlastgen pools found in kMate: {len(fg) - len(missing)}/{len(fg)}")
    if missing:
        print(f"  MISSING ({len(missing)}): {missing[:20]}")
        fg = fg[~fg.pool.isin(missing)].reset_index(drop=True)

    out_path = f"{PM}/pool_gen{OUT_GEN}_{KIND}_af.npy"
    out = np.lib.format.open_memmap(out_path, mode="w+", dtype=np.float32,
                                    shape=(len(fg), n_rec))
    meta_rows = []
    for i, r in fg.iterrows():
        g = r.gen
        row = grow[g][r.pool]
        out[i, :] = gmat[g][row, :]
        meta_rows.append(gmeta[g].iloc[row])
        if (i + 1) % 50 == 0:
            print(f"  merged {i + 1}/{len(fg)} pools", flush=True)
    out.flush()
    meta = pd.DataFrame(meta_rows).reset_index(drop=True)
    meta.to_csv(f"{PM}/pool_gen{OUT_GEN}_{KIND}.meta.csv", index=False)
    print(f"\n-> {out_path}  [{out.shape[0]} x {out.shape[1]:,}]")
    print(f"-> {PM}/pool_gen{OUT_GEN}_{KIND}.meta.csv  ({len(meta)} pools, "
          f"bio1 range {meta.bio1.min():.1f}-{meta.bio1.max():.1f})")


if __name__ == "__main__":
    main()
