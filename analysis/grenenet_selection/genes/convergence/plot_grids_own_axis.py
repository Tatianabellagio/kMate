#!/usr/bin/env python
"""Re-render the per-garden trajectory grid for screen candidates, ordering the gardens by
the climate axis the GEA actually found each candidate on (`gea_best_axis`), instead of
the bio1/bio12 choice made by `screen_3criteria.py`.

Why: 3 of the 89 candidates were found on bio1 and none on bio12. Ordering a bio6 or bio15
hit by precipitation scatters its gradient across the page and the figure cannot be used
to judge the climate criterion (see screen_c3_all_axes.py). The locus/LD figure does not
depend on the axis and is not re-rendered.

Output -> results/plots/screen_ownaxis/<sym>_garden_trajectories.{png,pdf}. Kept separate
from plots/screen/ so the two orderings are not confused.

Usage:
    python plot_grids_own_axis.py --tiers A B C D      # from screen_visual_review.csv
    python plot_grids_own_axis.py --genes AT5G10130 AT2G27030
Env: kmate. Compute node. ~30-60 s per gene (reads the three gen matrices).
"""
from __future__ import annotations
import os
import sys
import json
import shutil
import argparse
import subprocess
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(HERE, "..", "..")))
import lib                                                       # noqa: E402

PY = sys.executable
GRID = f"{lib.GEA}/genes/dissection/plot_variant_garden_grid.py"
SRC = f"{lib.GEA}/genes/dissection/results/plots"
OUT = f"{HERE}/results"
DST = f"{OUT}/plots/screen_ownaxis"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tiers", nargs="*", default=None)
    ap.add_argument("--genes", nargs="*", default=None)
    a = ap.parse_args()
    os.makedirs(DST, exist_ok=True)

    T = pd.read_csv(f"{OUT}/screen_visual_review.csv", dtype=str).fillna("")
    M = pd.read_csv(f"{OUT}/master_candidate_genes.csv").set_index("gene")
    if a.genes:
        T = T[T.target_gene.isin(a.genes)]
    elif a.tiers:
        T = T[T.final_tier.str[0].isin(a.tiers)]
    print(f"{len(T)} grids -> {DST}\n", flush=True)

    ok, bad = [], []
    for r in T.itertuples():
        axis = M.loc[r.target_gene, "gea_best_axis"] if r.target_gene in M.index else "bio1"
        axis = axis if isinstance(axis, str) and axis else "bio1"
        rl, al = int(r.ref_len), int(r.alt_len)
        cfg = {"chrom": r.chrom, "pos": int(r.pos), "ref_len": rl, "alt_len": al,
               "sym": r.sym, "store_row": int(r.store_row), "climate": axis,
               "label": f"{abs(al - rl)} bp {r.region} {'deletion' if rl > al else 'insertion'}"}
        print(f"  {r.sym:12s} {r.target_gene}  axis={axis}", flush=True)
        p = subprocess.run([PY, GRID, json.dumps(cfg)], capture_output=True, text=True)
        if p.returncode != 0:
            bad.append((r.sym, (p.stderr or "").strip().splitlines()[-1:]))
            print(f"      FAILED: {bad[-1][1]}", flush=True)
            continue
        n = 0
        for e in (".png", ".pdf"):
            s = f"{SRC}/{r.sym}_garden_trajectories{e}"
            if os.path.exists(s):
                shutil.copy2(s, f"{DST}/{r.sym}_garden_trajectories{e}")
                n += 1
        (ok if n else bad).append(r.sym if n else (r.sym, "no figure written"))

    print(f"\nrendered {len(ok)}: {', '.join(ok)}")
    if bad:
        print(f"failed {len(bad)}: {bad}")


if __name__ == "__main__":
    main()
