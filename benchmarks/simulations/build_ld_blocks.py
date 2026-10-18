#!/usr/bin/env python
"""The r2 >= 0.1 LD blocks of each variant panel on Chr1, used by kMate's `--unit ld` fits (config.VARIANT_PANELS).

kMate's own partition (src/kmate/ld_partition.py: ld_blocks_tsv, window 100 records) on the panel's genotypes.
Usage: python build_ld_blocks.py [PANEL ...] [--out-dir DIR]   (default: every panel, written to its configured path;
--out-dir writes <DIR>/<panel>.tsv instead, to compare with the existing files). An existing output file is kept
(the package caches blocks per panel). Env: kmate. Compute node.
"""
import os, sys
sys.path.insert(0, "/global/scratch/projects/fc_moilab/tbellg/kmate/src")
from kmate.ld_partition import ld_blocks_tsv
import config as C

R2, WINDOW = 0.1, 100

if __name__ == "__main__":
    args = sys.argv[1:]
    out_dir = args[args.index("--out-dir") + 1] if "--out-dir" in args else None
    panels = [a for a in args if not a.startswith("--") and a != out_dir] or list(C.VARIANT_PANELS)
    for p in panels:
        out = f"{out_dir}/{p}.tsv" if out_dir else C.VARIANT_PANELS[p]["ld_blocks"]
        print(ld_blocks_tsv(C.VARIANT_PANELS[p]["prefix"], C.CHROM, R2, out, window=WINDOW), flush=True)
