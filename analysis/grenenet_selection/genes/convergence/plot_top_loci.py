#!/usr/bin/env python
"""Render the combined locus figure for the top convergence candidates.

Thin driver: it does not plot anything itself. It selects candidates from the master
table, assembles the JSON config each one needs, and calls
`dissection/plot_locus_combined.py` -- the existing 3-panel figure (LFMM Manhattan
coloured by climate delta-p + founder haplotype panel ordered cold->hot + founder-LD
triangle, sharing a genomic x-axis), plus its printed LD-confirm of the lead variant
against the gene's own variants. That LD number is the point: it says whether the
candidate variant is genuinely tied to the gene or CARK-style detached.

Selection, in order: the genes hit by BOTH scans, then the strongest GEA-only genes
carrying a theme. `--genes` overrides with an explicit list.

The plotter writes to `dissection/loci/`; figures are copied into this section's
`results/plots/loci/` so the section's outputs stay under its own results/ per LAYOUT.

env: kmate.  Compute node -- needs bcftools and the founder panel VCFs.
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
import lib                                                    # noqa: E402

OUT = f"{HERE}/results"
PLOTTER = os.path.join(os.path.dirname(HERE), "dissection", "plot_locus_combined.py")
SRC_LOCI = os.path.join(os.path.dirname(HERE), "dissection", "loci")
DST_LOCI = f"{OUT}/plots/loci"
PY = "/global/home/users/tbellg/miniforge3/envs/kmate/bin/python"


def pick(M: pd.DataFrame, n: int) -> pd.DataFrame:
    cross = M[M.L_cross.fillna(False)]
    themed = M[(M.n_independent >= 1) & M.themes.notna() & M.themes.ne("")
               & ~M.gene.isin(cross.gene)]
    themed = themed.sort_values(["gea_n_clusters", "gwas_n_gardens"],
                                ascending=False)
    return pd.concat([cross, themed]).head(n)


def config_for(r, V: pd.DataFrame) -> dict | None:
    """Build the plotter's JSON config; the lead variant supplies ref_len/alt_len."""
    ch, pos = r.lead_variant.split(":")
    pos = int(pos)
    v = V[(V.chrom == ch) & (V.pos == pos)]
    if not len(v):
        print(f"  {r.gene}: lead variant {r.lead_variant} not in the variant table")
        return None
    v = v.iloc[0]
    axis = r.gea_best_axis if isinstance(r.gea_best_axis, str) and r.gea_best_axis \
        else "bio1"
    sym = r.symbol if isinstance(r.symbol, str) and r.symbol else r.gene
    return {"gene": r.gene, "sym": sym.replace("/", "_"), "chrom": ch,
            "gstart": int(r.gene_start), "gend": int(r.gene_end),
            "vpos": pos, "ref_len": int(v.ref_len), "alt_len": int(v.alt_len),
            "axis": axis, "pad": 8000,
            # the plotter annotates the lead with its sub-genic region
            "region": str(v.region) if isinstance(v.region, str) else "intergenic"}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=6)
    ap.add_argument("--genes", nargs="*", default=None)
    a = ap.parse_args()

    M = pd.read_csv(f"{OUT}/master_candidate_genes.csv")
    V = pd.read_csv(f"{OUT}/variants_classified.csv")
    sel = M[M.gene.isin(a.genes)] if a.genes else pick(M, a.n)
    os.makedirs(DST_LOCI, exist_ok=True)

    ok, bad = [], []
    for r in sel.itertuples():
        cfg = config_for(r, V)
        if cfg is None:
            bad.append((r.gene, "no lead variant")); continue
        print(f"\n===== {cfg['sym']} ({r.gene}) {cfg['chrom']}:{cfg['vpos']} "
              f"axis={cfg['axis']} =====", flush=True)
        p = subprocess.run([PY, PLOTTER, json.dumps(cfg)],
                           capture_output=True, text=True)
        print(p.stdout[-1500:] if p.stdout else "")
        if p.returncode != 0:
            print(p.stderr[-1200:])
            bad.append((r.gene, "plotter failed")); continue
        moved = False
        for ext in (".png", ".pdf"):
            s = f"{SRC_LOCI}/{cfg['sym']}_combined{ext}"
            if os.path.exists(s):
                shutil.copy2(s, f"{DST_LOCI}/{cfg['sym']}_combined{ext}")
                moved = True
        ok.append(cfg["sym"]) if moved else bad.append((r.gene, "no figure written"))

    print(f"\nrendered {len(ok)}: {', '.join(ok)}")
    if bad:
        print(f"failed {len(bad)}: {bad}")
    print(f"figures -> {DST_LOCI}")


if __name__ == "__main__":
    main()
