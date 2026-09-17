#!/usr/bin/env python
"""Render the two inspection figures for candidates that survive `screen_3criteria.py`.

For each candidate gene, one per variant lead:
  <sym>_combined[_gwas].{png,pdf}   dissection/plot_locus_combined.py -- Manhattan of the
                                    scan that FOUND it + founder panel + founder-LD triangle
  <sym>_garden_trajectories.{png,pdf}  dissection/plot_variant_garden_grid.py -- per-garden
                                    allele-frequency trajectory, cold -> warm

Both land in `results/plots/screen/` (kept separate from `plots/loci/`, which holds the
earlier themed-shortlist run, so the two selections do not get mixed).

Scan choice mirrors plot_top_loci.py: a GWAS-only candidate must be drawn against the GWAS
panel or its Manhattan is empty and reads as a null.

Usage:
    python plot_screen_candidates.py --n 12            # top N by screen rank
    python plot_screen_candidates.py --genes AT5G04170 AT4G23260
Env: kmate. Compute node.
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
DIS = f"{lib.GEA}/genes/dissection"
LOCUS = f"{DIS}/plot_locus_combined.py"
GRID = f"{DIS}/plot_variant_garden_grid.py"
# plot_locus_combined.py writes to dissection/loci/ (NOT results/loci/, which holds the
# older curated GEA set + CANDIDATE_VERDICTS.md); plot_variant_garden_grid.py writes to
# dissection/results/plots/. Verified against each script's own output path.
SRC_LOCI, SRC_PLOTS = f"{DIS}/loci", f"{DIS}/results/plots"
OUT = f"{HERE}/results"
DST = f"{OUT}/plots/screen"


def run(script: str, cfg: dict, label: str) -> bool:
    p = subprocess.run([PY, script, json.dumps(cfg)], capture_output=True, text=True)
    tail = (p.stdout or "").strip().splitlines()
    for ln in tail[-6:]:
        print(f"      {ln}")
    if p.returncode != 0:
        print(f"      FAILED ({label}): {(p.stderr or '').strip().splitlines()[-1:]}")
        return False
    return True


def collect(src_dir: str, stem: str) -> int:
    n = 0
    for e in (".png", ".pdf"):
        s = f"{src_dir}/{stem}{e}"
        if os.path.exists(s):
            shutil.copy2(s, f"{DST}/{stem}{e}")
            n += 1
    return n


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=12)
    ap.add_argument("--genes", nargs="*", default=None)
    ap.add_argument("--min-pass", type=int, default=3,
                    help="only candidates meeting at least this many of the 3 criteria")
    a = ap.parse_args()
    os.makedirs(DST, exist_ok=True)

    # Drive off screen_top_loci.csv -- the SAME row the ranking table quotes. Re-deriving a
    # representative here (sort + drop_duplicates over screen_3criteria) picked a DIFFERENT
    # af_store record for 5 of 18 genes, because ambiguous keys expand to several records per
    # gene. Two of those five had opposite-sign climate correlations (AT5G40855 -0.62 vs
    # +0.29; AT2G05840 -0.71 vs +0.37), so the figure contradicted the table it was filed
    # under. One table, one row, one figure.
    S = pd.read_csv(f"{OUT}/screen_top_loci.csv")
    M = pd.read_csv(f"{OUT}/master_candidate_genes.csv").set_index("gene")

    sel = S[S.target_gene.isin(a.genes)] if a.genes else S
    sel = sel.head(a.n)
    print(f"{len(sel)} candidate genes to render -> {DST}\n")

    ok, bad = [], []
    for r in sel.itertuples():
        g = r.target_gene
        if g not in M.index:
            bad.append((g, "not in master table")); continue
        m = M.loc[g]
        # Some master-table rows have no gene coordinates (NaN gene_start/gene_end). Skip
        # those instead of dying: an uncaught ValueError here aborted a 14-gene batch on
        # gene 7 and silently lost the remaining 8.
        if not (pd.notna(m.gene_start) and pd.notna(m.gene_end)):
            bad.append((g, "no gene coordinates in master table")); continue
        sym = str(m.symbol) if isinstance(m.symbol, str) and m.symbol else g
        sym = sym.replace("/", "_")
        scan = "gwas" if str(m.found_by) == "GWAS" else "gea"
        axis = m.gea_best_axis if isinstance(m.gea_best_axis, str) and m.gea_best_axis else "bio1"
        print(f"  {sym} ({g})  {r.chrom}:{r.pos}  {r.ftier}  scan={scan}  "
              f"C3 {r.C3_best_axis} rho={r.C3_best_rho:+.2f} p={r.C3_best_p:.1e}")

        cfg = {"gene": g, "sym": sym, "chrom": r.chrom,
               "gstart": int(m.gene_start), "gend": int(m.gene_end),
               "vpos": int(r.pos), "ref_len": int(r.ref_len), "alt_len": int(r.alt_len),
               "axis": axis, "pad": 8000, "scan": scan,
               "region": str(r.region) if isinstance(r.region, str) else "intergenic"}
        got = 0
        if run(LOCUS, cfg, "locus"):
            got += collect(SRC_LOCI, f"{sym}_combined{'_gwas' if scan == 'gwas' else ''}")

        # store_row pins the EXACT af_store record: the four-field key is not unique
        # (14.1% of non-SNP records collide), so without it the grid can plot a different
        # variant than the screen row it came from.
        grid = {"chrom": r.chrom, "pos": int(r.pos), "ref_len": int(r.ref_len),
                "alt_len": int(r.alt_len), "sym": sym, "store_row": int(r.store_row),
                "climate": str(r.C3_best_axis),
                "label": f"{abs(int(r.alt_len) - int(r.ref_len))} bp {r.region} "
                         f"{'deletion' if r.ref_len > r.alt_len else 'insertion'}"}
        if run(GRID, grid, "grid"):
            got += collect(SRC_PLOTS, f"{sym}_garden_trajectories")
        (ok if got else bad).append((g, sym) if got else (g, "no figure written"))

    print(f"\nrendered {len(ok)}: {', '.join(s for _, s in ok)}")
    if bad:
        print(f"failed {len(bad)}: {bad}")


if __name__ == "__main__":
    main()
