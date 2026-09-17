#!/usr/bin/env python
"""Render locus/LD figure + own-axis garden grid for every gene in a candidate table.

Generalises `plot_screen_candidates.py` (hard-wired to screen_top_loci.csv and to the
bio1/bio12 axis) so the same two figures can be produced for `screen_top_loci_ownaxis.csv`
and `screen_gwas_rescreen.csv`. Differences:

  * gardens in the grid are ordered by the axis the table says the hit was found on
    (`C3own_axis`, or `C3_axis_used` for the GWAS table), never bio1/bio12 by default
  * genes with no coordinates in the master table (ncRNAs, TE genes) fall back to the
    TAIR10 genes+transposons GFF instead of being skipped
  * a `--worker i --nworkers n` split so several processes can share the node; the two
    dissection scripts write per-symbol files, so workers on different genes do not collide
  * a symbol that is already used by a DIFFERENT gene in the output folder gets the locus id
    appended -- the symbol column is a display label, not an identifier (two genes are
    called DFC, for instance)

Outputs land in results/plots/screen/ as
  <sym>_combined[_gwas].{png,pdf}           locus, founder panel, LD triangle
  <sym>_garden_trajectories_ownaxis.{png,pdf}

Usage:
  python plot_shortlist.py --table results/screen_top_loci_ownaxis.csv --worker 0 --nworkers 3
  python plot_shortlist.py --table results/screen_gwas_rescreen.csv --only-pass
Env: kmate. Compute node. ~60-80 s per gene.
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
LOCUS, GRID = f"{DIS}/plot_locus_combined.py", f"{DIS}/plot_variant_garden_grid.py"
SRC_LOCI, SRC_PLOTS = f"{DIS}/loci", f"{DIS}/results/plots"
OUT = f"{HERE}/results"
DST = f"{OUT}/plots/screen"


def gff_coords(gene: str):
    """(start, end) from the genes+transposons GFF for ids missing from the master table."""
    for path in (lib.TAIR10_GENES_TE, lib.TAIR10_GENES):
        try:
            with open(path) as fh:
                for ln in fh:
                    if f"ID={gene};" in ln and "\tgene\t" in ln or (
                            f"ID={gene};" in ln and "\ttransposable_element_gene\t" in ln):
                        f = ln.split("\t")
                        return int(f[3]), int(f[4])
        except FileNotFoundError:
            continue
    return None


def run(script, cfg, label):
    p = subprocess.run([PY, script, json.dumps(cfg)], capture_output=True, text=True)
    if p.returncode != 0:
        err = (p.stderr or "").strip().splitlines()[-1:] or ["?"]
        print(f"      FAILED ({label}): {err[0][:160]}", flush=True)
        return False
    return True


def collect(src, stem, dst_stem):
    n = 0
    for e in (".png", ".pdf"):
        s = f"{src}/{stem}{e}"
        if os.path.exists(s):
            shutil.copy2(s, f"{DST}/{dst_stem}{e}")
            n += 1
    return n


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--table", required=True)
    ap.add_argument("--worker", type=int, default=0)
    ap.add_argument("--nworkers", type=int, default=1)
    ap.add_argument("--only-pass", action="store_true",
                    help="GWAS table: only rows with PASS_new (else every row)")
    ap.add_argument("--genes", nargs="*", default=None)
    a = ap.parse_args()
    os.makedirs(DST, exist_ok=True)

    T = pd.read_csv(a.table)
    if a.only_pass and "PASS_new" in T.columns:
        T = T[T.PASS_new]
    if a.genes:
        T = T[T.target_gene.isin(a.genes)]
    T = T.reset_index(drop=True)
    T = T.iloc[a.worker::a.nworkers]
    M = pd.read_csv(f"{OUT}/master_candidate_genes.csv").set_index("gene")
    # symbol -> gene map of what is already in the folder, so a re-used symbol is disambiguated
    prev = pd.read_csv(f"{OUT}/screen_visual_review.csv", dtype=str)
    used = dict(zip(prev.sym, prev.target_gene))
    axis_col = "C3own_axis" if "C3own_axis" in T.columns else "C3_axis_used"
    print(f"worker {a.worker}/{a.nworkers}: {len(T)} genes from {os.path.basename(a.table)}",
          flush=True)

    ok, bad = [], []
    for r in T.itertuples():
        g = r.target_gene
        sym = str(r.symbol) if isinstance(r.symbol, str) and r.symbol and r.symbol != "nan" else g
        sym = sym.replace("/", "_").replace(" ", "_")
        if used.get(sym, g) != g:
            sym = f"{sym}_{g}"
        found = str(M.found_by.get(g, "GEA"))
        scan = "gwas" if found == "GWAS" else "gea"
        axis = str(getattr(r, axis_col)) if isinstance(getattr(r, axis_col), str) else "bio1"
        gea_axis = M.gea_best_axis.get(g) if g in M.index else None
        gea_axis = gea_axis if isinstance(gea_axis, str) else axis
        coords = None
        if g in M.index and pd.notna(M.gene_start.get(g)) and pd.notna(M.gene_end.get(g)):
            coords = (int(M.gene_start[g]), int(M.gene_end[g]))
        else:
            coords = gff_coords(g)
        rl, al = int(r.ref_len), int(r.alt_len)
        region = str(r.region) if isinstance(r.region, str) else "intergenic"
        print(f"  {sym:14s} {g}  {r.chrom}:{r.pos}  {r.ftier}  axis={axis}  scan={scan}"
              f"{'' if coords else '  (no gene coords: locus plot skipped)'}", flush=True)
        got = 0
        if coords:
            cfg = {"gene": g, "sym": sym, "chrom": r.chrom, "gstart": coords[0],
                   "gend": coords[1], "vpos": int(r.pos), "ref_len": rl, "alt_len": al,
                   "axis": gea_axis if gea_axis.startswith("bio") else "bio1",
                   "pad": 8000, "scan": scan, "region": region}
            if run(LOCUS, cfg, "locus"):
                stem = f"{sym}_combined{'_gwas' if scan == 'gwas' else ''}"
                got += collect(SRC_LOCI, stem, stem)
        grid = {"chrom": r.chrom, "pos": int(r.pos), "ref_len": rl, "alt_len": al,
                "sym": sym, "store_row": int(r.store_row), "climate": axis,
                "label": f"{abs(al - rl)} bp {region} {'deletion' if rl > al else 'insertion'}"}
        if run(GRID, grid, "grid"):
            got += collect(SRC_PLOTS, f"{sym}_garden_trajectories",
                           f"{sym}_garden_trajectories_ownaxis")
        (ok if got else bad).append(sym if got else (sym, "nothing written"))

    print(f"\nworker {a.worker}: rendered {len(ok)}, failed {len(bad)}: {bad}")


if __name__ == "__main__":
    main()
